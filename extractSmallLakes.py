"""Extract small lakes (area < ``MAX_LAKE_AREA_KM2``) from every NDWI
raster produced by ``cloudMaskNDWI.py``, rejecting shadowed slopes
using the Copernicus DEM mosaic in ``Data/DEM/``.

For each ``DD-MM-YYYY_<TILE>_NDWI.tif`` in ``Data/Sentinel/NDWI/`` the
script:

1. Thresholds the raster with ``NDWI > NDWI_THRESHOLD`` to obtain a
   binary water mask (NaN pixels - cloud, cirrus, nodata - are treated
   as non-water and simply excluded).
2. Labels 8-connected components of that mask
   (:func:`scipy.ndimage.label`).
3. Reprojects the merged Copernicus DEM onto the scene's UTM grid,
   computes slope in degrees from central differences, and for every
   labelled component computes the mean slope over its pixels.
4. Keeps only components whose contiguous area falls in
   ``[MIN_LAKE_AREA_KM2, MAX_LAKE_AREA_KM2]`` (default 0.01 to 2 km^2)
   AND whose mean slope is below ``SLOPE_MAX_DEG`` (default 10 deg).
   The slope test rejects the terrain-shadow false positives that
   would otherwise show up on steep, north-facing valley walls.
5. Polygonises the retained components with
   :func:`rasterio.features.shapes` (8-connected) and writes each
   polygon as a feature with ``area_km2``, ``n_pixels``,
   ``mean_slope_deg``, ``date`` and ``tile`` attributes.

Outputs
-------
* ``Data/Lakes/DD-MM-YYYY_<TILE>_smallLakes.gpkg`` - one GeoPackage per
  scene, in the scene's native UTM CRS.
* ``Data/Lakes/all_small_lakes.gpkg`` - all polygons across all scenes
  reprojected to EPSG:4326 for easy overlay / QGIS visualisation.

Requirements::

    pip install rasterio geopandas numpy scipy shapely tqdm
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import shapes
from rasterio.io import MemoryFile
from rasterio.merge import merge as rio_merge
from rasterio.vrt import WarpedVRT
from scipy import ndimage
from shapely.geometry import shape
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NDWI_DIR = Path("../Data/Sentinel/NDWI")
DEM_DIR = Path("../Data/DEM")
OUTPUT_DIR = Path("../Data/Lakes")
COMBINED_OUTPUT = OUTPUT_DIR / "all_small_lakes.gpkg"

# NDWI > NDWI_THRESHOLD is treated as water. 0.20 is McFeeters'
# canonical value; drop it to ~0.05 to also capture partially
# ice-covered or turbid lake pixels (at the cost of more false
# positives from wet vegetation and shadow).
NDWI_THRESHOLD = 0.20

# Retain a polygon only if its area is in [MIN, MAX] km^2. MIN removes
# single-pixel speckle; MAX is the "small lakes only" upper bound.
MIN_LAKE_AREA_KM2 = 0.001   # ~10 pixels at 10 m
MAX_LAKE_AREA_KM2 = 2.00   # user-requested upper bound
PIXEL_AREA_M2 = 100.0      # Sentinel-2 output grid (10 m x 10 m)
MIN_LAKE_PIXELS = int(round(MIN_LAKE_AREA_KM2 * 1_000_000.0 / PIXEL_AREA_M2))
MAX_LAKE_PIXELS = int(round(MAX_LAKE_AREA_KM2 * 1_000_000.0 / PIXEL_AREA_M2))

# Reject polygons whose mean terrain slope (from the Copernicus DEM)
# is >= this value. Removes shadowed valley walls that spuriously look
# water-like in NDWI. Set to a large number (e.g. 90.0) to disable.
SLOPE_MAX_DEG = 20.0

# CRS of the combined output layer. WGS84 is the most portable choice
# for a multi-tile dataset; set to None to keep the scene UTM CRS(s).
COMBINED_CRS = "EPSG:4326"

# NDWI file naming produced by cloudMaskNDWI.py:
#   DD-MM-YYYY_<MGRS-tile>_NDWI.tif  (e.g. 04-07-2025_44SKA_NDWI.tif)
# The optional leading "T" tolerates the alternative Sentinel-2 style
# "T44SKA" if the upstream script is ever changed.
NAME_REGEX = re.compile(
    r"(?P<date>\d{2}-\d{2}-\d{4})_T?(?P<tile>\d{2}[A-Z]{3})_NDWI\.tif$"
)


# ---------------------------------------------------------------------------
# DEM mosaic and slope
# ---------------------------------------------------------------------------
def load_dem_mosaic(dem_dir: Path) -> MemoryFile | None:
    """Merge every ``*.tif`` DEM tile in ``dem_dir`` into an in-memory
    raster. Returns ``None`` if the directory does not exist or is empty."""
    if not dem_dir.exists():
        return None
    tif_paths = sorted(dem_dir.glob("*.tif"))
    if not tif_paths:
        return None

    sources = [rasterio.open(p) for p in tif_paths]
    try:
        mosaic, transform = rio_merge(sources)
        profile = sources[0].profile.copy()
        profile.update(
            driver="GTiff",
            height=mosaic.shape[1],
            width=mosaic.shape[2],
            count=1,
            transform=transform,
        )
    finally:
        for src in sources:
            src.close()

    memfile = MemoryFile()
    with memfile.open(**profile) as dst:
        dst.write(mosaic)
    return memfile


def slope_on_target_grid(
    dem_memfile: MemoryFile,
    target_crs,
    target_transform,
    target_shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Reproject the DEM mosaic onto the target grid (typically Sentinel-2
    10 m UTM) and return ``(slope_deg, missing_mask)``.

    Slope is computed from central differences of the reprojected DEM.
    Because the target is a projected CRS with square metric pixels the
    calculation is straightforward.
    """
    with dem_memfile.open() as dem_src:
        with WarpedVRT(
            dem_src,
            crs=target_crs,
            transform=target_transform,
            width=target_shape[1],
            height=target_shape[0],
            resampling=Resampling.bilinear,
        ) as vrt:
            dem_masked = vrt.read(1, masked=True)

    missing_mask = np.ma.getmaskarray(dem_masked)
    dem = np.ascontiguousarray(dem_masked.filled(0.0), dtype=np.float32)

    px = abs(float(target_transform.a))
    py = abs(float(target_transform.e))
    grad_y, grad_x = np.gradient(dem, py, px)
    slope_deg = np.degrees(
        np.arctan(np.sqrt(grad_x * grad_x + grad_y * grad_y))
    ).astype(np.float32)
    return slope_deg, missing_mask


def mean_slope_per_label(
    slope_deg: np.ndarray,
    missing_mask: np.ndarray,
    labels: np.ndarray,
    n_labels: int,
) -> np.ndarray:
    """Return the mean slope (deg) for every label in ``[0, n_labels]``,
    ignoring pixels where the DEM was missing. NaN for labels with no
    valid DEM pixels."""
    label_ids = np.arange(n_labels + 1)
    slope_filled = np.where(missing_mask, 0.0, slope_deg).astype(np.float32)
    valid = (~missing_mask).astype(np.float32)
    slope_sum = ndimage.sum(slope_filled, labels=labels, index=label_ids)
    valid_count = ndimage.sum(valid, labels=labels, index=label_ids)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(valid_count > 0, slope_sum / valid_count, np.nan)
    return mean.astype(np.float32)


# ---------------------------------------------------------------------------
# Per-scene extraction
# ---------------------------------------------------------------------------
def extract_small_lakes(
    ndwi_path: Path,
    output_dir: Path,
    dem_memfile: MemoryFile | None = None,
) -> gpd.GeoDataFrame | None:
    """Extract polygons of small water bodies from a single NDWI GeoTIFF.

    Returns the resulting GeoDataFrame (in the raster's native CRS) or
    ``None`` if the scene name is not recognised or nothing survived
    the size / slope filters. The GeoDataFrame is also written to a
    per-scene GeoPackage inside ``output_dir``.
    """
    match = NAME_REGEX.search(ndwi_path.name)
    if match is None:
        tqdm.write(f"  Skipping (name not recognised): {ndwi_path.name}")
        return None

    date_uk = match.group("date")
    tile = match.group("tile")
    out_path = output_dir / f"{date_uk}_{tile}_smallLakes.gpkg"

    with rasterio.open(ndwi_path) as src:
        ndwi = src.read(1)
        transform = src.transform
        crs = src.crs

    water = np.isfinite(ndwi) & (ndwi > NDWI_THRESHOLD)
    if not water.any():
        tqdm.write(f"  {ndwi_path.name}: no water pixels above NDWI > {NDWI_THRESHOLD:g}")
        return None

    # 8-connectivity so diagonally-touching pixels join one lake.
    structure = np.ones((3, 3), dtype=bool)
    labels, n_labels = ndimage.label(water, structure=structure)
    sizes = np.bincount(labels.ravel(), minlength=n_labels + 1)
    sizes[0] = 0  # label 0 is background

    keep_size = (sizes >= MIN_LAKE_PIXELS) & (sizes <= MAX_LAKE_PIXELS)
    keep_size[0] = False
    if not keep_size.any():
        tqdm.write(
            f"  {ndwi_path.name}: 0 lakes in "
            f"{MIN_LAKE_AREA_KM2}-{MAX_LAKE_AREA_KM2} km^2 window "
            f"(largest component = {sizes.max() * PIXEL_AREA_M2 / 1e6:.3f} km^2)"
        )
        return None

    # ---- Slope filter ----------------------------------------------------
    if dem_memfile is not None:
        slope_deg, dem_missing = slope_on_target_grid(
            dem_memfile,
            target_crs=crs,
            target_transform=transform,
            target_shape=ndwi.shape,
        )
        mean_slope = mean_slope_per_label(
            slope_deg, dem_missing, labels, n_labels,
        )
        # NaN mean slope = no DEM coverage; reject conservatively.
        keep_slope = np.isfinite(mean_slope) & (mean_slope < SLOPE_MAX_DEG)
    else:
        mean_slope = np.full(len(sizes), np.nan, dtype=np.float32)
        keep_slope = np.ones(len(sizes), dtype=bool)

    keep = keep_size & keep_slope
    n_size = int(keep_size.sum())
    n_final = int(keep.sum())
    if n_final == 0:
        tqdm.write(
            f"  {ndwi_path.name}: {n_size} lake(s) matched size filter but "
            f"none had mean slope < {SLOPE_MAX_DEG:g} deg"
        )
        return None

    kept_mask = keep[labels]
    kept_labels = np.where(kept_mask, labels, 0).astype(np.int32)

    records: list[dict] = []
    for geom, value in shapes(
        kept_labels,
        mask=kept_mask,
        transform=transform,
        connectivity=8,
    ):
        label_id = int(value)
        if label_id == 0:
            continue
        polygon = shape(geom)
        n_pixels = int(sizes[label_id])
        slope_val = float(mean_slope[label_id])
        records.append(
            {
                "label_id": label_id,
                "n_pixels": n_pixels,
                "area_km2": n_pixels * PIXEL_AREA_M2 / 1_000_000.0,
                "mean_slope_deg": slope_val if np.isfinite(slope_val) else None,
                "date": date_uk,
                "tile": tile,
                "source": ndwi_path.name,
                "geometry": polygon,
            }
        )

    if not records:
        tqdm.write(f"  {ndwi_path.name}: no polygons produced from filtered labels")
        return None

    gdf = gpd.GeoDataFrame(records, crs=crs)
    output_dir.mkdir(parents=True, exist_ok=True)
    gdf.to_file(out_path, driver="GPKG", layer="small_lakes")

    dropped = n_size - n_final
    slope_note = (
        f" ({dropped} dropped by slope>{SLOPE_MAX_DEG:g}deg)"
        if dem_memfile is not None and dropped > 0
        else ""
    )
    tqdm.write(
        f"  {out_path.name}: {len(gdf):4d} lakes{slope_note} | "
        f"area {gdf['area_km2'].min():.3f} - {gdf['area_km2'].max():.3f} km^2 | "
        f"total {gdf['area_km2'].sum():.2f} km^2"
    )
    return gdf


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ndwi_files = sorted(NDWI_DIR.glob("*_NDWI.tif"))
    if not ndwi_files:
        sys.exit(f"No NDWI rasters found in {NDWI_DIR.resolve()}")

    print(f"Found {len(ndwi_files)} NDWI raster(s) in {NDWI_DIR.resolve()}")
    print(
        f"Water criterion: NDWI > {NDWI_THRESHOLD:g}, "
        f"area in [{MIN_LAKE_AREA_KM2:g}, {MAX_LAKE_AREA_KM2:g}] km^2 "
        f"({MIN_LAKE_PIXELS} - {MAX_LAKE_PIXELS} pixels @ 10 m), "
        f"mean slope < {SLOPE_MAX_DEG:g} deg"
    )
    print(f"Per-scene outputs -> {OUTPUT_DIR.resolve()}")

    print(f"Loading DEM tiles from {DEM_DIR.resolve()} ...")
    dem_memfile = load_dem_mosaic(DEM_DIR)
    if dem_memfile is None:
        print(
            f"  No DEM tiles found in {DEM_DIR}; slope filter (< "
            f"{SLOPE_MAX_DEG:g} deg) will be skipped."
        )
    else:
        print(f"  Slope filter active: keeping mean slope < {SLOPE_MAX_DEG:g} deg.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    per_scene: list[gpd.GeoDataFrame] = []
    try:
        for path in tqdm(ndwi_files, desc="Lakes", unit="scene"):
            try:
                gdf = extract_small_lakes(
                    path, OUTPUT_DIR, dem_memfile=dem_memfile,
                )
            except Exception as exc:  # noqa: BLE001
                tqdm.write(f"  Failed on {path.name}: {exc}")
                continue
            if gdf is not None and not gdf.empty:
                per_scene.append(gdf)
    finally:
        if dem_memfile is not None:
            dem_memfile.close()

    if not per_scene:
        print("No lakes matched the size + slope filters in any scene.")
        return

    # ---- Combine everything into one WGS84 GeoPackage --------------------
    reprojected = [
        (gdf.to_crs(COMBINED_CRS) if COMBINED_CRS else gdf) for gdf in per_scene
    ]
    combined = gpd.GeoDataFrame(
        pd.concat(reprojected, ignore_index=True),
        crs=reprojected[0].crs,
    )
    combined.to_file(COMBINED_OUTPUT, driver="GPKG", layer="small_lakes")

    print(
        f"\nCombined {len(combined):,} lake polygons across "
        f"{len(per_scene)} scene(s) -> {COMBINED_OUTPUT.resolve()}"
    )
    print(
        f"  Area: min={combined['area_km2'].min():.3f}  "
        f"median={combined['area_km2'].median():.3f}  "
        f"mean={combined['area_km2'].mean():.3f}  "
        f"max={combined['area_km2'].max():.3f}  "
        f"total={combined['area_km2'].sum():.2f} km^2"
    )
    print("Done.")


if __name__ == "__main__":
    main()
