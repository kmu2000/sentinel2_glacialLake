"""Compute cloud-masked NDWI rasters for every Sentinel-2 L1C scene
below ``Data/Sentinel/`` and report NDWI statistics for pixels falling
inside the Zhang et al. (2022) Himalayan glacial-lake polygons.

At start-up the script prompts for the state / AOI name; the resulting
NDWI GeoTIFFs are written into ``Data/Sentinel/NDWI/<state>/`` (the
state name is sanitised so it is always a safe folder name).

For each ``S2X_MSIL1C_*.zip`` archive the script:

1. Opens the Green (B03) and NIR (B08) bands directly from the zip
   using GDAL's ``/vsizip/`` virtual filesystem and converts digital
   numbers to top-of-atmosphere reflectance using the per-band
   ``RADIO_ADD_OFFSET`` values declared in ``MTD_MSIL1C.xml`` (offsets
   were introduced in processing baseline 04.00).
2. Computes McFeeters' NDWI = (B03 - B08) / (B03 + B08).
3. Applies the L1C ``MSK_CLASSI_B00.jp2`` cloud/cirrus mask: cloudy
   pixels, along with L1C nodata pixels, are written as NaN.
4. Reprojects the Zhang 2022 lake polygons to the scene UTM grid,
   selects the polygons intersecting the tile, rasterises them onto the
   NDWI grid, and prints summary statistics (count, min, max, mean,
   std, median, and a handful of percentiles) for all valid NDWI pixels
   that fall inside those polygons.
5. Accumulates the per-scene "in-lake" pixels and prints an overall
   summary across every processed scene at the end.

Outputs
-------
* ``Data/Sentinel/NDWI/<state>/DD-MM-YYYY_<TILE>_NDWI.tif`` - Float32
  NDWI with NaN wherever the pixel is nodata or flagged as cloud /
  cirrus.
* Per-scene and overall NDWI-inside-lake statistics printed to stdout.

Requirements::

    pip install rasterio numpy geopandas tqdm
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Callable

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.transform import array_bounds
from shapely.geometry import box
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INPUT_DIR = Path("../Data/Sentinel")
# NDWI GeoTIFFs are written under ``NDWI_PARENT_DIR / <state>/`` where
# ``<state>`` is the sanitised name entered at start-up (see
# ``prompt_state_folder`` below).
NDWI_PARENT_DIR = INPUT_DIR / "NDWI"

# Reference glacial-lake polygons (Zhang et al., 2022, Himalayan region).
# The stats block is skipped silently if this file is not present.
LAKES_SHP = Path("../Data/Zhang2022_Data/GL/Himalayan_Glacial_Lakes_2020.shp")

GREEN_BAND = "B03"  # Sentinel-2 Green, native 10 m
NIR_BAND = "B08"    # Sentinel-2 NIR,   native 10 m

QUANTIFICATION_VALUE = 10000.0  # DN -> reflectance scale factor
DEFAULT_OFFSET = 0.0            # DN offset for pre-baseline-04 products

# Zero-based band_id used in the RADIO_ADD_OFFSET nodes of MTD_MSIL1C.xml.
BAND_INDEX = {
    "B01": 0, "B02": 1, "B03": 2, "B04": 3, "B05": 4, "B06": 5,
    "B07": 6, "B08": 7, "B8A": 8, "B09": 9, "B10": 10, "B11": 11,
    "B12": 12,
}

NAME_REGEX = re.compile(
    r"S2[ABC]_MSIL1C_(?P<date>\d{8})T\d{6}_.*_T(?P<tile>\d{2}[A-Z]{3})_"
)


# ---------------------------------------------------------------------------
# Zip helpers
# ---------------------------------------------------------------------------
def _list_zip(zip_path: Path) -> list[str]:
    with zipfile.ZipFile(zip_path) as zf:
        return zf.namelist()


def _find_member(names: list[str], predicate: Callable[[str], bool]) -> str | None:
    for name in names:
        if predicate(name):
            return name
    return None


def _vsizip(zip_path: Path, inner: str) -> str:
    """Build a GDAL ``/vsizip/`` URI for a member inside the archive."""
    return f"/vsizip/{zip_path.resolve().as_posix()}/{inner}"


# ---------------------------------------------------------------------------
# Radiometric metadata
# ---------------------------------------------------------------------------
def read_radiometric_offsets(zip_path: Path, product_xml: str) -> dict[int, float]:
    """Return the ``RADIO_ADD_OFFSET`` values keyed by ``band_id`` (0-based)."""
    with zipfile.ZipFile(zip_path) as zf, zf.open(product_xml) as fh:
        tree = ET.parse(fh)
    offsets: dict[int, float] = {}
    for node in tree.iter():
        if not node.tag.endswith("RADIO_ADD_OFFSET"):
            continue
        try:
            band_id = int(node.attrib.get("band_id", "-1"))
            offsets[band_id] = float((node.text or "0").strip())
        except (TypeError, ValueError):
            continue
    return offsets


# ---------------------------------------------------------------------------
# Raster IO
# ---------------------------------------------------------------------------
def write_float_geotiff(
    array: np.ndarray,
    path: Path,
    transform,
    crs,
    description: str,
) -> None:
    """Write ``array`` as a DEFLATE-compressed Float32 single-band GeoTIFF."""
    path.parent.mkdir(parents=True, exist_ok=True)
    height, width = array.shape
    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "count": 1,
        "nodata": np.nan,
        "width": width,
        "height": height,
        "transform": transform,
        "crs": crs,
        "compress": "deflate",
        "predictor": 3,  # floating-point predictor
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype(np.float32), 1)
        dst.set_band_description(1, description)


# ---------------------------------------------------------------------------
# Lake-polygon NDWI statistics
# ---------------------------------------------------------------------------
def load_lake_polygons(shp_path: Path) -> gpd.GeoDataFrame | None:
    """Load the reference lake polygons; return ``None`` if missing."""
    if not shp_path.exists():
        return None
    gdf = gpd.read_file(shp_path)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
    if gdf.crs is None:
        tqdm.write(
            f"  Warning: {shp_path.name} has no CRS; assuming EPSG:4326."
        )
        gdf = gdf.set_crs("EPSG:4326")
    return gdf.reset_index(drop=True)


def extract_ndwi_in_lakes(
    ndwi: np.ndarray,
    transform,
    crs,
    lakes: gpd.GeoDataFrame,
) -> tuple[np.ndarray, int]:
    """Return (finite NDWI values inside lake polygons, number of polygons
    that intersect the raster extent).

    Polygons are reprojected to the raster CRS, filtered to those
    intersecting the raster bbox, then rasterised onto the raster grid
    using ``rasterio.features.geometry_mask``.
    """
    height, width = ndwi.shape
    west, south, east, north = array_bounds(height, width, transform)
    tile_box = box(west, south, east, north)

    lakes_local = lakes.to_crs(crs)
    intersecting = lakes_local[lakes_local.geometry.intersects(tile_box)]
    if intersecting.empty:
        return np.empty(0, dtype=np.float32), 0

    outside = geometry_mask(
        intersecting.geometry.values,
        out_shape=(height, width),
        transform=transform,
        invert=False,
        all_touched=False,
    )
    inside = ~outside
    values = ndwi[inside]
    return values[np.isfinite(values)], len(intersecting)


def format_ndwi_stats(label: str, values: np.ndarray) -> str:
    """Format min/mean/median/percentile summary of NDWI ``values``."""
    if values.size == 0:
        return f"{label} NDWI-in-lakes: no valid pixels"
    p1, p5, p25, p50, p75, p95, p99 = np.percentile(
        values, [1, 5, 25, 50, 75, 95, 99]
    )
    return (
        f"{label} NDWI-in-lakes: n={values.size:,}  "
        f"min={float(values.min()):+.3f}  P1={p1:+.3f}  P5={p5:+.3f}  "
        f"P25={p25:+.3f}  median={p50:+.3f}  mean={float(values.mean()):+.3f}  "
        f"P75={p75:+.3f}  P95={p95:+.3f}  P99={p99:+.3f}  "
        f"max={float(values.max()):+.3f}  std={float(values.std()):.3f}"
    )


# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------
def process_scene(
    zip_path: Path,
    output_dir: Path,
    lakes: gpd.GeoDataFrame | None = None,
) -> tuple[Path | None, np.ndarray]:
    """Compute cloud-masked NDWI for a single L1C zip and write a GeoTIFF.

    Returns ``(output_path, in_lake_ndwi_values)``; the second element is
    an empty array if ``lakes`` is None or nothing overlaps this tile.
    """
    empty = np.empty(0, dtype=np.float32)

    match = NAME_REGEX.search(zip_path.name)
    if match is None:
        tqdm.write(f"  Skipping (name not recognised): {zip_path.name}")
        return None, empty

    date_stamp = match.group("date")  # YYYYMMDD
    tile = match.group("tile")
    date_uk = f"{date_stamp[6:8]}-{date_stamp[4:6]}-{date_stamp[0:4]}"
    output_path = output_dir / f"{date_uk}_{tile}_NDWI.tif"
    if output_path.exists():
        tqdm.write(f"  Skipping (exists): {output_path.name}")
        lake_values = empty
        if lakes is not None:
            with rasterio.open(output_path) as src:
                lake_values, n_polys = extract_ndwi_in_lakes(
                    src.read(1), src.transform, src.crs, lakes,
                )
            tqdm.write(
                "  " + format_ndwi_stats(
                    f"{output_path.name} ({n_polys} polygon(s))", lake_values
                )
            )
        return output_path, lake_values

    names = _list_zip(zip_path)
    product_xml = _find_member(names, lambda n: n.endswith("/MTD_MSIL1C.xml"))
    green_jp2 = _find_member(
        names,
        lambda n: "/IMG_DATA/" in n and n.endswith(f"_{GREEN_BAND}.jp2"),
    )
    nir_jp2 = _find_member(
        names,
        lambda n: "/IMG_DATA/" in n and n.endswith(f"_{NIR_BAND}.jp2"),
    )
    mask_jp2 = _find_member(names, lambda n: n.endswith("MSK_CLASSI_B00.jp2"))

    if green_jp2 is None or nir_jp2 is None:
        tqdm.write(f"  Missing B03/B08 in {zip_path.name}; skipping.")
        return None, empty

    offsets: dict[int, float] = {}
    if product_xml is not None:
        try:
            offsets = read_radiometric_offsets(zip_path, product_xml)
        except Exception as exc:  # noqa: BLE001
            tqdm.write(
                f"  Warning: could not read radiometric offsets for "
                f"{zip_path.name} ({exc}); assuming 0."
            )
    green_off = offsets.get(BAND_INDEX[GREEN_BAND], DEFAULT_OFFSET)
    nir_off = offsets.get(BAND_INDEX[NIR_BAND], DEFAULT_OFFSET)

    # Read the Green band (10 m) and use its grid as the reference.
    with rasterio.open(_vsizip(zip_path, green_jp2)) as src:
        green = src.read(1).astype(np.float32)
        target_shape = green.shape
        target_transform = src.transform
        target_crs = src.crs

    # Read the NIR band (native 10 m, same grid as Green).
    with rasterio.open(_vsizip(zip_path, nir_jp2)) as src:
        nir = src.read(1).astype(np.float32)

    nodata_mask = (green == 0) | (nir == 0)

    green_ref = (green + green_off) / QUANTIFICATION_VALUE
    nir_ref = (nir + nir_off) / QUANTIFICATION_VALUE

    with np.errstate(divide="ignore", invalid="ignore"):
        denom = green_ref + nir_ref
        ndwi = np.where(
            np.abs(denom) > 0.0,
            (green_ref - nir_ref) / denom,
            np.nan,
        )

    # L1C opaque cloud + cirrus mask (60 m, resampled to the reference grid).
    if mask_jp2 is not None:
        with rasterio.open(_vsizip(zip_path, mask_jp2)) as src:
            classi_layers = src.read(
                [1, 2],  # opaque cloud, cirrus cloud
                out_shape=(2, *target_shape),
                resampling=Resampling.nearest,
            )
        cloud_mask = (classi_layers[0] > 0) | (classi_layers[1] > 0)
    else:
        tqdm.write(
            f"  No MSK_CLASSI mask in {zip_path.name}; NDWI will not be "
            "cloud-masked."
        )
        cloud_mask = np.zeros(target_shape, dtype=bool)

    invalid = nodata_mask | cloud_mask
    ndwi_out = np.where(invalid, np.nan, ndwi).astype(np.float32)

    write_float_geotiff(
        ndwi_out,
        output_path,
        transform=target_transform,
        crs=target_crs,
        description="NDWI = (B03 - B08) / (B03 + B08); cloud/cirrus = NaN",
    )

    n_valid = int(np.isfinite(ndwi_out).sum())
    pct_invalid = (ndwi_out.size - n_valid) / ndwi_out.size * 100.0
    tqdm.write(
        f"  {output_path.name}: {pct_invalid:5.2f}% invalid "
        "(nodata + cloud/cirrus)"
    )

    lake_values = np.empty(0, dtype=np.float32)
    if lakes is not None:
        lake_values, n_polys = extract_ndwi_in_lakes(
            ndwi_out, target_transform, target_crs, lakes,
        )
        tqdm.write(
            "  " + format_ndwi_stats(
                f"{output_path.name} ({n_polys} polygon(s))", lake_values
            )
        )

    return output_path, lake_values


# ---------------------------------------------------------------------------
# Interactive prompt
# ---------------------------------------------------------------------------
def prompt_state_folder(parent: Path) -> tuple[str, Path]:
    """Prompt for the state / AOI name; return ``(display_name, output_dir)``.

    ``display_name`` is the raw user input (used in printed messages).
    ``output_dir`` is ``parent / <sanitised name>`` where the name is
    sanitised to keep only ``[A-Za-z0-9_-]`` characters so it is always
    a safe folder name on any filesystem.
    """
    try:
        while True:
            name = input(
                f"Name of the state / AOI (sub-folder of {parent}): "
            ).strip()
            if not name:
                continue
            safe = re.sub(r"[^\w\-]+", "_", name).strip("_")
            if not safe:
                print("  Name contains no usable characters. Try again.")
                continue
            return name, parent / safe
    except (EOFError, KeyboardInterrupt):
        sys.exit("\nAborted: no output folder provided.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    state_name, output_dir = prompt_state_folder(NDWI_PARENT_DIR)

    # rglob so we pick up scenes both at the root of ``Data/Sentinel/``
    # (older layout) and inside an AOI sub-folder created by
    # ``downloadSentinel2.py`` (e.g. ``Data/Sentinel/<state>/*.zip``).
    scenes = sorted(INPUT_DIR.rglob("S2*_MSIL1C_*.zip"))
    if not scenes:
        sys.exit(
            f"No Sentinel-2 L1C zip archives found under {INPUT_DIR.resolve()}"
        )

    print(f"\nFound {len(scenes)} scene(s) under {INPUT_DIR.resolve()}")
    print(f"NDWI outputs ({state_name}) -> {output_dir.resolve()}")

    lakes = load_lake_polygons(LAKES_SHP)
    if lakes is None:
        print(
            f"Lake polygons {LAKES_SHP} not found; NDWI-in-lake statistics "
            "will be skipped."
        )
    else:
        print(
            f"Loaded {len(lakes)} lake polygon(s) from {LAKES_SHP} "
            f"(source CRS: {lakes.crs})"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    all_lake_values: list[np.ndarray] = []
    for scene in tqdm(scenes, desc="NDWI", unit="scene"):
        try:
            _path, lake_values = process_scene(scene, output_dir, lakes=lakes)
            if lake_values.size:
                all_lake_values.append(lake_values)
        except Exception as exc:  # noqa: BLE001
            tqdm.write(f"  Failed on {scene.name}: {exc}")

    if lakes is not None:
        combined = (
            np.concatenate(all_lake_values)
            if all_lake_values
            else np.empty(0, dtype=np.float32)
        )
        print()
        print(format_ndwi_stats("Overall (all scenes)", combined))

    print("Done.")


if __name__ == "__main__":
    main()
