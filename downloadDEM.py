"""Download the 30 m Copernicus DEM (GLO-30) tiles covering a
user-supplied Indian state.

At start-up the script prompts for the state / AOI name. That name is
used both to look up the state polygon in the GADM level-1 shapefile
and (in sanitised form) as the destination sub-folder under
``Data/DEM/``. Every 1-degree Copernicus DEM tile whose footprint
intersects the state polygon is downloaded from the public AWS Open
Data mirror:

    https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com/

No authentication is required.

Outputs
-------
* ``Data/DEM/<state>/Copernicus_DSM_COG_10_*_DEM.tif`` - one Cloud
  Optimised GeoTIFF per 1-degree tile, in geographic coordinates
  (EPSG:4326) with a nominal ground sampling of ~30 m at the equator.

Notes
-----
Copernicus DEM GLO-30 is derived from TanDEM-X and is distinct from the
NASA/USGS SRTM product. Ocean-only 1-degree cells are not published on
the AWS mirror and are silently skipped (HTTP 404).

Requirements::

    pip install geopandas shapely requests tqdm
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import geopandas as gpd
import requests
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SHAPEFILE = Path("../Data/gadm41_IND_shp/gadm41_IND_1.shp")
# Per-run tiles are written under ``OUTPUT_PARENT_DIR / <state>/`` where
# ``<state>`` is the sanitised name entered at start-up (see
# ``prompt_state`` below).
OUTPUT_PARENT_DIR = Path("../Data/DEM")

# Regional endpoint avoids the extra redirect through the global one.
BUCKET_URL = "https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com"


# ---------------------------------------------------------------------------
# Area of interest
# ---------------------------------------------------------------------------
def load_state_geometry(shapefile: Path, state_name: str) -> BaseGeometry:
    """Return the dissolved state boundary in WGS84 as a shapely geometry."""
    gdf = gpd.read_file(shapefile)

    name_field = next(
        (c for c in ("NAME_1", "NAME_2", "NAME_0") if c in gdf.columns),
        None,
    )
    if name_field is None:
        raise ValueError("Shapefile does not contain a recognised NAME_* column.")

    subset = gdf[gdf[name_field].str.casefold() == state_name.casefold()]
    if subset.empty:
        raise ValueError(
            f"State {state_name!r} not found in column {name_field!r} of "
            f"{shapefile}."
        )

    if subset.crs is None or subset.crs.to_epsg() != 4326:
        subset = subset.to_crs(epsg=4326)

    return unary_union(subset.geometry.values)


# ---------------------------------------------------------------------------
# Tile enumeration
# ---------------------------------------------------------------------------
def tile_stem(lat: int, lon: int) -> str:
    """Return the Copernicus DEM COG stem for the 1-degree tile whose
    south-west corner sits at (lat, lon)."""
    lat_prefix = "N" if lat >= 0 else "S"
    lon_prefix = "E" if lon >= 0 else "W"
    return (
        f"Copernicus_DSM_COG_10_{lat_prefix}{abs(lat):02d}_00_"
        f"{lon_prefix}{abs(lon):03d}_00_DEM"
    )


def tile_url(lat: int, lon: int) -> str:
    """Full HTTPS URL of the 1-degree Copernicus DEM tile on AWS."""
    stem = tile_stem(lat, lon)
    return f"{BUCKET_URL}/{stem}/{stem}.tif"


def tiles_intersecting(geometry: BaseGeometry) -> list[tuple[int, int]]:
    """Return all (lat, lon) 1-degree tiles that intersect the geometry."""
    minx, miny, maxx, maxy = geometry.bounds
    tiles: list[tuple[int, int]] = []
    for lat in range(math.floor(miny), math.ceil(maxy)):
        for lon in range(math.floor(minx), math.ceil(maxx)):
            if box(lon, lat, lon + 1, lat + 1).intersects(geometry):
                tiles.append((lat, lon))
    return sorted(tiles)


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
def download_tile(url: str, destination: Path) -> bool:
    """Stream one tile to disk. Returns False on 404 (tile not published)."""
    if destination.exists() and destination.stat().st_size > 0:
        tqdm.write(f"  Skipping (already downloaded): {destination.name}")
        return True

    with requests.get(url, stream=True, timeout=300) as response:
        if response.status_code == 404:
            tqdm.write(f"  Tile not published (void/ocean): {destination.stem}")
            return False
        response.raise_for_status()
        total = int(response.headers.get("Content-Length", 0))
        with open(destination, "wb") as fh, tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            desc=destination.name,
            leave=False,
        ) as bar:
            for chunk in response.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
                    bar.update(len(chunk))
    return True


# ---------------------------------------------------------------------------
# Interactive prompt
# ---------------------------------------------------------------------------
def prompt_state(
    shapefile: Path, output_parent: Path
) -> tuple[str, BaseGeometry, Path]:
    """Prompt for a state / AOI name; return ``(display_name, geometry,
    output_dir)``.

    ``display_name`` is the raw user input (used both to look up the
    state polygon in ``shapefile`` and in printed messages).
    ``output_dir`` is ``output_parent / <sanitised name>`` where the
    name is sanitised to keep only ``[A-Za-z0-9_-]`` characters so it
    is always a safe folder name on any filesystem. The loop re-asks
    if the entered name cannot be matched in the shapefile.
    """
    try:
        while True:
            name = input(
                f"Name of the state / AOI (sub-folder of {output_parent}): "
            ).strip()
            if not name:
                continue
            safe = re.sub(r"[^\w\-]+", "_", name).strip("_")
            if not safe:
                print("  Name contains no usable characters. Try again.")
                continue
            try:
                geometry = load_state_geometry(shapefile, name)
            except ValueError as exc:
                print(f"  {exc}")
                continue
            return name, geometry, output_parent / safe
    except (EOFError, KeyboardInterrupt):
        sys.exit("\nAborted: no state provided.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    state_name, aoi, output_dir = prompt_state(SHAPEFILE, OUTPUT_PARENT_DIR)
    minx, miny, maxx, maxy = aoi.bounds
    print(f"\nArea of interest: {state_name}")
    print(
        f"  Bounding box: {miny:.2f}-{maxy:.2f} N, {minx:.2f}-{maxx:.2f} E"
    )

    tiles = tiles_intersecting(aoi)
    if not tiles:
        sys.exit(
            f"No Copernicus DEM tiles intersect the requested state "
            f"({state_name})."
        )
    print(f"  {len(tiles)} 1-degree tile(s) to download:")
    for lat, lon in tiles:
        print(f"    {tile_stem(lat, lon)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nDownloading to {output_dir.resolve()} ...")
    for lat, lon in tqdm(tiles, desc="DEM tiles", unit="tile"):
        stem = tile_stem(lat, lon)
        destination = output_dir / f"{stem}.tif"
        try:
            download_tile(tile_url(lat, lon), destination)
        except requests.HTTPError as exc:
            tqdm.write(f"  Failed on {stem}: {exc}")

    print("Done.")


if __name__ == "__main__":
    main()
