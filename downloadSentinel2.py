"""Download the minimum set of least-cloudy Sentinel-2 scenes needed to
cover a rectangular area of interest defined by top-left and bottom-right
(longitude, latitude) corners.

Workflow
--------
1. Prompt the user for the name of the state / AOI (used both as a
   sub-folder under ``../Data/Sentinel/`` and in printed messages), and
   for the upper-left and lower-right corners of the area of interest
   (WGS84 lon/lat).
2. Query the Copernicus Data Space Ecosystem (CDSE) for Sentinel-2
   scenes within the configured date window whose cloud cover is <=
   ``MAX_CLOUD_COVER``, and keep the least-cloudy scene per MGRS tile.
3. Reduce the tile set to the minimum needed to cover the AOI via
   greedy set-cover, breaking ties by lower cloud cover.
4. Report the plan (number of scenes, per-tile cloud cover, coverage
   percentage of the AOI, and any remaining "holes"), skipping tiles
   already on disk.
5. Ask the user for confirmation before downloading.

If the reported coverage is not high enough, widen ``START_DATE`` /
``END_DATE`` or raise ``MAX_CLOUD_COVER`` manually and re-run.

Setup
-----
1. Register (free) at https://dataspace.copernicus.eu/.
2. Run the script; it will prompt for your CDSE username and password on
   stdin (the password is read via ``getpass`` so it is not echoed).
3. Install the required Python packages::

       pip install shapely pyproj requests tqdm
"""

from __future__ import annotations

import getpass
import re
import sys
from datetime import datetime
from pathlib import Path

import pyproj
import requests
from shapely import wkt as shapely_wkt
from shapely.geometry import box, shape as shapely_shape
from shapely.geometry.base import BaseGeometry
from shapely.ops import transform as shapely_transform, unary_union
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# The area of interest and the output sub-folder name are entered
# interactively at run-time (see ``prompt_bbox`` and
# ``prompt_output_folder`` below).

START_DATE = "2025-08-01"
END_DATE = "2025-09-30"
# Upper bound on cloudCover in the CDSE catalogue query. Kept at 100 by
# default so that no MGRS tile is silently dropped just because every
# scene of that tile happens to be cloudy; the per-tile "least cloudy"
# selection below still picks the cleanest available scene for each tile.
# Lower this only if you want to hard-exclude very cloudy candidates.
MAX_CLOUD_COVER = 100.0
# Warn (do not filter) when a tile's best scene is cloudier than this.
CLOUD_WARN_PCT = 30.0

# Inward buffer (metres) applied to the AOI polygon before the coverage
# check. Useful when the AOI comes from a coarser boundary (e.g. a
# generalised administrative shapefile) whose vertices do not line up
# with the Sentinel-2 tile footprints and would otherwise leave
# hair-thin slivers uncovered. For a manually specified rectangular
# bounding box this is not needed and can stay at 0.
AOI_INWARD_BUFFER_M = 0.0

PRODUCT_TYPE = "S2MSI1C"  # Sentinel-2 Level-1C: orthorectified TOA reflectance
# Parent directory under which each run creates (or reuses) a sub-folder
# named after the state / AOI entered at the prompt.
OUTPUT_PARENT_DIR = Path("../Data/Sentinel")

CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)
CDSE_CATALOG_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
CDSE_DOWNLOAD_URL = (
    "https://catalogue.dataspace.copernicus.eu/odata/v1/Products({pid})/$value"
)

# Sentinel-2 product names embed the MGRS tile as "_Tddccc_" (e.g. "_T43SGA_").
TILE_REGEX = re.compile(r"_T(\d{2}[A-Z]{3})_")


# ---------------------------------------------------------------------------
# Authentication
# ---------------------------------------------------------------------------
def prompt_credentials() -> tuple[str, str]:
    """Prompt the user for CDSE username and password on stdin.

    The password is read via :func:`getpass.getpass` so it is not echoed.
    Both prompts re-ask on empty input; a KeyboardInterrupt or EOF from
    the user aborts the whole script cleanly.
    """
    print(
        "Copernicus Data Space credentials "
        "(register free at https://dataspace.copernicus.eu/):"
    )
    try:
        username = input("  Username (email): ").strip()
        while not username:
            username = input("  Username cannot be empty. Try again: ").strip()
        password = getpass.getpass("  Password: ")
        while not password:
            password = getpass.getpass("  Password cannot be empty. Try again: ")
    except (EOFError, KeyboardInterrupt):
        sys.exit("\nAborted: no credentials provided.")
    return username, password


def get_access_token(username: str, password: str) -> str:
    """Request an OAuth2 access token from the Copernicus Data Space."""
    data = {
        "client_id": "cdse-public",
        "grant_type": "password",
        "username": username,
        "password": password,
    }
    response = requests.post(CDSE_TOKEN_URL, data=data, timeout=60)
    response.raise_for_status()
    return response.json()["access_token"]


# ---------------------------------------------------------------------------
# Area of interest
# ---------------------------------------------------------------------------
# CDSE truncates OData requests whose URL exceeds ~8 kB, so we keep the WKT
# well below that limit to leave headroom for the rest of the query string.
MAX_WKT_LEN = 3500


def load_bbox_geometry(
    top_left: tuple[float, float],
    bottom_right: tuple[float, float],
) -> BaseGeometry:
    """Return the AOI rectangle in WGS84 as a shapely polygon.

    ``top_left`` and ``bottom_right`` are ``(longitude, latitude)`` tuples
    in degrees. Top-left must have a smaller longitude and a larger
    latitude than bottom-right; the function raises ``ValueError``
    otherwise.
    """
    tl_lon, tl_lat = top_left
    br_lon, br_lat = bottom_right
    if not (-180.0 <= tl_lon <= 180.0 and -180.0 <= br_lon <= 180.0):
        raise ValueError(
            f"Longitudes must be in [-180, 180]; got tl={tl_lon}, br={br_lon}."
        )
    if not (-90.0 <= tl_lat <= 90.0 and -90.0 <= br_lat <= 90.0):
        raise ValueError(
            f"Latitudes must be in [-90, 90]; got tl={tl_lat}, br={br_lat}."
        )
    if tl_lon >= br_lon or tl_lat <= br_lat:
        raise ValueError(
            "Upper-left corner must be north-west of the lower-right "
            "corner: expected upper-left longitude < lower-right longitude "
            "and upper-left latitude > lower-right latitude; got "
            f"upper-left={top_left}, lower-right={bottom_right}."
        )
    return box(tl_lon, br_lat, br_lon, tl_lat)


def _read_corner(prompt: str) -> tuple[float, float]:
    """Read one ``'lon lat'`` (space- or comma-separated) pair from stdin.

    Re-asks on blank input, wrong token count, or unparseable numbers.
    ``EOFError`` / ``KeyboardInterrupt`` are propagated so callers can
    treat Ctrl-D / Ctrl-C as an abort.
    """
    while True:
        raw = input(prompt).strip()
        if not raw:
            continue
        parts = [p for p in re.split(r"[,\s]+", raw) if p]
        if len(parts) != 2:
            print("  Expected two numbers (longitude latitude). Try again.")
            continue
        try:
            return (float(parts[0]), float(parts[1]))
        except ValueError:
            print("  Both values must be numbers. Try again.")


def prompt_output_folder(parent: Path) -> tuple[str, Path]:
    """Prompt for the state / AOI name; return ``(display_name, output_dir)``.

    ``display_name`` is the raw user input (used in printed messages).
    ``output_dir`` is ``parent / <sanitised name>`` where the name is
    sanitised to keep only ``[A-Za-z0-9_-]`` characters, so it is always
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


def prompt_bbox() -> tuple[tuple[float, float], tuple[float, float]]:
    """Prompt the user for the AOI bounding box on stdin.

    Returns ``((upper_left_lon, upper_left_lat), (lower_right_lon,
    lower_right_lat))``. Re-asks until the pair passes the checks in
    :func:`load_bbox_geometry` (in-range coordinates that form a proper
    NW / SE rectangle).
    """
    print(
        "Area of interest bounding box (WGS84 lon/lat degrees; "
        "upper-left = NW corner, lower-right = SE corner):"
    )
    try:
        while True:
            top_left = _read_corner("  Upper-left  'lon lat': ")
            bottom_right = _read_corner("  Lower-right 'lon lat': ")
            try:
                load_bbox_geometry(top_left, bottom_right)
            except ValueError as exc:
                print(f"  {exc}\n  Please re-enter both corners.")
                continue
            return top_left, bottom_right
    except (EOFError, KeyboardInterrupt):
        sys.exit("\nAborted: no bounding box provided.")


def build_search_wkt(geometry: BaseGeometry, max_len: int = MAX_WKT_LEN) -> str:
    """Return a WKT polygon that fully contains ``geometry`` and is short
    enough to fit inside a CDSE OData request.

    Prefers the convex hull (which always contains the input) and falls
    back to the axis-aligned envelope if the hull's WKT still overflows
    the URL. Any extra tiles returned by the coarser polygon are trivially
    pruned by the greedy set-cover step.
    """
    candidates: list[tuple[str, BaseGeometry]] = [
        ("convex hull", geometry.convex_hull),
        ("bounding box", geometry.envelope),
    ]

    for label, candidate in candidates:
        wkt = candidate.wkt
        if len(wkt) <= max_len:
            print(f"  Using {label} for the search polygon ({len(wkt)} chars)")
            return wkt

    raise RuntimeError(
        "Could not reduce the area of interest to a WKT small enough for a "
        "CDSE OData query."
    )


# ---------------------------------------------------------------------------
# Catalogue search
# ---------------------------------------------------------------------------
def search_products(
    aoi_wkt: str,
    start_date: str,
    end_date: str,
    product_type: str,
    max_cloud_cover: float,
) -> list[dict]:
    """Page through the CDSE catalogue and return every matching product."""
    filter_expr = (
        "Collection/Name eq 'SENTINEL-2' "
        f"and OData.CSC.Intersects(area=geography'SRID=4326;{aoi_wkt}') "
        f"and ContentDate/Start ge {start_date}T00:00:00.000Z "
        f"and ContentDate/Start lt {end_date}T23:59:59.999Z "
        "and Attributes/OData.CSC.StringAttribute/any("
        "att:att/Name eq 'productType' "
        f"and att/OData.CSC.StringAttribute/Value eq '{product_type}') "
        "and Attributes/OData.CSC.DoubleAttribute/any("
        "att:att/Name eq 'cloudCover' "
        f"and att/OData.CSC.DoubleAttribute/Value le {max_cloud_cover:.2f})"
    )
    params = {
        "$filter": filter_expr,
        "$expand": "Attributes",
        "$orderby": "ContentDate/Start asc",
        "$top": 1000,
    }

    products: list[dict] = []
    url: str | None = CDSE_CATALOG_URL
    is_first_page = True
    while url:
        response = requests.get(
            url,
            params=params if is_first_page else None,
            timeout=120,
        )
        response.raise_for_status()
        payload = response.json()
        products.extend(payload.get("value", []))
        url = payload.get("@odata.nextLink")
        is_first_page = False

    return products


def _cloud_cover(product: dict) -> float:
    for attribute in product.get("Attributes", []):
        if attribute.get("Name") == "cloudCover":
            try:
                return float(attribute["Value"])
            except (TypeError, ValueError):
                pass
    return 100.0


def _tile_id(product: dict) -> str | None:
    match = TILE_REGEX.search(product.get("Name", ""))
    return match.group(1) if match else None


def select_least_cloudy_per_tile(products: list[dict]) -> list[dict]:
    """Return one scene per MGRS tile, keeping the lowest cloud cover."""
    best_by_tile: dict[str, dict] = {}
    for product in products:
        tile = _tile_id(product)
        if tile is None:
            continue
        current_best = best_by_tile.get(tile)
        if current_best is None or _cloud_cover(product) < _cloud_cover(current_best):
            best_by_tile[tile] = product
    return sorted(best_by_tile.values(), key=lambda p: _tile_id(p) or "")


# ---------------------------------------------------------------------------
# Coverage verification
# ---------------------------------------------------------------------------
# Equal-area projection used for area/percentage calculations. EPSG:6933 is
# the WGS 84 / NSIDC EASE-Grid 2.0 Global equal-area CRS; areas come out in
# square metres, which we convert to km^2 below.
_TO_EQUAL_AREA = pyproj.Transformer.from_crs(
    "EPSG:4326", "EPSG:6933", always_xy=True
).transform
_FROM_EQUAL_AREA = pyproj.Transformer.from_crs(
    "EPSG:6933", "EPSG:4326", always_xy=True
).transform


def _to_equal_area(geometry: BaseGeometry) -> BaseGeometry:
    return shapely_transform(_TO_EQUAL_AREA, geometry)


def _from_equal_area(geometry: BaseGeometry) -> BaseGeometry:
    return shapely_transform(_FROM_EQUAL_AREA, geometry)


def effective_aoi(geometry: BaseGeometry) -> BaseGeometry:
    """Return the AOI polygon shrunk inward by ``AOI_INWARD_BUFFER_M``
    metres, for use in the coverage check. Passthrough when the buffer is
    non-positive.
    """
    if AOI_INWARD_BUFFER_M <= 0:
        return geometry
    buffered_ea = _to_equal_area(geometry).buffer(-AOI_INWARD_BUFFER_M)
    if buffered_ea.is_empty:
        return geometry
    return _from_equal_area(buffered_ea)


def _parse_footprint(product: dict) -> BaseGeometry | None:
    """Return the WGS84 footprint of a CDSE product, or None if unavailable."""
    geo = product.get("GeoFootprint")
    if isinstance(geo, dict):
        try:
            return shapely_shape(geo)
        except Exception:  # noqa: BLE001
            pass
    fp_str = product.get("Footprint")
    if isinstance(fp_str, str):
        # CDSE returns strings like "geography'SRID=4326;POLYGON((...))'".
        # Strip the OData wrapper down to the bare WKT.
        wkt_match = re.search(r"(MULTIPOLYGON|POLYGON)\s*\(.*\)", fp_str, re.DOTALL)
        if wkt_match:
            try:
                return shapely_wkt.loads(wkt_match.group(0))
            except Exception:  # noqa: BLE001
                pass
    return None


def report_aoi_coverage(
    aoi_geometry: BaseGeometry,
    products: list[dict],
) -> BaseGeometry | None:
    """Print how much of the AOI is covered by the selected products.

    Returns the uncovered geometry (may be empty), or ``None`` if none of
    the products expose a usable footprint.
    """
    footprints = [fp for p in products if (fp := _parse_footprint(p)) is not None]
    missing = len(products) - len(footprints)
    if not footprints:
        print("  Coverage check skipped: no product footprints available.")
        return None
    if missing:
        print(
            f"  Note: {missing} of {len(products)} selected products had no "
            "footprint metadata; coverage below is a lower bound."
        )

    aoi = effective_aoi(aoi_geometry)
    if AOI_INWARD_BUFFER_M > 0:
        print(
            f"  Coverage measured against a {AOI_INWARD_BUFFER_M:g} m "
            "inward-buffered AOI polygon (to ignore boundary-precision noise)."
        )

    union = unary_union(footprints)
    aoi_ea = _to_equal_area(aoi)
    intersect_ea = _to_equal_area(aoi.intersection(union))
    uncovered = aoi.difference(union)
    uncovered_ea = _to_equal_area(uncovered)

    aoi_area_km2 = aoi_ea.area / 1_000_000.0
    covered_km2 = intersect_ea.area / 1_000_000.0
    uncovered_km2 = uncovered_ea.area / 1_000_000.0
    covered_pct = 100.0 * covered_km2 / aoi_area_km2 if aoi_area_km2 else 0.0

    print(
        f"  AOI area:   {aoi_area_km2:8,.1f} km^2   "
        f"covered: {covered_km2:8,.1f} km^2 ({covered_pct:5.2f}%)   "
        f"uncovered: {uncovered_km2:8,.1f} km^2"
    )

    if uncovered.is_empty or uncovered_km2 < 1e-3:
        print("  Full AOI coverage.")
        return uncovered

    gap_geoms = list(getattr(uncovered, "geoms", [uncovered]))
    # Rank gaps largest-first and print the top few.
    gap_geoms.sort(
        key=lambda g: _to_equal_area(g).area,
        reverse=True,
    )
    top = min(len(gap_geoms), 5)
    print(f"  {len(gap_geoms)} uncovered polygon(s); top {top} by area:")
    for gap in gap_geoms[:top]:
        area_km2 = _to_equal_area(gap).area / 1_000_000.0
        minx, miny, maxx, maxy = gap.bounds
        print(
            f"    {area_km2:8.2f} km^2  bbox=({miny:.3f},{minx:.3f}) - "
            f"({maxy:.3f},{maxx:.3f})"
        )
    return uncovered


def aoi_coverage_pct(
    aoi_geometry: BaseGeometry,
    products: list[dict],
) -> float:
    """Return the equal-area coverage percentage of the given products."""
    footprints = [fp for p in products if (fp := _parse_footprint(p)) is not None]
    if not footprints:
        return 0.0
    aoi = effective_aoi(aoi_geometry)
    aoi_ea = _to_equal_area(aoi)
    if aoi_ea.area <= 0:
        return 0.0
    union = unary_union(footprints)
    return 100.0 * _to_equal_area(aoi.intersection(union)).area / aoi_ea.area


def minimum_covering_scenes(
    aoi_geometry: BaseGeometry,
    candidates: list[dict],
) -> list[dict]:
    """Return the smallest subset of ``candidates`` whose footprints jointly
    cover the AOI geometry.

    Greedy set-cover: at each step the scene that adds the largest
    equal-area intersection with the remaining uncovered part of the AOI
    is picked; ties are broken by lower cloud cover.
    """
    parsed: list[tuple[dict, BaseGeometry]] = []
    for product in candidates:
        footprint = _parse_footprint(product)
        if footprint is not None:
            parsed.append((product, footprint))
    if not parsed:
        return []

    aoi = effective_aoi(aoi_geometry)
    aoi_area = _to_equal_area(aoi).area
    if aoi_area <= 0:
        return []

    selected: list[dict] = []
    remaining = aoi
    # Stop when the still-uncovered fraction is negligible (guards against
    # infinite loops when polygons intersect only in geodesic slivers).
    min_gain_frac = 1e-4
    while _to_equal_area(remaining).area / aoi_area > min_gain_frac and parsed:
        best_idx = -1
        best_new = 0.0
        best_cloud = float("inf")
        for idx, (product, footprint) in enumerate(parsed):
            new_area = _to_equal_area(remaining.intersection(footprint)).area
            if new_area <= 0:
                continue
            cloud = _cloud_cover(product)
            if (
                new_area > best_new
                or (new_area == best_new and cloud < best_cloud)
            ):
                best_idx = idx
                best_new = new_area
                best_cloud = cloud
        if best_idx < 0:
            break  # no candidate adds new coverage
        product, footprint = parsed.pop(best_idx)
        selected.append(product)
        remaining = remaining.difference(footprint)
    # Present in tile order for a stable, readable printout.
    return sorted(selected, key=lambda p: _tile_id(p) or "")


def confirm(prompt: str) -> bool:
    """Ask a yes/no question on stdin. Anything other than y/yes returns False."""
    try:
        reply = input(prompt).strip().lower()
    except EOFError:
        return False
    return reply in ("y", "yes")


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
def _open_download_stream(url: str, token: str) -> requests.Response:
    """Follow CDSE redirects manually, re-attaching the bearer token each hop.

    The requests library strips the ``Authorization`` header when a redirect
    crosses hosts (a security default). CDSE routinely redirects from
    ``catalogue.dataspace.copernicus.eu`` to ``download.dataspace.copernicus.eu``
    (and occasionally to a CreoDIAS S3 host), so we walk the redirect chain
    ourselves and only forward the token to trusted CDSE subdomains.
    """
    session = requests.Session()
    headers = {"Authorization": f"Bearer {token}"}
    max_hops = 10
    for _ in range(max_hops):
        response = session.get(
            url,
            headers=headers,
            stream=True,
            allow_redirects=False,
            timeout=600,
        )
        if response.status_code not in (301, 302, 303, 307, 308):
            return response
        next_url = response.headers.get("Location")
        response.close()
        if not next_url:
            raise requests.HTTPError(
                f"Redirect from {url} without a Location header",
                response=response,
            )
        url = requests.compat.urljoin(url, next_url)
        # Only forward credentials to CDSE-owned hosts.
        if ".dataspace.copernicus.eu" not in requests.utils.urlparse(url).netloc:
            headers = {}
    raise requests.HTTPError(f"Too many redirects while downloading {url}")


def product_local_path(product: dict, output_dir: Path) -> Path:
    """Return the on-disk path a product would be saved to."""
    filename = product["Name"].replace(".SAFE", "") + ".zip"
    return output_dir / filename


def product_already_downloaded(product: dict, output_dir: Path) -> bool:
    """True when the product's zip already exists on disk with non-zero size."""
    path = product_local_path(product, output_dir)
    return path.exists() and path.stat().st_size > 0


def existing_tiles_on_disk(output_dir: Path) -> dict[str, Path]:
    """Return a mapping of MGRS tile ID -> existing L1C zip on disk.

    Any Sentinel-2 L1C archive already present in ``output_dir`` counts as
    coverage for its MGRS tile, so a tile downloaded during a previous run
    (even with a different acquisition date) will not be re-downloaded.
    """
    tiles: dict[str, Path] = {}
    if not output_dir.exists():
        return tiles
    for path in sorted(output_dir.glob("S2*_MSIL1C_*.zip")):
        try:
            if path.stat().st_size == 0:
                continue
        except OSError:
            continue
        match = TILE_REGEX.search(path.name)
        if match is None:
            continue
        tiles.setdefault(match.group(1), path)
    return tiles


def download_product(product: dict, token: str, output_dir: Path) -> Path:
    """Stream a single product to disk. Skips files already downloaded."""
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = product_local_path(product, output_dir)
    if product_already_downloaded(product, output_dir):
        print(f"  Skipping (already downloaded): {destination.name}")
        return destination

    url = CDSE_DOWNLOAD_URL.format(pid=product["Id"])
    with _open_download_stream(url, token) as response:
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
    return destination


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _search_and_pick_per_tile(
    aoi_wkt: str,
    start_date: str,
    end_date: str,
) -> list[dict]:
    """Search the catalogue and return the least-cloudy scene per MGRS tile."""
    start_display = datetime.strptime(start_date, "%Y-%m-%d").strftime("%d/%m/%Y")
    end_display = datetime.strptime(end_date, "%Y-%m-%d").strftime("%d/%m/%Y")
    print(
        f"\nSearching Sentinel-2 {PRODUCT_TYPE} scenes between "
        f"{start_display} and {end_display} ..."
    )
    products = search_products(
        aoi_wkt=aoi_wkt,
        start_date=start_date,
        end_date=end_date,
        product_type=PRODUCT_TYPE,
        max_cloud_cover=MAX_CLOUD_COVER,
    )
    per_tile = select_least_cloudy_per_tile(products)
    print(
        f"  {len(products)} candidate scenes, "
        f"{len(per_tile)} unique MGRS tiles after least-cloudy per-tile filter"
    )
    return per_tile


def main() -> None:
    area_name, output_dir = prompt_output_folder(OUTPUT_PARENT_DIR)
    top_left, bottom_right = prompt_bbox()
    username, password = prompt_credentials()
    print(
        f"\nBuilding area of interest ({area_name}) from bounding box "
        f"upper-left={top_left}, lower-right={bottom_right} ..."
    )
    print(f"Downloads will go to {output_dir.resolve()}")
    aoi_geometry = load_bbox_geometry(top_left, bottom_right)
    aoi_wkt = build_search_wkt(aoi_geometry)

    per_tile = _search_and_pick_per_tile(aoi_wkt, START_DATE, END_DATE)
    coverage_pct = aoi_coverage_pct(aoi_geometry, per_tile)
    print(f"  Coverage of {area_name}: {coverage_pct:5.2f}%")

    if not per_tile:
        sys.exit("No candidate scenes found; nothing to do.")

    # ---- Minimum-scene set-cover ----------------------------------------
    print(
        "\nComputing minimum-scene cover (greedy set-cover, tie-break on "
        "cloud cover) ..."
    )
    min_cover = minimum_covering_scenes(aoi_geometry, per_tile)
    final_coverage = aoi_coverage_pct(aoi_geometry, min_cover)

    print(
        f"\nMinimum of {len(min_cover)} scene(s) needed to cover "
        f"{final_coverage:5.2f}% of {area_name} "
        f"(search window from {START_DATE} to {END_DATE}):"
    )
    for product in min_cover:
        cloud = _cloud_cover(product)
        marker = "  <-- cloudy" if cloud > CLOUD_WARN_PCT else ""
        print(
            f"    {_tile_id(product)}  "
            f"cloud={cloud:5.2f}%  "
            f"{product['Name']}{marker}"
        )

    cloudy = [p for p in min_cover if _cloud_cover(p) > CLOUD_WARN_PCT]
    if cloudy:
        print(
            f"\n  Note: {len(cloudy)} tile(s) have no scene cleaner than "
            f"{CLOUD_WARN_PCT:g}% cloud in this window. This is 'best "
            "available' rather than a cloud-free image."
        )

    print("\nDetailed coverage of the minimum set:")
    report_aoi_coverage(aoi_geometry, min_cover)

    # ---- Filter out anything already available on disk -------------------
    output_dir.mkdir(parents=True, exist_ok=True)
    existing_tiles = existing_tiles_on_disk(output_dir)
    to_download = [p for p in min_cover if _tile_id(p) not in existing_tiles]
    already_covered = [p for p in min_cover if _tile_id(p) in existing_tiles]

    if already_covered:
        print(
            f"\n{len(already_covered)} of {len(min_cover)} tile(s) already on "
            f"disk in {output_dir}; will not be re-downloaded:"
        )
        for product in already_covered:
            tile = _tile_id(product) or "?"
            existing_name = existing_tiles[tile].name
            proposed_name = product_local_path(product, output_dir).name
            if existing_name == proposed_name:
                print(f"    {tile}: {existing_name} (exact match)")
            else:
                print(
                    f"    {tile}: have {existing_name}, "
                    f"would-have downloaded {proposed_name}"
                )

    if not to_download:
        print("\nAll required tiles are already on disk. Nothing to download.")
        return

    total_bytes = sum(
        int(p.get("ContentLength", 0) or 0) for p in to_download
    )
    total_gb = total_bytes / (1024 ** 3) if total_bytes else None
    size_hint = f" (~{total_gb:.1f} GB)" if total_gb else ""

    print(f"\n{len(to_download)} new scene(s) still needed{size_hint}:")
    for product in to_download:
        print(
            f"    {_tile_id(product)}  "
            f"cloud={_cloud_cover(product):5.2f}%  "
            f"{product['Name']}"
        )

    if not confirm(
        f"\nProceed to download {len(to_download)} scene(s) to "
        f"{output_dir}? [y/N]: "
    ):
        print("Aborted by user.")
        return

    print(f"\nDownloading {len(to_download)} new scene(s) to {output_dir} ...")
    token = get_access_token(username, password)
    for product in to_download:
        try:
            download_product(product, token, output_dir)
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 401:
                # Token likely expired mid-run; refresh and retry once.
                token = get_access_token(username, password)
                download_product(product, token, output_dir)
            else:
                print(f"  Failed to download {product['Name']}: {exc}")

    print("Done.")


if __name__ == "__main__":
    main()
