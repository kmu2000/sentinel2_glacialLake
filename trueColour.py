"""Generate a true-colour RGB GeoTIFF for every Sentinel-2 L1C scene
below ``Data/Sentinel/``.

At start-up the script prompts for the state / AOI name; the resulting
TCI GeoTIFFs are written into ``Data/Sentinel/TCI/<state>/`` (the state
name is sanitised so it is always a safe folder name).

For each ``S2X_MSIL1C_*.zip`` archive the script:

1. Opens the Red (B04), Green (B03) and Blue (B02) bands - all native
   10 m - directly from the zip using GDAL's ``/vsizip/`` virtual
   filesystem.
2. Converts digital numbers to top-of-atmosphere reflectance using the
   per-band ``RADIO_ADD_OFFSET`` values declared in
   ``MTD_MSIL1C.xml`` (offsets were introduced in processing baseline
   04.00).
3. Applies a fixed reflectance -> 8-bit stretch:
   ``reflectance <= REF_MIN -> 0``, ``reflectance >= REF_MAX -> 255``.
   Optional gamma correction (``GAMMA > 1`` lifts shadows).
4. Writes a DEFLATE-compressed 3-band uint8 GeoTIFF with photometric
   interpretation RGB and an internal mask band so tile edges (L1C
   nodata) display as transparent in QGIS / rasterio.

Outputs
-------
* ``Data/Sentinel/TCI/<state>/DD-MM-YYYY_<TILE>_TCI.tif``

Requirements::

    pip install rasterio numpy tqdm
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Callable

import numpy as np
import rasterio
from rasterio.enums import ColorInterp
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INPUT_DIR = Path("../Data/Sentinel")
# TCI GeoTIFFs are written under ``TCI_PARENT_DIR / <state>/`` where
# ``<state>`` is the sanitised name entered at start-up (see
# ``prompt_state_folder`` below).
TCI_PARENT_DIR = INPUT_DIR / "TCI"

RED_BAND = "B04"
GREEN_BAND = "B03"
BLUE_BAND = "B02"

# Reflectance -> uint8 stretch. Values <= REF_MIN become 0, values >=
# REF_MAX become 255. REF_MAX = 0.30 is a good default for mixed
# land/water Sentinel-2 scenes; raise to ~0.4 for very bright scenes
# (fresh snow, deserts) or drop to ~0.20 to brighten dark scenes.
REF_MIN = 0.00
REF_MAX = 0.30
# Optional gamma. GAMMA > 1 brightens shadows (mid-tones lifted);
# GAMMA < 1 darkens them. Set to 1.0 for a strict linear stretch.
GAMMA = 1.2

QUANTIFICATION_VALUE = 10000.0  # DN -> reflectance scale factor
DEFAULT_OFFSET = 0.0            # for pre-baseline-04 products

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
# Stretch
# ---------------------------------------------------------------------------
def _dn_to_reflectance(dn: np.ndarray, offset: float) -> np.ndarray:
    return (dn.astype(np.float32) + offset) / QUANTIFICATION_VALUE


def _stretch_to_byte(reflectance: np.ndarray) -> np.ndarray:
    """Linear stretch of ``reflectance`` from [REF_MIN, REF_MAX] to [0, 255]
    with an optional gamma correction."""
    span = REF_MAX - REF_MIN
    if span <= 0:
        raise ValueError("REF_MAX must be strictly greater than REF_MIN")
    scaled = (reflectance - REF_MIN) / span
    np.clip(scaled, 0.0, 1.0, out=scaled)
    if GAMMA != 1.0:
        scaled = np.power(scaled, 1.0 / GAMMA, dtype=np.float32)
    return (scaled * 255.0).astype(np.uint8)


# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------
def process_scene(zip_path: Path, output_dir: Path) -> Path | None:
    """Build a true-colour RGB GeoTIFF for a single L1C zip archive."""
    match = NAME_REGEX.search(zip_path.name)
    if match is None:
        tqdm.write(f"  Skipping (name not recognised): {zip_path.name}")
        return None

    date_stamp = match.group("date")  # YYYYMMDD
    tile = match.group("tile")
    date_uk = f"{date_stamp[6:8]}-{date_stamp[4:6]}-{date_stamp[0:4]}"
    out_path = output_dir / f"{date_uk}_{tile}_TCI.tif"
    if out_path.exists():
        tqdm.write(f"  Skipping (exists): {out_path.name}")
        return out_path

    names = _list_zip(zip_path)
    product_xml = _find_member(names, lambda n: n.endswith("/MTD_MSIL1C.xml"))
    red_jp2 = _find_member(
        names,
        lambda n: "/IMG_DATA/" in n and n.endswith(f"_{RED_BAND}.jp2"),
    )
    green_jp2 = _find_member(
        names,
        lambda n: "/IMG_DATA/" in n and n.endswith(f"_{GREEN_BAND}.jp2"),
    )
    blue_jp2 = _find_member(
        names,
        lambda n: "/IMG_DATA/" in n and n.endswith(f"_{BLUE_BAND}.jp2"),
    )

    if not (red_jp2 and green_jp2 and blue_jp2):
        tqdm.write(f"  Missing B02/B03/B04 in {zip_path.name}; skipping.")
        return None

    offsets: dict[int, float] = {}
    if product_xml is not None:
        try:
            offsets = read_radiometric_offsets(zip_path, product_xml)
        except Exception as exc:  # noqa: BLE001
            tqdm.write(
                f"  Warning: could not read radiometric offsets for "
                f"{zip_path.name} ({exc}); assuming 0."
            )
    red_off = offsets.get(BAND_INDEX[RED_BAND], DEFAULT_OFFSET)
    green_off = offsets.get(BAND_INDEX[GREEN_BAND], DEFAULT_OFFSET)
    blue_off = offsets.get(BAND_INDEX[BLUE_BAND], DEFAULT_OFFSET)

    with rasterio.open(_vsizip(zip_path, red_jp2)) as src:
        red = src.read(1)
        target_transform = src.transform
        target_crs = src.crs
        height, width = red.shape
    with rasterio.open(_vsizip(zip_path, green_jp2)) as src:
        green = src.read(1)
    with rasterio.open(_vsizip(zip_path, blue_jp2)) as src:
        blue = src.read(1)

    nodata_mask = (red == 0) | (green == 0) | (blue == 0)

    r8 = _stretch_to_byte(_dn_to_reflectance(red, red_off))
    g8 = _stretch_to_byte(_dn_to_reflectance(green, green_off))
    b8 = _stretch_to_byte(_dn_to_reflectance(blue, blue_off))

    r8[nodata_mask] = 0
    g8[nodata_mask] = 0
    b8[nodata_mask] = 0

    profile = {
        "driver": "GTiff",
        "dtype": "uint8",
        "count": 3,
        "width": width,
        "height": height,
        "transform": target_transform,
        "crs": target_crs,
        "compress": "deflate",
        "predictor": 2,  # horizontal-differencing predictor for uint8
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "photometric": "RGB",
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.colorinterp = (
            ColorInterp.red,
            ColorInterp.green,
            ColorInterp.blue,
        )
        dst.write(r8, 1)
        dst.write(g8, 2)
        dst.write(b8, 3)
        # Internal mask band: 255 = valid, 0 = transparent (nodata).
        dst.write_mask(((~nodata_mask).astype(np.uint8)) * 255)
        dst.set_band_description(1, "Red (B04)")
        dst.set_band_description(2, "Green (B03)")
        dst.set_band_description(3, "Blue (B02)")

    pct_valid = float((~nodata_mask).sum()) / red.size * 100.0
    tqdm.write(f"  {out_path.name}: {pct_valid:5.2f}% valid pixels")
    return out_path


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
    state_name, output_dir = prompt_state_folder(TCI_PARENT_DIR)

    # rglob so we pick up scenes both at the root of ``Data/Sentinel/``
    # (older layout) and inside an AOI sub-folder created by
    # ``downloadSentinel2.py`` (e.g. ``Data/Sentinel/<state>/*.zip``).
    scenes = sorted(INPUT_DIR.rglob("S2*_MSIL1C_*.zip"))
    if not scenes:
        sys.exit(
            f"No Sentinel-2 L1C zip archives found under {INPUT_DIR.resolve()}"
        )

    print(f"\nFound {len(scenes)} scene(s) under {INPUT_DIR.resolve()}")
    print(f"True-colour outputs ({state_name}) -> {output_dir.resolve()}")
    print(
        f"Stretch: reflectance [{REF_MIN:.2f}, {REF_MAX:.2f}] -> [0, 255], "
        f"gamma = {GAMMA:g}"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    for scene in tqdm(scenes, desc="TCI", unit="scene"):
        try:
            process_scene(scene, output_dir)
        except Exception as exc:  # noqa: BLE001
            tqdm.write(f"  Failed on {scene.name}: {exc}")

    print("Done.")


if __name__ == "__main__":
    main()
