# sentinel2_glacialLake

Pipeline for mapping small glacial lakes from Sentinel-2 imagery. Given a
user-supplied bounding box (typically the extent of an Indian state), the
scripts download cloud-free Sentinel-2 L1C scenes and a matching
Copernicus DEM mosaic, generate cloud-masked NDWI rasters and true-colour
composites, extract small water bodies as polygons filtered by terrain
slope, and compare the result against the Zhang et al. (2022) reference
lake inventory.

## Pipeline overview

The scripts are designed to be run in order. Each writes into a
predictable folder under `../Data/` (see [Expected directory layout](#expected-directory-layout))
so the outputs of one step are the inputs of the next.

| Step | Script                | Reads from                          | Writes to                                    |
|-----:|-----------------------|-------------------------------------|----------------------------------------------|
| 1    | `downloadSentinel2.py`| CDSE catalogue                      | `../Data/Sentinel/<AOI>/*.zip`               |
| 2    | `downloadDEM.py`      | AWS Copernicus DEM open data        | `../Data/DEM/*.tif`                          |
| 3    | `trueColour.py`       | `../Data/Sentinel/**/*.zip`         | `../Data/Sentinel/TCI/<state>/*.tif`         |
| 4    | `cloudMaskNDWI.py`    | `../Data/Sentinel/**/*.zip`         | `../Data/Sentinel/NDWI/*.tif`                |
| 5    | `extractSmallLakes.py`| `../Data/Sentinel/NDWI/`, `../Data/DEM/` | `../Data/Lakes/*.gpkg`, `all_small_lakes.gpkg` |
| 6    | `PGDL.py`             | glacier + reference-lake shapefiles, `../Data/Lakes/all_small_lakes.gpkg` | `../Himachal_Pradesh_lakes/*.gpkg` |

### What each script does

- **`downloadSentinel2.py`** &mdash; Prompts for a state / AOI name, a
  bounding box (upper-left and lower-right lon/lat) and Copernicus Data
  Space credentials. Queries the CDSE OData catalogue, picks the
  least-cloudy Sentinel-2 L1C scene per MGRS tile, reduces that set to a
  minimum cover via greedy set-cover, reports coverage and cloud stats,
  and downloads any tiles not already on disk to
  `../Data/Sentinel/<AOI>/`.
- **`downloadDEM.py`** &mdash; Downloads 30 m Copernicus DEM (GLO-30) tiles
  covering the AOI from the AWS Open Data mirror into `../Data/DEM/`.
- **`trueColour.py`** &mdash; Prompts for the state / AOI name, then
  builds a stretched 8-bit RGB GeoTIFF (B04/B03/B02) per scene into
  `../Data/Sentinel/TCI/<state>/`. Scenes are discovered recursively
  under `../Data/Sentinel/`, so both the legacy flat layout and the AOI
  sub-folders created by `downloadSentinel2.py` are picked up.
- **`cloudMaskNDWI.py`** &mdash; Reads B03 (Green) and B08 (NIR) from each
  L1C zip, computes McFeeters NDWI, masks cloud + cirrus pixels using the
  L1C `MSK_CLASSI` layer (writes NaN), and optionally reports
  NDWI-in-lakes statistics against the Zhang et al. (2022) reference
  polygons.
- **`extractSmallLakes.py`** &mdash; Prompts for the lower and upper
  lake-area bounds (km²) and the maximum mean terrain slope (degrees),
  then thresholds each NDWI raster, does 8-connected component
  labelling, filters components by that area window and by mean slope
  computed from the Copernicus DEM, polygonises the survivors, and
  writes per-scene GeoPackages plus a combined `all_small_lakes.gpkg`
  in EPSG:4326. Defaults are `[0.001, 2.0] km²` and `< 20°`; press
  Enter at any prompt to keep the default.
- **`PGDL.py`** &mdash; Restricts the RGI glacier inventory, the Zhang
  reference lakes, and the computed lakes to a chosen state polygon, and
  further filters the computed lakes to those within 10 km of any
  glacier, writing per-dataset GeoPackages.

## Setup

### Prerequisites

- Python 3.10 or later.
- A [Copernicus Data Space](https://dataspace.copernicus.eu/) account
  (free) &mdash; needed by `downloadSentinel2.py`.
- The `libgdal-jp2openjpeg` GDAL driver, required by rasterio to open the
  Sentinel-2 JPEG-2000 bands directly from the zip archives. If you use
  conda:
  ```bash
  conda install -c conda-forge libgdal-jp2openjpeg
  ```

### Install the Python dependencies

```bash
python -m venv .venv
source .venv/bin/activate         # macOS / Linux
pip install -r requirements.txt
```

## Expected directory layout

The scripts use relative paths of the form `../Data/...`, so they expect
to be run from a `script/` folder sitting next to a `Data/` folder:

```
LakeMapping/
├── script/                      <-- this repository lives here
│   ├── downloadSentinel2.py
│   ├── downloadDEM.py
│   ├── trueColour.py
│   ├── cloudMaskNDWI.py
│   ├── extractSmallLakes.py
│   ├── PGDL.py
│   ├── requirements.txt
│   ├── README.md
│   └── data/                    (bundled reference shapefiles, in-repo)
│       ├── Zhang2022_Data/GL/Himalayan_Glacial_Lakes_2020.shp
│       └── gadm/gadm41_IND_1.shp    (GADM level-1 admin boundaries)
├── Data/                        (populated by the pipeline; user-owned)
│   ├── Sentinel/<AOI>/*.zip     (populated by downloadSentinel2.py)
│   ├── Sentinel/TCI/<state>/*.tif  (populated by trueColour.py)
│   ├── Sentinel/NDWI/*.tif      (populated by cloudMaskNDWI.py)
│   ├── DEM/*.tif                (populated by downloadDEM.py)
│   └── Lakes/*.gpkg             (populated by extractSmallLakes.py)
├── GlacierInventory/
│   └── 14_rgi60_SouthAsiaWest_India.shp
└── Himachal_Pradesh_lakes/      (populated by PGDL.py)
```

The Zhang 2022 glacial-lake inventory and the GADM state boundaries are
bundled inside `script/data/` for convenience. The `Data/`,
`GlacierInventory/` and `Himachal_Pradesh_lakes/` folders sitting
alongside `script/` are intentionally **not** part of this repository
&mdash; they are supplied by the user (RGI glaciers, downloaded
imagery) or produced by the pipeline.

> **Note.** The scripts currently read the reference shapefiles from
> `../Data/Zhang2022_Data/GL/...` and the GADM shapefile from a similar
> `../Data/gadm.../` location. If you cloned this repo fresh and want to
> use the bundled copies in `script/data/`, either update the
> `LAKES_SHP` / `GADM_SHP` paths in `cloudMaskNDWI.py` and `PGDL.py`, or
> symlink the folders into `../Data/` alongside the pipeline outputs.

### Reference data sources

- GADM administrative boundaries: <https://gadm.org/download_country.html>
- RGI 6.0 glacier outlines: <https://www.glims.org/RGI/rgi60_dl.html>
  (this project uses region `14_rgi60_SouthAsiaWest`).
- Zhang et al. (2022) Himalayan glacial-lake inventory: [Zenodo record](https://doi.org/10.5281/zenodo.5808110).
- Sentinel-2: [Copernicus Data Space](https://dataspace.copernicus.eu/).
- Copernicus DEM (GLO-30): [AWS Open Data](https://registry.opendata.aws/copernicus-dem/).

## Usage

Run the scripts in order from the `script/` folder. Every step is
idempotent &mdash; already-downloaded scenes, generated NDWI/TCI rasters,
and extracted lakes are detected on disk and skipped.

Three of the scripts are interactive:

- `downloadSentinel2.py` prompts for an AOI name (used as the sub-folder
  under `../Data/Sentinel/`), the upper-left and lower-right
  longitude/latitude of the bounding box, and your Copernicus Data
  Space username + password.
- `trueColour.py` prompts for a state / AOI name and writes the RGB
  GeoTIFFs into `../Data/Sentinel/TCI/<state>/`.
- `extractSmallLakes.py` prompts for the minimum and maximum lake area
  (km²) and the maximum mean slope (degrees). Press Enter at any prompt
  to accept the shown default.

```bash
cd script/

# 1. Download Sentinel-2 tiles for your AOI  (interactive)
python downloadSentinel2.py

# 2. Download the matching Copernicus DEM tiles
python downloadDEM.py

# 3. (Optional) build RGB composites for visual QC  (interactive)
python trueColour.py

# 4. Compute cloud-masked NDWI rasters
python cloudMaskNDWI.py

# 5. Extract small lakes (with slope filter)  (interactive)
python extractSmallLakes.py

# 6. Compare to the reference glacier + lake inventories
python PGDL.py
```

## Tuning knobs

Each script keeps its knobs in a short config block near the top:

- `downloadSentinel2.py`: `START_DATE`, `END_DATE`, `MAX_CLOUD_COVER`,
  `CLOUD_WARN_PCT`. (Bounding box, AOI name and credentials are asked
  for at runtime.)
- `cloudMaskNDWI.py`: `LAKES_SHP` (Zhang polygons used for the
  NDWI-in-lakes statistics).
- `extractSmallLakes.py`: `NDWI_THRESHOLD` (default 0.20). The area
  window and slope cap are asked for at runtime with defaults
  `DEFAULT_MIN_LAKE_AREA_KM2` = 0.001 km², `DEFAULT_MAX_LAKE_AREA_KM2`
  = 2.0 km² and `DEFAULT_SLOPE_MAX_DEG` = 20°; edit those constants to
  change what pressing Enter accepts.
- `trueColour.py`: `REF_MIN`, `REF_MAX`, `GAMMA` for the display stretch.
  (The state / AOI sub-folder is asked for at runtime.)

## Credits

- Sentinel-2 L1C data: Copernicus Sentinel data (ESA), retrieved via the
  Copernicus Data Space Ecosystem.
- Copernicus DEM (GLO-30): produced by Airbus, distributed via the
  [Copernicus Programme](https://spacedata.copernicus.eu/).
- Zhang, G. *et al.* (2022) glacial lake inventory.
- Randolph Glacier Inventory 6.0.
