import os
import geopandas as gpd
from shapely.ops import unary_union

outdir = "../Himachal_Pradesh_lakes"
os.makedirs(outdir, exist_ok=True)
glacier_file = "../GlacierInventory/14_rgi60_SouthAsiaWest_India.shp"
state_file = "../Data/gadm41_IND_shp/gadm41_IND_1.shp"
state_name = "Himachal Pradesh"

reference_lake_file = "../Data/Zhang2022_Data/GL/Himalayan_Glacial_Lakes_2020.shp"
lake_file = "../Data/Lakes/all_small_lakes.gpkg"
buffer = 10000  # metres

glacier_gdf = gpd.read_file(glacier_file)
state_gdf = gpd.read_file(state_file)
reference_lake_gdf = gpd.read_file(reference_lake_file)
lake_gdf = gpd.read_file(lake_file)

# Extract the state polygon in WGS84 (GADM level-1 uses the NAME_1 field).
state_subset = state_gdf[
    state_gdf["NAME_1"].str.casefold() == state_name.casefold()
]
if state_subset.empty:
    raise ValueError(f"State {state_name!r} not found in {state_file}")
state_subset = state_subset.to_crs("EPSG:4326")
state_geom = unary_union(list(state_subset.geometry))

# Subset each layer to features intersecting the state polygon; reproject
# to a common CRS first so the spatial predicate is well defined.
glacier_state = glacier_gdf.to_crs("EPSG:4326")
glacier_state = glacier_state[glacier_state.geometry.intersects(state_geom)].copy()

lake_state_ref = reference_lake_gdf.to_crs("EPSG:4326")
lake_state_ref = lake_state_ref[lake_state_ref.geometry.intersects(state_geom)].copy()

lake_state_calc = lake_gdf.to_crs("EPSG:4326")
lake_state_calc = lake_state_calc[lake_state_calc.geometry.intersects(state_geom)].copy()

print(f"Glaciers in {state_name}:       {len(glacier_state):5d}")
print(f"Reference lakes in {state_name}: {len(lake_state_ref):5d}")
print(f"Computed lakes in {state_name}:  {len(lake_state_calc):5d}")

# Keep only lakes lying within `buffer` metres of any glacier. Reproject
# to a local UTM CRS so `buffer` (in metres) is meaningful, then use a
# spatial join against the buffered glaciers (rtree spatial index makes
# this fast even with thousands of glaciers).
utm_crs = glacier_state.estimate_utm_crs()
glacier_utm = glacier_state.to_crs(utm_crs)
lake_utm = lake_state_calc.to_crs(utm_crs)

buffered = glacier_utm[["geometry"]].copy()
buffered["geometry"] = glacier_utm.buffer(buffer)

matched = gpd.sjoin(lake_utm, buffered, predicate="intersects", how="inner")
lake_state_filtered = lake_state_calc.loc[matched.index.unique()].copy()

print(
    f"Lakes within {buffer / 1000:g} km of any glacier: "
    f"{len(lake_state_filtered):5d}"
)

# Write the state-clipped subsets to the output directory. "fid" is a
# reserved column name in GeoPackage (used for the primary key), so it
# must be dropped before writing if it slipped in from an input .gpkg.
def _to_gpkg(gdf, path, layer):
    gdf.drop(columns=["fid"], errors="ignore").to_file(
        path, driver="GPKG", layer=layer
    )


glacier_out = os.path.join(outdir, "glacier_state.gpkg")
lake_ref_out = os.path.join(outdir, "lake_state_ref.gpkg")
lake_calc_out = os.path.join(outdir, "lake_state_calc.gpkg")
lake_filtered_out = os.path.join(outdir, "lake_state_filtered.gpkg")

_to_gpkg(glacier_state, glacier_out, "glacier_state")
_to_gpkg(lake_state_ref, lake_ref_out, "lake_state_ref")
_to_gpkg(lake_state_calc, lake_calc_out, "lake_state_calc")
_to_gpkg(lake_state_filtered, lake_filtered_out, "lake_state_filtered")

print(f"\nWrote outputs to {os.path.abspath(outdir)}:")
print(f"  {os.path.basename(glacier_out)}")
print(f"  {os.path.basename(lake_ref_out)}")
print(f"  {os.path.basename(lake_calc_out)}")
print(f"  {os.path.basename(lake_filtered_out)}")
