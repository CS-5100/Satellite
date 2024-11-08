import pandas as pd
import pyproj
import geodataframe_processing as gdfp
from search_functions import hill_climbing
from custom_plotting import plot_initial_final_satellites

# Global Parameters
BUFFER_RADIUS = 12065
PERTURB_DISTANCE_KM = (
    500  # we may want to make this something that decays exponentially
)
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# initial check what the bounds are for the given projection
eq_area = pyproj.CRS.from_epsg(EQUAL_AREA_EPSG)
transformer = pyproj.Transformer.from_crs(eq_area.geodetic_crs, eq_area, always_xy=True)
print(transformer.transform_bounds(*eq_area.area_of_use.bounds))

# Load existing satellite data from TLEs and flag new satellites
existing_satellites_gdf = gdfp.load_existing_satellites(EPSG=EQUAL_AREA_EPSG, show_head=True)
existing_satellites_gdf["new_satellite"] = False  # Mark existing satellites

# Load map data
land_map, ocean_map = gdfp.load_land_ocean_data(EPSG=EQUAL_AREA_EPSG)

# Generate new satellites, create a copy of them for plotting, and combine with existing ones
new_satellites_gdf = gdfp.generate_new_satellites(
    EPSG=EQUAL_AREA_EPSG,
    num_satellites=30,
    input_map=land_map,
    true_random=False,
    sample_separate=False,
)

initial_satellites_gdf = new_satellites_gdf.copy()
satellite_gdf = pd.concat([existing_satellites_gdf, new_satellites_gdf])

# Run local search with perturbation targeted to new satellites
optimized_gdf, optimized_land_coverage = hill_climbing(
    satellite_gdf=satellite_gdf,
    map=land_map,
    new_satellite_column="new_satellite",
    buffer_radius=BUFFER_RADIUS,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG,
    equal_area_epsg=EQUAL_AREA_EPSG,
    perturbation_distance=PERTURB_DISTANCE_KM,
    num_iterations=30,
)

# extracting the new satellites from the GeoDataFrame for plotting purposes
final_satellites_gdf = optimized_gdf[optimized_gdf["new_satellite"]].copy()

# re-project data onto a geodetic reference system for printing
# and updating coordinates
optimized_gdf_geodetic = optimized_gdf.to_crs(epsg=4326)
optimized_gdf_geodetic_coordinates = optimized_gdf_geodetic.get_coordinates()
optimized_gdf_geodetic["Latitude"] = optimized_gdf_geodetic_coordinates.y
optimized_gdf_geodetic["Longitude"] = optimized_gdf_geodetic_coordinates.x

# Print out details of the newly generated satellites
print("\nNew Satellites' Final Positions After Optimization:")
print(
    optimized_gdf_geodetic[optimized_gdf_geodetic["new_satellite"]][
        ["Satellite", "Longitude", "Latitude", "geometry"]
    ]
)

plot_initial_final_satellites(
    existing_satellites=existing_satellites_gdf,
    initial_positions=initial_satellites_gdf,
    final_positions=final_satellites_gdf,
    land=land_map,
    ocean=ocean_map,
)
