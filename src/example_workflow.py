import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from geodataframe_processing import (
    load_satellites_from_file,
    load_land_ocean_data,
    generate_new_satellites,
)
from search_functions import random_restart_simulated_annealing
from custom_plotting import plot_initial_final_satellites

# Global Parameter Definition
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933
BASE_BUFFER = 12065
TLE_TIME = datetime(2024, 11, 14, 16, 12, 32, 202983)
PERTURB_DISTANCE_KM = 100
LOCAL_SEARCH_ITERATIONS = 200


# Importing a set of TLEs as the existing satellites
# and copy the GeoDataFrame for plotting positions
base_gdf = load_satellites_from_file(
    EPSG=EQUAL_AREA_EPSG,
    input_filename="parameter_sensitivity_analysis_tles.txt",
    time=TLE_TIME,
)
existing_gdf = base_gdf.copy()

# Import the land map data
land_map, ocean_map = load_land_ocean_data(EPSG=EQUAL_AREA_EPSG)

# Generate an initial set of satellite candidates
# and make a copy for plotting
initial_candidate_gdf = generate_new_satellites(
    EPSG=EQUAL_AREA_EPSG,
    num_satellites=30,
    input_map=land_map,
    true_random=False,
    sample_separate=False,
)
initial_gdf = initial_candidate_gdf.copy()

# Concatenate a GeoDataFrame of existing satellites and candidate satellites
satellite_gdf = pd.concat([base_gdf, initial_candidate_gdf])

# Use one of the two search functions to optimize the land coverage area
# and return the optimal satellites
optimized_gdf, _, _, _ = random_restart_simulated_annealing(
    satellite_gdf=satellite_gdf,
    map=land_map,
    new_satellite_column="new_satellite",
    buffer_radius=BASE_BUFFER,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG,
    equal_area_epsg=EQUAL_AREA_EPSG,
    perturbation_distance=PERTURB_DISTANCE_KM,
    num_iterations=100,
)

# get the optimized candidate satellite positions
final_gdf = optimized_gdf[optimized_gdf["new_satellite"] == True]

# make and show the before and after plot
before_after_plot = plot_initial_final_satellites(existing_satellites=existing_gdf,
                                                  initial_positions=initial_gdf,
                                                  final_positions=final_gdf,
                                                  land=land_map,
                                                  ocean=ocean_map)
plt.show()