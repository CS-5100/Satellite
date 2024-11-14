import pandas as pd
import matplotlib.pyplot as plt
import geodataframe_processing as gdfp
from search_functions import random_restart_simulated_annealing
from custom_plotting import plot_initial_final_satellites

# Constants
EARTH_SURFACE_AREA_SQ_KM = 509600000
EARTH_LAND_AREA_SQ_KM = 148326000
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Global Parameters
BUFFER_RADIUS = 12065 # 121065
PERTURB_DISTANCE_KM = 50

# Load existing satellite data from TLEs and flag new satellites
existing_satellites_gdf = gdfp.load_existing_satellites(EPSG=EQUAL_AREA_EPSG)
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

# Run the simulated annealing optimization with random restarts
optimized_gdf, optimized_land_coverage, iteration_times, best_land_coverage_list = (
    random_restart_simulated_annealing(
        satellite_gdf=satellite_gdf,
        map=land_map,
        new_satellite_column="new_satellite",
        buffer_radius=BUFFER_RADIUS,
        equal_distance_epsg=EQUAL_DISTANCE_EPSG,
        equal_area_epsg=EQUAL_AREA_EPSG,
        perturbation_distance=PERTURB_DISTANCE_KM,
        num_iterations=100,
    )
)

# extracting the new satellites from the GeoDataFrame for plotting purposes
final_satellites_gdf = optimized_gdf[optimized_gdf["new_satellite"]].copy()

# Plot best land coverage over iterations
plt.figure(figsize=(10, 6))
plt.plot(
    range(1, len(best_land_coverage_list) + 1),
    best_land_coverage_list,
    marker="o",
    color="g",
)
plt.xlabel("Iteration Number")
plt.ylabel("Best Land Coverage (km²)")
plt.title("Best Land Coverage Over Iterations with Simulated Annealing")
plt.grid(True)
plt.show()

# Plot number of iterations vs time
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(iteration_times) + 1), iteration_times, marker="o", color="b")
plt.xlabel("Iteration Number")
plt.ylabel("Time (seconds)")
plt.title("Iteration vs. Time Taken for Local Search Optimization")
plt.grid(True)
plt.show()

plot_initial_final_satellites(
    existing_satellites=existing_satellites_gdf,
    initial_positions=initial_satellites_gdf,
    final_positions=final_satellites_gdf,
    land=land_map,
    ocean=ocean_map,
    buffer=None,
)