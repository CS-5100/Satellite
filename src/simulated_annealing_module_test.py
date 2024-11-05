# Importing necessary liraries
import time
import math
import random
from random import randint
import matplotlib.pyplot as plt

from pathlib import Path
import pandas as pd

import geodataframe_processing as gdfp
from custom_plotting import plot_initial_final_satellites

# Constants
EARTH_SURFACE_AREA_SQ_KM = 509600000
EARTH_LAND_AREA_SQ_KM = 148326000
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Global Parameters
BUFFER_RADIUS = 121065
PERTURB_DISTANCE_KM = 500

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

def simulated_annealing(satellite_gdf, map, new_satellite_column, buffer_radius,
                        num_iterations=100,
                        initial_temp=1000, cooling_rate=0.95, restart_threshold=10):
    best_gdf = satellite_gdf.copy()
    best_land_coverage = gdfp.calculate_land_coverage(best_gdf, map, buffer_radius)
    current_gdf = best_gdf.copy()
    current_land_coverage = best_land_coverage
    temperature = initial_temp

    # Tracking iteration times and land coverage for plotting
    iteration_times = []
    best_land_coverage_list = []

    # Restart counter
    stagnation_counter = 0

    for iteration in range(num_iterations):
        start_time = time.time()  # Start timing the iteration
        print(f"Iteration {iteration + 1}/{num_iterations}, Temperature: {temperature:.2f}")

        # Perturb positions of new satellites only
        new_gdf = current_gdf.copy()
        gdfp.perturb_positions(
            new_gdf,
            new_satellite_column,
            eq_dist_epsg=EQUAL_DISTANCE_EPSG,
            eq_area_epsg=EQUAL_AREA_EPSG,
            max_shift_km=PERTURB_DISTANCE_KM,
            random_state=randint(0, 10000),
        )

        # Calculate new land coverage area
        new_land_coverage = new_land_coverage = gdfp.calculate_land_coverage(
            new_gdf, map, buffer_radius
        )

        # Calculate the improvement or deterioration
        delta_coverage = new_land_coverage - current_land_coverage

        # Acceptance probability for worse solutions
        acceptance_probability = math.exp(delta_coverage / temperature) if delta_coverage < 0 else 1

        if delta_coverage > 0 or random.uniform(0, 1) < acceptance_probability:
            print(f"Accepted new configuration with land coverage: {new_land_coverage:.2f} km²")
            current_gdf, current_land_coverage = new_gdf, new_land_coverage
            stagnation_counter = 0  # Reset stagnation counter
            # Update best solution if this is the best encountered so far
            if new_land_coverage > best_land_coverage:
                best_gdf, best_land_coverage = new_gdf, new_land_coverage
        else:
            print("Retained current configuration.")

        # Cooling step
        temperature *= cooling_rate

        # Increment stagnation counter if no improvement
        stagnation_counter += 1

        # Check if we should restart due to stagnation
        if stagnation_counter >= restart_threshold:
            print("Random restart due to stagnation.")
            current_gdf = satellite_gdf.copy()
            current_land_coverage = gdfp.calculate_land_coverage(current_gdf, map, buffer_radius)
            stagnation_counter = 0  # Reset the stagnation counter
            temperature = initial_temp  # Reset temperature for new start

        # Track iteration times and best land coverage
        elapsed_time = time.time() - start_time
        iteration_times.append(elapsed_time)
        best_land_coverage_list.append(best_land_coverage)

    print(f"Optimized land coverage area after simulated annealing: {best_land_coverage:.2f} km²")
    return best_gdf, best_land_coverage, iteration_times, best_land_coverage_list

# Run the simulated annealing optimization with random restarts
optimized_gdf, optimized_land_coverage, iteration_times, best_land_coverage_list = simulated_annealing(
    satellite_gdf, land_map, 'new_satellite', BUFFER_RADIUS, num_iterations=50
)

# extracting the new satellites from the GeoDataFrame for plotting purposes
final_satellites_gdf = optimized_gdf[optimized_gdf["new_satellite"]].copy()

# # Plot best land coverage over iterations
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(best_land_coverage_list) + 1), best_land_coverage_list, marker='o', color='g')
# plt.xlabel("Iteration Number")
# plt.ylabel("Best Land Coverage (km²)")
# plt.title("Best Land Coverage Over Iterations with Simulated Annealing")
# plt.grid(True)
# plt.show()

# # Plot number of iterations vs time
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(iteration_times) + 1), iteration_times, marker='o', color='b')
# plt.xlabel("Iteration Number")
# plt.ylabel("Time (seconds)")
# plt.title("Iteration vs. Time Taken for Local Search Optimization")
# plt.grid(True)
# plt.show()

plot_initial_final_satellites(existing_satellites=existing_satellites_gdf,
                                  initial_positions=initial_satellites_gdf,
                                  final_positions=final_satellites_gdf,
                                  land=land_map,
                                  ocean=ocean_map,
                                  buffer=None)

# # Define the output file path
# output_filepath = Path.cwd() / "data/optimized_satellite_positions.txt"

# # Open the file in write mode and save the results
# with open(output_filepath, "w") as file:
#     # Write the optimized land coverage area
#     file.write(f"Optimized land coverage area: {optimized_land_coverage:.2f} km²\n\n")

#     # Write the header for new satellite positions
#     file.write("New Satellites' Final Positions After Optimization:\n")

#     # Write each new satellite's final position
#     for _, row in optimized_gdf[optimized_gdf['new_satellite']].iterrows():
#         file.write(f"{row['Satellite Name']}: Latitude {row['Latitude']:.2f}, "
#                    f"Longitude {row['Longitude']:.2f}, Geometry {row['geometry']}\n")

# print(f"Output saved to {output_filepath}")