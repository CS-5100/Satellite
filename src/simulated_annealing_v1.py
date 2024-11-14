# Importing necessary liraries
import time
import math
import random
from random import randint
import matplotlib.pyplot as plt

from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd

from shapely.ops import unary_union
from shapely.geometry import Point
from shapely.affinity import translate

import tle_processing as tlp

# Constants
EARTH_SURFACE_AREA_SQ_KM = 509600000
EARTH_LAND_AREA_SQ_KM = 148326000
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Global Parameters
PLOT = True
BUFFER_RADIUS = 121065
PERTURB_DISTANCE_KM = 500
BUFFER_PLOTS = True

# Path to data
dirpath = Path(__file__).parent.resolve() / ".." / "data"

def load_map_data():
    # Load land and ocean map data (adjust paths as necessary)
    land_filepath = dirpath / "map_data" / "ne_10m_land_scale_rank.zip"
    ocean_filepath = dirpath / "map_data" / "ne_10m_ocean_scale_rank.zip"
    land_map = gpd.read_file(filename=land_filepath).to_crs(epsg=EQUAL_AREA_EPSG)
    ocean_map = gpd.read_file(filename=ocean_filepath).to_crs(epsg=EQUAL_AREA_EPSG)
    return land_map, ocean_map

def load_existing_satellites():
    # Download and convert TLE data to GeoDataFrame for existing satellites
    starlink_current_tle_list = tlp.download_current_tles_as_list()
    starlink_current = tlp.tles_to_dataframe(raw_tle_list=starlink_current_tle_list)
    starlink_gdf = tlp.tle_dataframe_to_geodataframe(starlink_current)
    starlink_gdf = starlink_gdf.to_crs(epsg=EQUAL_AREA_EPSG)
    print(starlink_gdf.head())
    return starlink_gdf

def generate_new_satellites(num_satellites=60):
    # Bounds to center satellites near common land areas (adjust as needed)
    lat_range = (-60, 70)  # Latitude range covering most inhabited regions
    lon_range = (-180, 180)  # Longitude range covering major land areas

    # Generate satellite data with coordinates closer to land
    new_satellite_data = {
        'Satellite Name': [f'NewSat{i+1}' for i in range(num_satellites)],
        'Latitude': [random.uniform(*lat_range) for _ in range(num_satellites)],
        'Longitude': [random.uniform(*lon_range) for _ in range(num_satellites)]
    }
    new_satellite_data['geometry'] = [Point(lon, lat) for lon, lat in zip(new_satellite_data['Longitude'], new_satellite_data['Latitude'])]

    # Construct the GeoDataFrame with geometry explicitly and re-project
    new_satellites_gdf = gpd.GeoDataFrame(new_satellite_data, geometry='geometry', crs="EPSG:4326").to_crs(epsg=EQUAL_AREA_EPSG)
    new_satellites_gdf['new_satellite'] = True  # Mark new satellites
    return new_satellites_gdf

def perturb_positions(gdf, new_satellite_column, max_shift_km=500, random_state=None):
    """
    Perturbs positions of a subset of new satellites randomly by up to max_shift_km.

    Args:
        gdf (GeoDataFrame): The GeoDataFrame containing satellite positions.
        new_satellite_column (str): Column name indicating new satellites.
        max_shift_km (float): Maximum shift in kilometers for each satellite.
        random_state (int): Random state for reproducibility.
    """
    # Define the max shift in meters
    max_shift_m = max_shift_km * 1000  # Convert km to meters

    # Select a subset of new satellites
    subset = gdf[gdf[new_satellite_column]].sample(frac=0.4, random_state=random_state)

    # Temporarily project to an equal-distance projection for accurate perturbation
    subset = subset.to_crs(epsg=EQUAL_DISTANCE_EPSG)

    # Apply random translations
    gdf.loc[subset.index, 'geometry'] = subset['geometry'].apply(
        lambda geom: translate(geom, xoff=np.random.uniform(-max_shift_m, max_shift_m),
                                      yoff=np.random.uniform(-max_shift_m, max_shift_m))
    )

    # Project back to the equal-area projection
    gdf = gdf.to_crs(epsg=EQUAL_AREA_EPSG)

def calculate_land_coverage(gdf, land_map):
    # Apply a buffer around each satellite to simulate coverage
    buffered_gdf = gdf.copy()
    buffered_gdf['geometry'] = buffered_gdf['geometry'].buffer(121065)  # Buffer of ~121.065 km radius

    # Find intersection area between satellite coverage and land area
    intersections = buffered_gdf.overlay(land_map, how='intersection')
    total_land_coverage = intersections['geometry'].area.sum() / (1000**2)  # Sum of non-overlapping areas in km²
    print(f"Current land coverage: {total_land_coverage:.2f} km²")  # Track coverage each iteration
    return total_land_coverage

# Load existing satellite data from TLEs and flag new satellites
existing_satellites_gdf = load_existing_satellites()
existing_satellites_gdf['new_satellite'] = False  # Mark existing satellites

# Generate new satellites and combine with existing ones
new_satellites_gdf = generate_new_satellites(num_satellites=60)
satellite_gdf = pd.concat([existing_satellites_gdf, new_satellites_gdf])

# Load map data
land_map, ocean_map = load_map_data()

def simulated_annealing(satellite_gdf, land_map, new_satellite_column, num_iterations=100,
                        initial_temp=1000, cooling_rate=0.95, restart_threshold=10,
                        max_steps_without_improvement=15, energy_difference_threshold=0.1):
    best_gdf = satellite_gdf.copy()
    best_land_coverage = calculate_land_coverage(best_gdf, land_map)
    current_gdf = best_gdf.copy()
    current_land_coverage = best_land_coverage
    temperature = initial_temp

    # Tracking iteration times and land coverage for plotting
    iteration_times = []
    best_land_coverage_list = []

    # List for storing improved solutions
    improvement_solutions = [(best_land_coverage, best_gdf.copy())]

    # Counters for restart criteria
    steps_without_improvement = 0

    for iteration in range(num_iterations):
        start_time = time.time()  # Start timing the iteration
        print(f"Iteration {iteration + 1}/{num_iterations}, Temperature: {temperature:.2f}")

        # Perturb positions of new satellites only
        new_gdf = current_gdf.copy()
        perturb_positions(new_gdf, new_satellite_column, random_state=random.randint(0, 10000))

        # Calculate new land coverage area
        new_land_coverage = calculate_land_coverage(new_gdf, land_map)

        # Calculate the improvement or deterioration
        delta_coverage = new_land_coverage - current_land_coverage

        # Acceptance probability for worse solutions
        acceptance_probability = math.exp(delta_coverage / temperature) if delta_coverage < 0 else 1

        if delta_coverage > 0 or random.uniform(0, 1) < acceptance_probability:
            print(f"Accepted new configuration with land coverage: {new_land_coverage:.2f} km²")
            current_gdf, current_land_coverage = new_gdf, new_land_coverage
            steps_without_improvement = 0  # Reset steps without improvement

            # Update the best solution if this configuration is better
            if new_land_coverage > best_land_coverage:
                best_gdf, best_land_coverage = new_gdf, new_land_coverage

            # Store new solution if it's an improvement and unique
            if not any(abs(new_land_coverage - cov) < 1e-6 for cov, _ in improvement_solutions):
                improvement_solutions.append((new_land_coverage, new_gdf.copy()))
        else:
            print("Retained current configuration.")
            steps_without_improvement += 1  # Increment if no improvement

        # Cooling step
        temperature *= cooling_rate

        # Restart criteria: if no improvement over max_steps_without_improvement or significant drop in energy
        if steps_without_improvement >= max_steps_without_improvement or \
           (best_land_coverage - current_land_coverage) / best_land_coverage > energy_difference_threshold:
            print("Restarting due to lack of improvement or significant drop in coverage.")

            # Randomly select an improved solution to restart from
            if improvement_solutions:
                selected_coverage, selected_gdf = random.choice(improvement_solutions)
                current_gdf = selected_gdf.copy()  # Restart from the chosen solution
                current_land_coverage = selected_coverage
            else:
                current_gdf = best_gdf.copy()  # Default to the best solution if no other improvements exist

            temperature = initial_temp  # Reset temperature
            steps_without_improvement = 0  # Reset the stagnation counter

        # Track iteration times and best land coverage
        elapsed_time = time.time() - start_time
        iteration_times.append(elapsed_time)
        best_land_coverage_list.append(best_land_coverage)

    print(f"Optimized land coverage area after simulated annealing: {best_land_coverage:.2f} km²")
    return best_gdf, best_land_coverage, iteration_times, best_land_coverage_list

# Run the simulated annealing optimization with random restarts
optimized_gdf, optimized_land_coverage, iteration_times, best_land_coverage_list = simulated_annealing(
    satellite_gdf, land_map, 'new_satellite', num_iterations=100
)

# Plot best land coverage over iterations
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(best_land_coverage_list) + 1), best_land_coverage_list, marker='o', color='g')
plt.xlabel("Iteration Number")
plt.ylabel("Best Land Coverage (km²)")
plt.title("Best Land Coverage Over Iterations with Simulated Annealing")
plt.grid(True)
plt.show()

# Plot number of iterations vs time
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(iteration_times) + 1), iteration_times, marker='o', color='b')
plt.xlabel("Iteration Number")
plt.ylabel("Time (seconds)")
plt.title("Iteration vs. Time Taken for Local Search Optimization")
plt.grid(True)
plt.show()

# Define the output file path
output_filepath = Path.cwd() / "data/optimized_satellite_positions.txt"

# Open the file in write mode and save the results
with open(output_filepath, "w") as file:
    # Write the optimized land coverage area
    file.write(f"Optimized land coverage area: {optimized_land_coverage:.2f} km²\n\n")

    # Write the header for new satellite positions
    file.write("New Satellites' Final Positions After Optimization:\n")

    # Write each new satellite's final position
    for _, row in optimized_gdf[optimized_gdf['new_satellite']].iterrows():
        file.write(f"{row['Satellite Name']}: Latitude {row['Latitude']:.2f}, "
                   f"Longitude {row['Longitude']:.2f}, Geometry {row['geometry']}\n")

print(f"Output saved to {output_filepath}")