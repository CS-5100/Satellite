import pandas as pd
import numpy as np
from pathlib import Path
import geopandas as gpd
from shapely.geometry import Point
import tle_processing as tlp 
from random import uniform
from search_functions import random_restart_simulated_annealing  
from geodataframe_processing import calculate_land_coverage  

# Constants
BUFFER_RADIUS = 121065  # Approximate coverage radius of a satellite (in meters)
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Helper Functions
def load_map_data():
    """
    Load land and ocean map data, reprojecting to an equal-area EPSG code.
    """
    dirpath = Path(__file__).parent.resolve() / ".." / "data"
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

def generate_new_satellites(num_satellites=30):
    """
    Generate a new set of satellites with random positions.
    """
    lat_range = (-60, 70)
    lon_range = (-180, 180)
    new_satellite_data = {
        'Satellite Name': [f'NewSat{i+1}' for i in range(num_satellites)],
        'Latitude': [uniform(*lat_range) for _ in range(num_satellites)],
        'Longitude': [uniform(*lon_range) for _ in range(num_satellites)]
    }
    new_satellite_data['geometry'] = [Point(lon, lat) for lon, lat in zip(new_satellite_data['Longitude'], new_satellite_data['Latitude'])]
    new_satellites_gdf = gpd.GeoDataFrame(new_satellite_data, geometry='geometry', crs="EPSG:4326").to_crs(epsg=EQUAL_AREA_EPSG)
    new_satellites_gdf['new_satellite'] = True
    return new_satellites_gdf

# Main Function
def calculate_average_added_coverage_simulated_annealing(
    num_runs=10,
    simulated_annealing_iterations=100,
    initial_temperature=1000,
    cooling_rate=0.95,
    restart_threshold=10,
    buffer_radius=BUFFER_RADIUS,
    equal_distance_epsg=EQUAL_DISTANCE_EPSG, 
    equal_area_epsg=EQUAL_AREA_EPSG, 
    perturbation_distance=5,
    random_seed=None
):
    """
    Calculate the average added coverage across multiple runs of the simulated annealing function.
    """
    # Load map and existing satellite data
    land_map, _ = load_map_data()
    existing_satellites_gdf = load_existing_satellites()
    existing_satellites_gdf["new_satellite"] = False  # Mark existing satellites

    added_coverages = []

    for run in range(num_runs):
        print(f"Run {run + 1}/{num_runs}")

        # Generate new satellites and combine with existing satellites
        new_satellites_gdf = generate_new_satellites(num_satellites=30)
        satellite_gdf = pd.concat([existing_satellites_gdf, new_satellites_gdf], ignore_index=True)

        # Perform simulated annealing optimization
        _, optimized_land_coverage, _, _ = random_restart_simulated_annealing(
            satellite_gdf=satellite_gdf,
            map=land_map,
            new_satellite_column="new_satellite",
            buffer_radius=buffer_radius,
            equal_distance_epsg=equal_distance_epsg,
            equal_area_epsg=equal_area_epsg,
            perturbation_distance=perturbation_distance,
            num_iterations=simulated_annealing_iterations,
            initial_temp=initial_temperature,
            cooling_rate=cooling_rate,
            restart_threshold=restart_threshold,
            random_seed=random_seed,
        )

        # Calculate the added coverage
        initial_coverage = calculate_land_coverage(existing_satellites_gdf, land_map, buffer_radius)
        added_coverage = optimized_land_coverage - initial_coverage
        added_coverages.append(added_coverage)

    # Calculate the average added coverage
    average_added_coverage = np.mean(added_coverages)
    print(f"\n--- Average Added Coverage ---")
    print(f"Average Added Coverage over {num_runs} runs: {average_added_coverage:.2f} km²")
    return average_added_coverage

average_coverage = calculate_average_added_coverage_simulated_annealing()
print(f"Calculated Average Added Coverage with Simulated Annealing: {average_coverage:.2f} km²")
