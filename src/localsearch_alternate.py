from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import tle_processing as tlp
from shapely.ops import unary_union
from shapely.geometry import Point
from shapely.affinity import translate
from random import randint

# Constants
EARTH_SURFACE_AREA_SQ_KM = 509600000
EARTH_LAND_AREA_SQ_KM = 148326000
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933
PLOT_INITIAL = True
PLOT_FINAL = True

def load_map_data():
    # Load land and ocean map data (adjust paths as necessary)
    dirpath = Path(__file__).parent.resolve() / ".." / "data"

    land_filepath = dirpath / "map_data" / "ne_10m_land_scale_rank.zip"
    ocean_filepath = dirpath / "map_data" / "ne_10m_ocean_scale_rank.zip"
    land_map = gpd.read_file(filename=land_filepath).to_crs(epsg=EQUAL_AREA_EPSG)
    ocean_map = gpd.read_file(filename=ocean_filepath).to_crs(epsg=EQUAL_AREA_EPSG)
    return land_map, ocean_map

def load_existing_satellites(show_head=False):
    # Download and convert TLE data to GeoDataFrame for existing satellites
    starlink_current_tle_list = tlp.download_current_tles_as_list()
    starlink_current = tlp.tles_to_dataframe(raw_tle_list=starlink_current_tle_list)
    starlink_gdf = tlp.tle_dataframe_to_geodataframe(starlink_current)
    starlink_gdf = starlink_gdf.to_crs(epsg=EQUAL_AREA_EPSG)
    # we don't necessarily always want the head of the dataframe to print
    if show_head:
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
    subset = gdf[gdf[new_satellite_column]].sample(frac=0.5, random_state=random_state)

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

def local_search_optimization(satellite_gdf, land_map, new_satellite_column, num_iterations=10):
    best_gdf = satellite_gdf.copy()
    best_land_coverage = calculate_land_coverage(best_gdf, land_map)

    for iteration in range(num_iterations):
        print(f"Iteration {iteration + 1}/{num_iterations}")
        
        # Perturb positions of new satellites only
        new_gdf = best_gdf.copy()
        perturb_positions(new_gdf, new_satellite_column, random_state=randint(0, 10000))
        
        # Calculate new land coverage area
        new_land_coverage = calculate_land_coverage(new_gdf, land_map)
        
        if new_land_coverage > best_land_coverage:
            print(f"Improvement found: {new_land_coverage:.2f} km² vs {best_land_coverage:.2f} km²")
            best_gdf, best_land_coverage = new_gdf, new_land_coverage
        else:
            print("No improvement in land coverage.")
    
    print(f"Optimized land coverage area: {best_land_coverage:.2f} km²")
    return best_gdf, best_land_coverage

# Load existing satellite data from TLEs and flag new satellites
existing_satellites_gdf = load_existing_satellites()
existing_satellites_gdf['new_satellite'] = False  # Mark existing satellites

# Generate new satellites and combine with existing ones
new_satellites_gdf = generate_new_satellites(num_satellites=30)
satellite_gdf = pd.concat([existing_satellites_gdf, new_satellites_gdf])

# Load map data
land_map, ocean_map = load_map_data()

# plot current set of satellites in addition to the randomly initialized satellites,
# in point form
if PLOT_INITIAL:
    initial_fig, initial_ax = plt.subplots()
    print(existing_satellites_gdf.crs)
    print(new_satellites_gdf.crs)
    print(land_map.crs)
    print(ocean_map.crs)
    land_map.plot(ax=initial_ax, color="#228B22")
    ocean_map.plot(ax=initial_ax, color="#246BCE")
    existing_satellites_gdf.plot(ax=initial_ax, color="#C51E3A", markersize=1)
    # the below color is apparently known as International Orange (Aerospace) and used in the aerospace industry
    new_satellites_gdf.plot(ax=initial_ax, color="#FF4F00", markersize=1)
    plt.show()

# Run local search with perturbation targeted to new satellites
optimized_gdf, optimized_land_coverage = local_search_optimization(
    satellite_gdf, land_map, 'new_satellite', num_iterations=20
)

# Print out details of the newly generated satellites
print("\nNew Satellites' Final Positions After Optimization:")
print(optimized_gdf[optimized_gdf['new_satellite']][['Satellite Name', 'Latitude', 'Longitude', 'geometry']])
