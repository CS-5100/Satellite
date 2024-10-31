from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import random
import pyproj
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

# Global Parameters
PLOT = True
BUFFER_RADIUS = 12065
PERTURB_DISTANCE_KM = 500 # we may want to make this something that decays exponentially
BUFFER_PLOTS = False

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

def generate_new_satellites(num_satellites=60, input_map = None, true_random = True):
    
    # if the user wants to sample points from a predetermined map,
    # but no map is available to sample from, tell the user to provide one and fail
    if input_map is None and true_random is False:
        print("you must specify a map GeoDataFrame to sample from.")
        return None
    
    # copy the map so that it does not get mutated
    map = input_map.copy()
    
    # to sample latitude and longitude, the map needs to be in a geodetic coordinate system
    # so coerce the map into a worldwide geodetic coordinate reference system
    if map is not None and map.crs != pyproj.CRS.from_epsg(4326):
        map = map.to_crs(epsg=4326)
    
    if true_random:
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
    
    else:
        # sample from a random set of points drawn from each defined map geometry
            
        # must have at least one point per geometry sampled,
        # but can have more if we want a lot of initial satellites
        num_points_per_polygon = int(np.ceil(float(num_satellites) / len(map)))
            
        # sample an appropriate number of points from each geometry
        all_new_points = map.sample_points(num_points_per_polygon)
            
        # sample a number of satellites from all the generated satellite positions
        new_points = all_new_points.sample(num_satellites)
            
        # once we have a set of generated points as a GeoSeries, we need to build the GeoDataFrame
        new_coordinates = new_points.get_coordinates()
        
        # Generate satellite data with coordinates closer to land
        new_satellite_data = {'Satellite Name': [f'NewSat{i+1}' for i in range(num_satellites)],
                                'Latitude': new_coordinates['y'],
                                'Longitude': new_coordinates['x'],
                                'geometry': new_points
                                }
        
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

def calculate_land_coverage(gdf, map, buffer_radius=BUFFER_RADIUS):
    # Apply a buffer around each satellite to simulate coverage
    buffered_gdf = gdf.copy()
    buffered_gdf['geometry'] = buffered_gdf['geometry'].buffer(buffer_radius)  # Buffer of ~121.065 km radius

    # Find intersection area between satellite coverage and land area
    intersections = buffered_gdf.overlay(map, how='intersection')
    total_land_coverage = intersections['geometry'].area.sum() / (1000**2)  # Sum of non-overlapping areas in km²
    print(f"Current land coverage: {total_land_coverage:.2f} km²")  # Track coverage each iteration
    return total_land_coverage

def local_search_optimization(satellite_gdf, land_map, new_satellite_column, buffer_radius = BUFFER_RADIUS, num_iterations=10):
    best_gdf = satellite_gdf.copy()
    best_land_coverage = calculate_land_coverage(best_gdf, land_map, buffer_radius=buffer_radius)

    for iteration in range(num_iterations):
        print(f"Iteration {iteration + 1}/{num_iterations}")
        
        # Perturb positions of new satellites only
        new_gdf = best_gdf.copy()
        perturb_positions(new_gdf, new_satellite_column, random_state=randint(0, 10000))
        
        # Calculate new land coverage area
        new_land_coverage = calculate_land_coverage(new_gdf, land_map, buffer_radius=buffer_radius)
        
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

# Load map data
land_map, ocean_map = load_map_data()

# Generate new satellites, create a copy of them for plotting, and combine with existing ones
new_satellites_gdf = generate_new_satellites(num_satellites=30, input_map=land_map, true_random=False)
initial_satellites_gdf = new_satellites_gdf.copy()
satellite_gdf = pd.concat([existing_satellites_gdf, new_satellites_gdf])

# Run local search with perturbation targeted to new satellites
optimized_gdf, optimized_land_coverage = local_search_optimization(
    satellite_gdf, land_map, 'new_satellite', num_iterations=30
)

# extracting the new satellites from the GeoDataFrame for plotting purposes
final_satellites_gdf = optimized_gdf[optimized_gdf['new_satellite']].copy()

# Print out details of the newly generated satellites
# print("\nNew Satellites' Final Positions After Optimization:")
# print(optimized_gdf[optimized_gdf['new_satellite']][['Satellite Name', 'Latitude', 'Longitude', 'geometry']])

# plot current set of satellites in addition to the randomly initialized satellites,
# in point form
if PLOT:
    
    # create initial figure
    initial_fig, initial_ax = plt.subplots(2, 1)
    before_ax = initial_ax[0]
    after_ax = initial_ax[1]
    
    # add map data to both plots
    land_map.plot(ax=before_ax, color="#228B22")
    ocean_map.plot(ax=before_ax, color="#246BCE")
    land_map.plot(ax=after_ax, color="#228B22")
    ocean_map.plot(ax=after_ax, color="#246BCE")
    
    # if a plot of buffered circles is requested
    if BUFFER_PLOTS:
        
        # copy current GeoDataFrames to preserve their initial state
        existing_satellites_gdf_buffered = existing_satellites_gdf.copy()
        initial_satellites_gdf_buffered = initial_satellites_gdf.copy()
        final_satellites_gdf_buffered = final_satellites_gdf.copy()
        
        # buffer the Point objects they create by the pre-specified buffer
        existing_satellites_gdf_buffered['geometry'] = existing_satellites_gdf_buffered['geometry'].buffer(BUFFER_RADIUS)
        initial_satellites_gdf_buffered['geometry'] = initial_satellites_gdf_buffered['geometry'].buffer(BUFFER_RADIUS)
        final_satellites_gdf_buffered['geometry'] = final_satellites_gdf_buffered['geometry'].buffer(BUFFER_RADIUS)
        
        # add the initial satellites to the initial map
        existing_satellites_gdf_buffered.plot(ax=before_ax, color="#0E0E10") # Jet Black
        initial_satellites_gdf_buffered.plot(ax=before_ax, color="#FF4F00")
        
        # add the final satellites to the final map
        existing_satellites_gdf_buffered.plot(ax=after_ax, color="#0E0E10") # Jet Black
        final_satellites_gdf_buffered.plot(ax=after_ax, color="#FF4F00")
        
    # otherwise only add points and set their markersize to be 1 for easy visibility
    else:
        existing_satellites_gdf.plot(ax=before_ax, color="#0E0E10", markersize=1) # Jet Black
        # the below color is apparently known as International Orange (Aerospace) and used in the aerospace industry
        initial_satellites_gdf.plot(ax=before_ax, color="#FF4F00", markersize=1)
        
        existing_satellites_gdf.plot(ax=after_ax, color="#0E0E10", markersize=1) # Jet Black
        # the below color is apparently known as International Orange (Aerospace) and used in the aerospace industry
        final_satellites_gdf.plot(ax=after_ax, color="#FF4F00", markersize=1)
        
    # show the plot
    plt.show()
