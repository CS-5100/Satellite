import numpy as np
import pandas as pd
import geopandas as gpd
from pyorbital import tlefile
from pyorbital.orbital import Orbital
from datetime import datetime
import download_tle_to_df as dtd
import calculateUniqueCA as cu  # Assuming this file contains helper functions for unique coverage area

# Constants
EARTH_RADIUS_KM = 6371  # Earth's radius in kilometers

# Step 1: Download existing satellite TLE data
def download_tles():
    # URL for the satellite positions (Starlink TLE data)
    url = "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle"
    return dtd.download_tle_data_to_dataframe(url)

# Step 2: Convert TLE data to Latitude, Longitude, and Altitude using pyorbital
def tle_to_lat_lon_alt(tle_df):
    """Convert TLE DataFrame into latitude, longitude, and altitude."""
    tle_df['Latitude'] = None
    tle_df['Longitude'] = None
    tle_df['Altitude'] = None

    for i, row in tle_df.iterrows():
        tle_line1 = row['TLE Line 1'].strip()
        tle_line2 = row['TLE Line 2'].strip()

        if len(tle_line1) == 69 and len(tle_line2) == 69:
            satellite_name = row['Satellite Name']
            orb = Orbital(satellite_name, line1=tle_line1, line2=tle_line2)
            
            lat, lon, alt = orb.get_lonlatalt(datetime.utcnow())
            
            tle_df.at[i, 'Latitude'] = lat
            tle_df.at[i, 'Longitude'] = lon
            tle_df.at[i, 'Altitude'] = alt
        else:
            print(f"Skipping invalid TLE for satellite {row['Satellite Name']}")

    return tle_df

# Step 3: Convert the DataFrame to a GeoDataFrame
def tle_dataframe_to_geodataframe(tle_df):
    tle_df = tle_to_lat_lon_alt(tle_df)
    tle_df['geometry'] = gpd.points_from_xy(tle_df['Longitude'], tle_df['Latitude'])
    return gpd.GeoDataFrame(tle_df, crs="EPSG:4326")

# Step 4: Generate new satellite positions
def generate_new_satellites(num_satellites=60):
    """Generate random positions for new satellites."""
    new_satellites = {
        'Satellite Name': [f'NEW-SAT-{i}' for i in range(num_satellites)],
        'Latitude': np.random.uniform(-90, 90, size=num_satellites),
        'Longitude': np.random.uniform(-180, 180, size=num_satellites),
        'Altitude': np.random.uniform(500, 2000, size=num_satellites),  # Mock altitudes for LEO satellites
    }
    new_satellites_df = pd.DataFrame(new_satellites)
    new_satellites_df['geometry'] = gpd.points_from_xy(new_satellites_df['Longitude'], new_satellites_df['Latitude'])
    return gpd.GeoDataFrame(new_satellites_df, crs="EPSG:4326")

# Function to calculate the coverage area based on altitude
def calculate_coverage_area(altitude_km):
    """Calculate the coverage area of a satellite based on its altitude."""
    return 2 * np.pi * EARTH_RADIUS_KM**2 * (1 - EARTH_RADIUS_KM / (EARTH_RADIUS_KM + altitude_km))

# Step 5: Add the coverage area to the GeoDataFrame
def add_coverage_area_to_gdf(geo_df):
    """Add a coverage area column to the GeoDataFrame based on satellite altitude."""
    geo_df['Coverage Area (km^2)'] = geo_df['Altitude'].apply(calculate_coverage_area)
    return geo_df

# Step 6: Objective function based on unique coverage area
def objective_function(new_satellites_gdf, existing_satellites_gdf):
    """Calculate the total unique area coverage of the satellite network."""
    # Combine the existing and new satellites
    combined_gdf = pd.concat([new_satellites_gdf, existing_satellites_gdf])
    
    # Add coverage area to the combined GeoDataFrame
    combined_gdf = add_coverage_area_to_gdf(combined_gdf)
    
    # Call the helper from calculateUniqueCA to calculate the coverage
    total_unique_coverage = cu.calculate_unique_coverage_area(combined_gdf)
    
    return total_unique_coverage

# Step 7: Local search algorithm to optimize the placement of new satellites
def local_search(existing_satellites_gdf, iterations=100, step_size=0.01):
    # Generate initial random placement of 60 new satellites
    best_new_satellites = generate_new_satellites()
    best_coverage = objective_function(best_new_satellites, existing_satellites_gdf)

    for i in range(iterations):
        # Perturb the satellite positions slightly to create new solutions
        new_satellites = generate_new_satellites()
        new_coverage = objective_function(new_satellites, existing_satellites_gdf)

        if new_coverage > best_coverage:
            best_new_satellites = new_satellites
            best_coverage = new_coverage
            print(f"Iteration {i}: New best coverage = {best_coverage}")

    return best_new_satellites

# Step 8: Prepare the mock TLE data and run the local search
def prepare_data_with_pyorbital():
    tle_df = download_tles()  # Download existing TLE data
    existing_satellites_gdf = tle_dataframe_to_geodataframe(tle_df)  # Convert to GeoDataFrame
    return existing_satellites_gdf

# Step 9: Run the local search to optimize new satellite placement
existing_satellites_gdf = prepare_data_with_pyorbital()  # Load existing satellite data
print (existing_satellites_gdf)

# Perform local search optimization to place 60 new satellites
optimized_new_satellites = local_search(existing_satellites_gdf, iterations=10)

# Step 10: Output the optimized satellite positions
print(optimized_new_satellites[['Satellite Name', 'Latitude', 'Longitude', 'Altitude']])

