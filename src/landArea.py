import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.validation import explain_validity
from shapely.ops import unary_union
from joblib import Parallel, delayed

# Load the shapefile into GeoPandas
shapefile_path = r'C:\Users\danny\PycharmProjects\Satellite\data\ne_110m_admin_0_countries.shp'
world = gpd.read_file(shapefile_path)
world = world.to_crs(epsg=6933)  # Reproject to equal-area projection (EPSG:6933)

# Validate and fix world geometries
def validate_world_geometry(geometry):
    """Check if world geometry is valid; if not, fix it using buffer(0)."""
    if not geometry.is_valid:
        return geometry.buffer(0)  # Fix invalid geometries
    return geometry

world['geometry'] = world['geometry'].apply(validate_world_geometry)

# Load the satellite data
file_path = r'C:\Users\danny\PycharmProjects\Satellite\data\satellite_coverage.csv'
satellite_data = pd.read_csv(file_path)

# Calculate the radius of coverage for each satellite
satellite_data['Coverage Radius (km)'] = np.sqrt(satellite_data['Coverage Area (km^2)'] / np.pi)

# Function to create a circular coverage area around the satellite
def create_coverage_polygon(lat, lon, radius_km, crs):
    """Creates a circular polygon representing the satellite's coverage area."""
    num_points = 100
    coverage_radius_deg = radius_km / 110.574  # Convert km to degrees (latitude)
    angles = np.linspace(0, 2 * np.pi, num_points)

    # Compute the points of the circle (latitude/longitude)
    coverage_circle_lat = lat + coverage_radius_deg * np.cos(angles)
    coverage_circle_lon = lon + coverage_radius_deg * np.sin(angles)

    # Create a polygon using the lat/lon points
    points = [Point(lon, lat) for lon, lat in zip(coverage_circle_lon, coverage_circle_lat)]
    coverage_polygon = Polygon([[p.x, p.y] for p in points])

    # Convert coverage_polygon to GeoSeries and reproject to the specified CRS
    coverage_gdf = gpd.GeoSeries([coverage_polygon], crs='EPSG:4326')
    coverage_gdf = coverage_gdf.to_crs(crs)  # Reproject for consistent area calculations

    return coverage_gdf[0]

# Validate and fix geometries
def validate_geometry(geometry, simplify_tolerance=100):
    """Check if the geometry is valid; if not, attempt to fix it using buffer(0)."""
    if not geometry.is_valid:
        print(f"Invalid geometry: {explain_validity(geometry)}")
        geometry = geometry.buffer(0)  # Fix invalid geometries
        if geometry.is_valid:
            print(f"Geometry fixed using buffer(0)")
        else:
            print("Still invalid after buffer(0)")

    # Simplify the geometry to avoid overly complex areas (reduce the tolerance)
    geometry = geometry.simplify(simplify_tolerance)
    return geometry

# Parallel function to process each satellite
def process_satellite(i, sat, crs):
    lat = sat['Latitude']
    lon = sat['Longitude']
    radius_km = sat['Coverage Radius (km)']

    # Create coverage polygon
    coverage_polygon = create_coverage_polygon(lat, lon, radius_km, crs)

    # Validate the geometry
    valid_coverage = validate_geometry(coverage_polygon)
    if valid_coverage.is_empty or valid_coverage is None:
        print(f"Skipping invalid geometry for satellite index {i}")
        return None

    return valid_coverage

# Use parallel processing to handle multiple satellites at once
crs = world.crs  # Use the CRS of the world data (EPSG:6933)
all_valid_coverages = Parallel(n_jobs=-1)(delayed(process_satellite)(i, sat, crs) for i, sat in satellite_data.iterrows())

# Filter out None values
all_valid_coverages = [coverage for coverage in all_valid_coverages if coverage is not None]

# Step 1: Combine all satellite coverages to avoid overlap
combined_coverage = unary_union(all_valid_coverages)  # Combine all valid coverages into a single geometry

# Step 2: Calculate the land and water coverage by intersecting the combined coverage with land
land_coverage = world.intersection(combined_coverage).area.sum() / 1e6  # Convert from m² to km²

# Step 3: Calculate total area of the combined satellite coverage
total_combined_area = combined_coverage.area / 1e6  # Total coverage area in km²

# Step 4: Calculate water coverage as the remaining coverage
water_coverage = total_combined_area - land_coverage

# Output the results
print(f"Total Unique Land Coverage: {land_coverage:.2f} km²")
print(f"Total Unique Water Coverage: {water_coverage:.2f} km²")
