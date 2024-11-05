import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.validation import explain_validity
from shapely.ops import unary_union
from joblib import Parallel, delayed

def calculate_total_coverage(satellite_data: pd.DataFrame, world: gpd.GeoDataFrame, crs='EPSG:6933', simplify_tolerance=100):
    """
    Function to calculate total unique land and water coverage from satellite data.

    Args:
    - satellite_data (pd.DataFrame): DataFrame containing satellite data with latitude, longitude, and coverage radius.
    - world (gpd.GeoDataFrame): GeoDataFrame representing the world landmass geometries.
    - crs (str): The coordinate reference system to use for area calculations (default is EPSG:6933).
    - simplify_tolerance (int): Tolerance level for simplifying geometries to avoid excessive complexity.

    Returns:
    - land_coverage (float): Total land coverage in square kilometers.
    - water_coverage (float): Total water coverage in square kilometers.
    """
    # Reproject world data to the specified CRS
    world = world.to_crs(crs)

    # Validate and fix world geometries
    world['geometry'] = world['geometry'].apply(lambda geom: geom.buffer(0) if not geom.is_valid else geom)

    # Function to create a circular coverage area around the satellite
    def create_coverage_polygon(lat, lon, radius_km, crs):
        num_points = 100
        coverage_radius_deg = radius_km / 110.574  # Convert km to degrees (latitude)
        angles = np.linspace(0, 2 * np.pi, num_points)
        coverage_circle_lat = lat + coverage_radius_deg * np.cos(angles)
        coverage_circle_lon = lon + coverage_radius_deg * np.sin(angles)
        points = [Point(lon, lat) for lon, lat in zip(coverage_circle_lon, coverage_circle_lat)]
        coverage_polygon = Polygon([[p.x, p.y] for p in points])
        coverage_gdf = gpd.GeoSeries([coverage_polygon], crs='EPSG:4326')
        return coverage_gdf.to_crs(crs)[0]

    # Validate and fix geometries
    def validate_geometry(geometry, simplify_tolerance):
        if not geometry.is_valid:
            geometry = geometry.buffer(0)  # Fix invalid geometries
        return geometry.simplify(simplify_tolerance)

    # Process each satellite and create valid coverage polygons
    def process_satellite(i, sat, crs):
        lat = sat['Latitude']
        lon = sat['Longitude']
        radius_km = sat['Coverage Radius (km)']
        coverage_polygon = create_coverage_polygon(lat, lon, radius_km, crs)
        valid_coverage = validate_geometry(coverage_polygon, simplify_tolerance)
        if valid_coverage.is_empty or valid_coverage is None:
            return None
        return valid_coverage

    # Parallel processing to handle satellite data
    all_valid_coverages = Parallel(n_jobs=-1)(delayed(process_satellite)(i, sat, crs) for i, sat in satellite_data.iterrows())
    all_valid_coverages = [coverage for coverage in all_valid_coverages if coverage is not None]

    # Combine satellite coverages to avoid overlap
    combined_coverage = unary_union(all_valid_coverages)

    # Calculate land and water coverage
    land_coverage = world.intersection(combined_coverage).area.sum() / 1e6  # m² to km²
    total_combined_area = combined_coverage.area / 1e6  # Total coverage in km²
    water_coverage = total_combined_area - land_coverage

    return land_coverage, water_coverage

# Example usage:
# Load the world data and satellite data
shapefile_path = r'C:\Users\danny\PycharmProjects\Satellite\data\ne_110m_admin_0_countries.shp'
world = gpd.read_file(shapefile_path)
file_path = r'C:\Users\danny\PycharmProjects\Satellite\data\satellite_coverage.csv'
satellite_data = pd.read_csv(file_path)
satellite_data['Coverage Radius (km)'] = np.sqrt(satellite_data['Coverage Area (km^2)'] / np.pi)

# Call the function
land_coverage, water_coverage = calculate_total_coverage(satellite_data, world)

# Output results
print(f"Total Unique Land Coverage: {land_coverage:.2f} km²")
print(f"Total Unique Water Coverage: {water_coverage:.2f} km²")
