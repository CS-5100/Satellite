from pathlib import Path
import geopandas as gpd
import numpy as np
import random
import pyproj
import tle_processing as tlp
from shapely import Point
from shapely.affinity import translate


def load_map_data(map_file_name: str, EPSG: int | str):
    # directory where the map files are stored will be static
    dirpath = (
        Path(__file__).parent.resolve() / ".." / "data" / "map_data" / map_file_name
    )
    output_map = gpd.read_file(filename=dirpath).to_crs(epsg=EPSG)
    return output_map


def load_land_ocean_data(EPSG: int | str):
    land = load_map_data("ne_10m_land_scale_rank.zip", EPSG=EPSG)
    ocean = load_map_data("ne_10m_ocean_scale_rank.zip", EPSG=EPSG)
    return land, ocean


def load_existing_satellites(EPSG: int | str, show_head=False):
    # Download and convert TLE data to GeoDataFrame for existing satellites
    starlink_current_tle_list = tlp.download_current_tles_as_list()
    starlink_current = tlp.tles_to_dataframe(raw_tle_list=starlink_current_tle_list)
    starlink_gdf = tlp.tle_dataframe_to_geodataframe(starlink_current)
    starlink_gdf = starlink_gdf.to_crs(epsg=EPSG)
    # we don't necessarily always want the head of the dataframe to print
    if show_head:
        print(starlink_gdf.head())
    return starlink_gdf


def generate_new_satellites(
    EPSG: int | str,
    num_satellites=60,
    input_map=None,
    true_random=True,
    sample_separate=True,
):
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
            "Satellite Name": [f"NewSat{i+1}" for i in range(num_satellites)],
            "Latitude": [random.uniform(*lat_range) for _ in range(num_satellites)],
            "Longitude": [random.uniform(*lon_range) for _ in range(num_satellites)],
        }
        new_satellite_data["geometry"] = [
            Point(lon, lat)
            for lon, lat in zip(
                new_satellite_data["Longitude"], new_satellite_data["Latitude"]
            )
        ]

        # Construct the GeoDataFrame with geometry explicitly and re-project
        new_satellites_gdf = gpd.GeoDataFrame(
            new_satellite_data, geometry="geometry", crs="EPSG:4326"
        ).to_crs(epsg=EPSG)
        new_satellites_gdf["new_satellite"] = True  # Mark new satellites
        return new_satellites_gdf

    else:
        if sample_separate:
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
            new_satellite_data = {
                "Satellite Name": [f"NewSat{i+1}" for i in range(num_satellites)],
                "Latitude": new_coordinates.y,
                "Longitude": new_coordinates.x,
                "geometry": new_points,
            }

            # Construct the GeoDataFrame with geometry explicitly and re-project
            new_satellites_gdf = gpd.GeoDataFrame(
                new_satellite_data, crs="EPSG:4326"
            ).to_crs(epsg=EPSG)
            new_satellites_gdf["new_satellite"] = True  # Mark new satellites
            return new_satellites_gdf

        else:

            # flatten map
            aggregate_map = map.dissolve()

            # sample from flattened map
            aggregate_map_sample = aggregate_map.sample_points(num_satellites)

            # separate MultiPoint object into distinct points in a GeoSeries
            aggregate_map_sample_points = aggregate_map_sample.explode(
                ignore_index=True
            )

            # get the coordinates from the points
            aggregate_map_coordinates = aggregate_map_sample_points.get_coordinates()

            # Generate satellite data with coordinates closer to land
            new_satellite_data = {
                "Satellite Name": [f"NewSat{i+1}" for i in range(num_satellites)],
                "Latitude": aggregate_map_coordinates.y,
                "Longitude": aggregate_map_coordinates.x,
                "geometry": aggregate_map_sample_points,
            }

            # Construct the GeoDataFrame with geometry explicitly and re-project
            new_satellites_gdf = gpd.GeoDataFrame(
                new_satellite_data, crs="EPSG:4326"
            ).to_crs(epsg=EPSG)
            new_satellites_gdf["new_satellite"] = True  # Mark new satellites
            return new_satellites_gdf


# create function to wrap a point around the globe if necessary
def wrap(point: Point, EPSG: int | str):

    # get the projected bounds of the CRS in cartesian coordinates
    eq_dist_crs = pyproj.CRS.from_epsg(EPSG)
    transformer = pyproj.Transformer.from_crs(
        eq_dist_crs.geodetic_crs, eq_dist_crs, always_xy=True
    )
    min_lon, min_lat, max_lon, max_lat = transformer.transform_bounds(
        *eq_dist_crs.area_of_use.bounds
    )

    # get the total length of the map along each direction
    map_length_x = max_lon - min_lon
    map_length_y = max_lat - min_lat

    # if the point is off the left edge of the earth,
    # shift the coordinate to the right by a map length
    if point.x < min_lon:
        point = translate(point, xoff=map_length_x)

    # if the point is off the right edge of the earth,
    # shift the coordinate to the left by a map length
    if point.x > max_lon:
        point = translate(point, xoff=-map_length_x)

    # if the point is off the bottom edge of the earth,
    # shift the coordinate upwards by a map length
    if point.y < min_lat:
        point = translate(point, yoff=map_length_y)

    # if the point is off the top edge of the earth
    # shift the coordinate downwards by a map length
    if point.y > max_lat:
        point = translate(point, yoff=-map_length_y)

    return point


def perturb_positions(
    gdf,
    new_satellite_column,
    eq_dist_epsg: int | str,
    eq_area_epsg: int | str,
    max_shift_km=500,
    random_state=None,
):
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
    # and get the projected bounds
    subset = subset.to_crs(epsg=eq_dist_epsg)

    # Apply random translations
    gdf.loc[subset.index, "geometry"] = subset["geometry"].apply(
        lambda geom: translate(
            geom,
            xoff=np.random.uniform(-max_shift_m, max_shift_m),
            yoff=np.random.uniform(-max_shift_m, max_shift_m),
        )
    )

    # Wrap the points around the earth to make sure they stay on the globe
    gdf.loc[subset.index, "geometry"] = subset["geometry"].apply(
        lambda geom: wrap(geom, EPSG=eq_dist_epsg)
    )

    # Project back to the equal-area projection
    gdf = gdf.to_crs(epsg=eq_area_epsg)


def calculate_land_coverage(gdf, map, buffer_radius):
    # Apply a buffer around each satellite to simulate coverage
    buffered_gdf = gdf.copy()
    buffered_gdf["geometry"] = buffered_gdf["geometry"].buffer(
        buffer_radius
    )  # Buffer of ~121.065 km radius

    # Find intersection area between satellite coverage and land area
    intersections = buffered_gdf.overlay(map, how="intersection")
    total_land_coverage = intersections["geometry"].area.sum() / (
        1000**2
    )  # Sum of non-overlapping areas in km²
    print(
        f"Current land coverage: {total_land_coverage:.2f} km²"
    )  # Track coverage each iteration
    return total_land_coverage
