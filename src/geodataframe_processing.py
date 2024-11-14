from pathlib import Path
from datetime import datetime, timezone
import geopandas as gpd
import numpy as np
import pyproj
from shapely.errors import GEOSException
from shapely import Point
from shapely.affinity import translate
import tle_processing as tlp


def load_map_data(map_file_name: str, EPSG: int | str):
    """Load a set of map data in the map_data project directory into the environment
    as a geopandas GeoDataFrame

    Args:
        map_file_name (str): the file name of the map data in the map_data directory
        EPSG (int | str): the EPSG code specifying the coordinate reference system to import the map data with

    Returns:
        gpd.GeoDataFrame: a GeoDataFrame with the given map data projected onto a coordinate reference system specified by the input EPSG code
    """
    # directory where the map files are stored will be static
    filepath = (
        Path(__file__).parent.resolve() / ".." / "data" / "map_data" / map_file_name
    )
    output_map = gpd.read_file(filename=filepath).to_crs(epsg=EPSG)
    return output_map


def load_land_ocean_data(EPSG: int | str):
    """Load land and ocean data from the Natural Earth website (https://www.naturalearthdata.com/downloads/10m-physical-vectors/)
    into the environment as a geopandas GeoDataFrame

    Args:
        EPSG (int | str): the EPSG code specifying the coordinate reference system to import the land and ocean data with

    Returns:
        tuple(gpd.GeoDataFrame, gpd.GeoDataFrame):
            a tuple with a GeoDataFrame containing data describing the land masses of the earth as vector polygons as well as
            a GeoDataFrame containing data describing the oceans of the earth as vector polygons
    """
    land = load_map_data("ne_10m_land_scale_rank.zip", EPSG=EPSG)
    ocean = load_map_data("ne_10m_ocean_scale_rank.zip", EPSG=EPSG)
    return land, ocean


def load_existing_satellites(EPSG: int | str, show_head=False):
    """A function creating an automated pipeline to download TLE data from the CelesTrak website
    (https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle) directly into
    a GeoDataFrame object

    Args:
        EPSG (int | str): the EPSG code specifying the coordinate reference system to project the satellite data to
        show_head (bool, optional): whether or not to print the first few entries of the GeoDataFrame to standard output. Defaults to False.

    Returns:
        gpd.GeoDataFrame:
            a geopandas GeoDataFrame with the following parameters
                    - Satellite (str): the name of the satellite
                    - Longitude (float): the longitude of the satellite in decimal degrees
                    - Latitude (float): the latitude of the satellite in decimal degrees
                    - Altitude (float): the altitude of the satellite in kilometers
                    - geometry (shapely Point): the shapely Point object representing the underlying satellite position
    """
    # Download TLEs from the internet into a list
    starlink_current_tle_list = tlp.download_current_tles_as_list()

    # Convert the TLE list to a DataFrame
    starlink_current = tlp.tles_to_dataframe(raw_tle_list=starlink_current_tle_list)

    # Convert the DataFrame into a GeoDataFrame
    starlink_gdf = tlp.tle_dataframe_to_geodataframe(starlink_current)

    # Re-project the GeoDataFrame onto the specified coordinate reference system
    starlink_gdf = starlink_gdf.to_crs(epsg=EPSG)

    # we don't necessarily always want the head of the dataframe to print
    if show_head:
        print(starlink_gdf.head())
    return starlink_gdf


def load_satellites_from_file(
    EPSG: int | str,
    input_filename: str,
    time=datetime.now(tz=timezone.utc),
    show_head=False,
):
    """A function creating an automated pipeline to take a text file with TLE data in the tle_text_files
    project directory, collected at a particular datetime, and load it into the environment directly as a
    GeoDataFrame object

    Args:
        EPSG (int | str): the EPSG code specifying the coordinate reference system to project the satellite data to
        input_filename (str): the name of the file in the tle_text_files directory to import as a GeoDataFrame
        time (datetime, optional): the datetime that the data was collected. Defaults to datetime.now(tz=timezone.utc).
        show_head (bool, optional): whether or not to print the first few entries of the GeoDataFrame to standard output. Defaults to False.

    Returns:
        gpd.GeoDataFrame:
            a geopandas GeoDataFrame with the following parameters
                    - Satellite (str): the name of the satellite
                    - Longitude (float): the longitude of the satellite in decimal degrees
                    - Latitude (float): the latitude of the satellite in decimal degrees
                    - Altitude (float): the altitude of the satellite in kilometers
                    - geometry (shapely Point): the shapely Point object representing the underlying satellite position
    """
    # current schema is designed so that TLE text files are in this particular directory
    filepath = (
        Path(__file__).parent.resolve()
        / ".."
        / "data"
        / "tle_text_files"
        / input_filename
    )
    # Download and convert TLE data to GeoDataFrame for existing satellites
    starlink_tle_list = tlp.convert_TLE_text_file_to_list(input_file_path=filepath)
    starlink_current = tlp.tles_to_dataframe(raw_tle_list=starlink_tle_list, time=time)
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
    random_seed=None,
):
    # TODO: documentation and addressing inconsistencies between import and generation GeoDataFrame columns

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

        # create random number generator
        rng = np.random.default_rng(seed=random_seed)

        # Bounds to center satellites near common land areas (adjust as needed)
        lat_range = (-60, 70)  # Latitude range covering most inhabited regions
        lon_range = (-180, 180)  # Longitude range covering major land areas

        # Generate satellite data with coordinates closer to land
        new_satellite_data = {
            "Satellite": [f"NewSat{i+1}" for i in range(num_satellites)],
            "Longitude": [rng.uniform(*lon_range) for _ in range(num_satellites)],
            "Latitude": [rng.uniform(*lat_range) for _ in range(num_satellites)],
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
                "Satellite": [f"NewSat{i+1}" for i in range(num_satellites)],
                "Longitude": new_coordinates.x,
                "Latitude": new_coordinates.y,
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
                "Satellite": [f"NewSat{i+1}" for i in range(num_satellites)],
                "Longitude": aggregate_map_coordinates.x,
                "Latitude": aggregate_map_coordinates.y,
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
    """takes a shapely Point object and re-maps its x and y coordinates to lie within the bounds
    of a Cartesian coordinate system defined by an EPSG-specified coordinate reference system

    Args:
        point (Point): a shapely Point object describing a point in a Cartesian coordinate system
        EPSG (int | str): the EPSG code of the coordinate reference system to re-map the point within

    Returns:
        point (Point): a shapely Point object describing a valid point in a Cartesian coordinate system within
        the bounds of a given coordinate reference system
    """
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
    
    try:
        x_coordinate = point.x
        y_coordinate = point.y
    except GEOSException:
        print("point had a NaN entry or was empty")
        return point

    # if the point is off the left edge of the earth,
    # shift the coordinate to the right by a map length
    if x_coordinate < min_lon:
        point = translate(point, xoff=map_length_x)

    # if the point is off the right edge of the earth,
    # shift the coordinate to the left by a map length
    if x_coordinate > max_lon:
        point = translate(point, xoff=-map_length_x)

    # if the point is off the bottom edge of the earth,
    # shift the coordinate upwards by a map length
    if y_coordinate < min_lat:
        point = translate(point, yoff=map_length_y)

    # if the point is off the top edge of the earth
    # shift the coordinate downwards by a map length
    if y_coordinate > max_lat:
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


        random_state (int): Random state for reproducibility.
    """
    """Takes an input GeoDataFrame, projects it to a coordinate reference system with a distance preserving projection,
    perturbs the coordinates of the given satellites randomly by an amount up to max_shift_km kilometers, re-maps any points
    that may not specify valid coordinates, and then re-projects the points back to an equal area coordinate reference system

    Args:
        gdf (gpd.GeoDataFrame): The GeoDataFrame containing satellite positions.
        new_satellite_column (str): Column name indicating new satellites.
        eq_dist_epsg (int | str): an equal distance coordinate reference system specified by an EPSG code
        eq_area_epsg (int | str): an equal area coordinate reference system specified by an EPSG code
        max_shift_km (float): Maximum shift in kilometers for each satellite.
        random_state (int, optional): The number specifying a random seed for a random number generator. Defaults to None.
    """
    # create random number generator
    rng = np.random.default_rng(seed=random_state)

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
            xoff=rng.uniform(-max_shift_m, max_shift_m),
            yoff=rng.uniform(-max_shift_m, max_shift_m),
        )
    )

    # Wrap the points around the earth to make sure they stay on the globe
    gdf.loc[subset.index, "geometry"] = subset["geometry"].apply(
        lambda geom: wrap(geom, EPSG=eq_dist_epsg)
    )

    # Project back to the equal-area projection
    gdf = gdf.to_crs(epsg=eq_area_epsg)


def calculate_land_coverage(
    gdf: gpd.GeoDataFrame,
    map: gpd.GeoDataFrame,
    buffer_radius: float,
    print_result=False,
):
    """Overlays the coverage areas of the satellites given in a GeoDataFrame onto a map and
    calculates the area of the intersection of all the objects

    Args:
        gdf (gpd.GeoDataFrame): a GeoDataFrame defining a set of satellites
        map (gpd.GeoDataFrame): a GeoDataFrame defining a map to calculate intersection area for
        buffer_radius (float): The radius of the circle covered by each satellite in meters?
        print_result (bool, optional): whether or not to print the results of the calculation to standard output. Defaults to False.

    Returns:
        float: The total map area covered by the given satellites in square kilometers
    """
    # copy the initial array to prevent mutation
    buffered_gdf = gdf.copy()

    # filter out any points that have NaN coordinates, as they should not contribute
    buffered_gdf = buffered_gdf[buffered_gdf["geometry"].is_valid]

    # Apply a buffer around each satellite to simulate coverage
    buffered_gdf["geometry"] = buffered_gdf["geometry"].buffer(buffer_radius)

    # Find intersection area between satellite coverage and land area
    intersections = buffered_gdf.overlay(map, how="intersection")
    total_land_coverage = intersections["geometry"].area.sum() / (
        1000**2
    )  # Sum of non-overlapping areas in km²
    if print_result:
        print(
            f"Current land coverage: {total_land_coverage:.2f} km²"
        )  # Track coverage each iteration
    return total_land_coverage
