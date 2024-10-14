# select import statements
from pathlib import Path
from datetime import datetime, timezone
from pyorbital.orbital import Orbital


# full package import statements
import pandas as pd
import numpy as np
import pyproj as pyj
import geopandas as gpd
import shapely as sp


def tle_lines_to_lists(lines: list, filter_dtc=False, time=datetime.now(timezone.utc)):
    """Takes a list of strings containing the lines of TLE data and converts them
    to a set of identifiers, coordinates, and thrown errors

    Args:
        lines (list): The TLE data in the form of a list of strings
        filter_dtc (bool, optional): Whether or not to exclude [DTC] labeled satellites.
            Defaults to False.
        time (datetime, optional): The time to be used for calculation of orbital parameters.
            Defaults to datetime.now(timezone.utc).

    Returns:
        tuple(satellites, longitudes, latitudes, altitudes, not_implemented_errors, crashed_errors):
            a tuple of lists with the following parameters
                - satellites (list(str)): the name of each satellite
                - longitudes (list(float)): the longitude of each satellite in decimal degrees
                - latitudes (list(float)): the latitude of each satellite in decimal degrees
                - altitudes (list(float)): the altitude of each satellite in kilometers
                - not_implemented_errors (list(str)): the name of each satellite that calculations fail for
                - crashed_errors (list(str)): the name of each satellite that is calculated to have crashed
    """
    # Initialize lists for satellite data
    satellites = []  # parameters directly pulled from TLE lines
    longitudes, latitudes, altitudes = (
        [],
        [],
        [],
    )  # parameters calculated with pyorbital
    not_implemented_errors, crashed_errors = [], []  # error code

    for i in range(0, len(lines), 3):

        name = lines[i]
        tle_line_1 = lines[i + 1]
        tle_line_2 = lines[i + 2]

        if filter_dtc and "DTC" in name:
            continue

        # I kept getting NotImplementedError:
        # 'Mode "Near-space, simplified equations" not implemented', so
        # I needed this try/except block to filter out error-throwing
        # satellites

        # I also got some satellites throwing errors that they have been calculated to
        # have crashed, so I put a branch handling that eventuality

        try:
            # create an Orbital object
            orbital_object = Orbital(satellite=name, line1=tle_line_1, line2=tle_line_2)
            # extract the longitude, latitude, and altitude from the object
            current_lon, current_lat, current_alt = orbital_object.get_lonlatalt(time)

            # append the values to the appropriate lists
            satellites.append(name)
            longitudes.append(current_lon)
            latitudes.append(current_lat)
            altitudes.append(current_alt)

        except NotImplementedError:
            not_implemented_errors.append(
                name
            )  # Collect errors for logging if necessary
        except Exception as e:
            if "crash" in str(e):
                crashed_errors.append(name)

    return (
        satellites,
        longitudes,
        latitudes,
        altitudes,
        not_implemented_errors,
        crashed_errors,
    )


# Function to process TLEs into a DataFrame object
def tles_to_dataframe(
    input_file_name: str, time: datetime, calculate_radii_areas=False, angle=45.0
):
    """Take a downloaded text file of TLE data generated at a specific datetime and
    generate a pandas DataFrame from the data as longitude, latitude, altitude triplets

    Args:
        input_file_name (str): _description_
        time (datetime): _description_
        calculate_radii_areas (bool, optional): _description_. Defaults to False.
        angle (float, optional): _description_. Defaults to 45.0.

    Returns:
        _type_: _description_
    """
    # Define data directory and input file
    dirpath = Path(__file__).parent.resolve() / ".." / "data" / "tle_text_files"
    input_file_path = dirpath / input_file_name

    # Extract the TLEs from the file and strip whitespace
    with open(input_file_path, "r") as input_file:
        tle_lines = [line.strip() for line in input_file.readlines()]

    # Ensure each TLE has 3 lines
    if len(tle_lines) % 3 != 0:
        print(
            "The TLE file is unbalanced and does not have the appropriate number of lines"
        )
        return None

    # turn the lines from the text file into orbital parameters
    (
        satellites,
        longitudes,
        latitudes,
        altitudes,
        not_implemented_errors,
        crashed_errors,
    ) = tle_lines_to_lists(tle_lines, filter_dtc=True, time=time)

    # create the output DataFrame Object
    data = {
        "Satellites": satellites,
        "Longitudes": longitudes,
        "Latitudes": latitudes,
        "Altitudes": altitudes,
    }

    # if requested calculate radii and projected areas and add them to the data
    if calculate_radii_areas:
        radii = [altitude * np.tan(angle / 2.0) for altitude in altitudes]
        areas = [np.pi * radius**2 for radius in radii]
        data["Radii"] = radii
        data["Areas"] = areas

    output = pd.DataFrame(data)

    # If there were any satellites with NotImplementedError, print them
    if not_implemented_errors:
        print(
            f"NotImplementedErrors encountered for satellites: {', '.join(not_implemented_errors)}"
        )

    # If there were any satellites that are calculated to have crashed, print them
    if crashed_errors:
        print(
            f"The following satellites were calculated to have crashed: {', '.join(crashed_errors)}"
        )

    return output


# def tles_to_geodataframe(input_file_name: str, time: datetime,
#                          buffer_points = True, angle = 45.0):
#     """Take a downloaded text file of TLE data generated at a specific datetime and
#     generate a pandas DataFrame from the data as longitude, latitude, altitude triplets

#     Args:
#         input_file_name (str): _description_
#         time (datetime): _description_
#         buffer_points (bool, optional): _description_. Defaults to True.
#         angle (float, optional): _description_. Defaults to 45.0.

#     Returns:
#         _type_: _description_
#     """
#     # create the DataFrame object from the TLE file
#     df = tles_to_dataframe(input_file_name=input_file_name,
#                            time=time,
#                            angle=angle)

#     # adding the Shapely object geometries from the data
#     if buffer_points:

#         # approach 3
#         df["geometry"] = gpd.points_from_xy(x=df["Longitude"],
#                                             y=df["Latitude"])


#         # make the GeoDataFrame object and then do the buffering
#         gdf = gpd.GeoDataFrame(df)

#         # approach 4
#         gdf["geometry"] = gdf["geometry"].buffer(distance=gdf["Radius"])

#         # # approach 3
#         # gdf["geometry"] = gdf["geometry"].combine(df["Radius"],
#         #                                           lambda x, y: x.buffer(y))

#         # # create the point objects as a pandas Series object of Shapely Point objects
#         # df["points"] = gpd.points_from_xy(x=df["Longitude"],
#         #                                 y=df["Latitude"])

#         # # approach 1
#         # # define an anonymous function for generating a buffer Shapely Polygon object
#         # # from each point with a given buffer and use that within a pd.Combine call
#         # df["geometry"] = df["points"].combine(df["Radius"],
#         #                                       lambda x, y: x.buffer(y))

#         # # approach 2
#         # geometry_list = []
#         # for i in range(len(df["points"])):
#         #     geometry_list.append(sp.buffer(df["points"].iloc[i], df["Radius"].iloc[i]))

#         # df["geometry"] = geometry_list
#     else:
#         df["geometry"] = gpd.points_from_xy(x=df["Longitude"],
#                                         y=df["Latitude"])

#         gdf = gpd.GeoDataFrame(df)

#     return gdf

# simple impromptu test code
# test = tles_to_dataframe("example_starlink_tle.txt", datetime.now(timezone.utc))
# print(test.head())
# test = tles_to_dataframe("starlink_tle_06OCT2024_21_22.txt", datetime.now(timezone.utc))
# test = tles_to_geodataframe(input_file_name="starlink_tle_06OCT2024_21_22.txt",
#                             time=datetime.now(timezone.utc))
# print(test.head())
# print(len(test))

# EARTH_SURFACE_AREA_SQ_KM = 509600000
# EARTH_LAND_AREA_SQ_KM = 148326000
# print("Maximum Possible Fraction of Earth's Surface Area Covered by Starlink: ",
#       sum(test["Area"])/(EARTH_SURFACE_AREA_SQ_KM))
# print("Maximum Possible Fraction of Earth's Land Area Covered by Starlink: ",
#       sum(test["Area"])/(EARTH_LAND_AREA_SQ_KM))
# print("Median Area Covered by a Satellite: ", np.median(test["Area"]), "square kilometers")
# print("Median Area Covered by a Satellite in terms of Earth Surface Area Coverage: ",
#       (np.median(test["Area"])/EARTH_SURFACE_AREA_SQ_KM) * 100, "percent")
