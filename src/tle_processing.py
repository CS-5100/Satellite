# select import statements
from pathlib import Path
from datetime import datetime, timezone
from pyorbital.orbital import Orbital
import requests


# full package import statements
import pandas as pd
import numpy as np
import pyproj as pyj
import geopandas as gpd
import shapely as sp

def download_current_tles_as_list(
    url="https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle",
):
    """Downloads TLE data from the internet and returns them as a list that can be processed by
    tle_lines_to_lists

    Args:
        url (str, optional): The url to download TLE data from.
            Defaults to "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle".

    Returns:
        list[str] | None: Returns the TLE data as a list of strings (each TLE is a set of 3 strings) or None if the GET request failed
    """
    try:
        # Send a GET request to the URL
        response = requests.get(
            url
        )  # may want to put a timeout argument here so that the program does not hang
        response.raise_for_status()  # Raise an exception for HTTP errors

        tle_split_text = response.text.split("\n")  # split the text into lines
        raw_tle_lines = [
            line.strip() for line in tle_split_text
        ]  # strip the white space from each line

        return raw_tle_lines[
            :-1
        ]  # return all the lines but the last one, which is a trailing newline

    # error handling
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")

        return None
    
def convert_TLE_text_file_to_list(input_file_path: Path):
    """Takes a text file of raw TLE lines and returns a list that can be processed by
    tle_lines_to_lists
    
    Args:
        input_file_path (Path): the path to the text file containing the raw TLE lines

    Returns:
        list[str] | None: returns a list of raw TLE lines where each 3 lines defines a TLE coordinate of a satellite
    """
    # test if the file exists / file path is valid
    if not input_file_path.exists():
        print("This file does not seem to exist. Please use a valid pathlib Path object to an existing data file.")
        return None
    
    # Extract the TLEs from the file and strip whitespace
    with open(input_file_path, "r") as input_file:
        tle_lines = [line.strip() for line in input_file.readlines()]
        
    # sometimes the text file ends with a newline and therefore
    # the last string is empty and needs removal
    if tle_lines[-1] == '':
        tle_lines.pop(-1)
        
    return tle_lines
    

def tle_lines_to_lists(lines: list, filter_dtc=False, time=datetime.now(timezone.utc)):
    """Takes a list of strings containing the raw lines of TLE data and converts them
    to a set of identifiers, coordinates, and thrown errors

    Args:
        lines (list[str]): The TLE data in the form of a list of strings
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
    raw_tle_list: list, time = datetime.now(timezone.utc), filter_dtc = True, calculate_radii_areas=False, angle=45.0
):
    """Takes a list of raw TLE strings generated with either convert_TLE_text_file_to_list or download_current_tles_as_list,
    processes them with tle_lines_to_lists to extract meaningful satellite parameters, and then outputs a Pandas DataFrame object
    with the data as longitude, latitude, altitude triplets along with radii and areas if requested.

    Args:
        raw_tle_list (list): _description_
        time (_type_, optional): _description_. Defaults to datetime.now(timezone.utc).
        filter_dtc (bool, optional): _description_. Defaults to True.
        calculate_radii_areas (bool, optional): _description_. Defaults to False.
        angle (float, optional): _description_. Defaults to 45.0.

    Returns:
        _type_: _description_
    """
    # Ensure each TLE has 3 lines
    if len(raw_tle_list) % 3 != 0:
        print(
            "The TLE list is unbalanced and does not have the appropriate number of lines"
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
    ) = tle_lines_to_lists(raw_tle_list, filter_dtc=filter_dtc, time=time)

    # create the output DataFrame Object
    data = {
        "Satellite": satellites,
        "Longitude": longitudes,
        "Latitude": latitudes,
        "Altitude": altitudes,
    }

    # if requested calculate radii and projected areas and add them to the data
    if calculate_radii_areas:
        radii = [altitude * np.tan(angle / 2.0) for altitude in altitudes]
        areas = [np.pi * radius**2 for radius in radii]
        data["Radius"] = radii
        data["Area"] = areas

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


# def download_current_tles_to_dataframe(filter_dtc = True):
#     """Downloads TLE data from the internet and returns them as a pandas DataFrame
    
#     Args:
#         filter_dtc (bool, optional): Whether or not to filter out DTC satellites. defaults to True

#     Returns:
#         pandas.DataFrame | None: Returns the TLE data as pandas DataFrame object or None if the GET request failed.
#             The columns of the DataFrame object are "Satellite" , "Longitude", "Latitude", and "Altitude".
#             The latitude and longitude are in decimal degrees, the altitude is in kilometers
#     """

#     # Get the TLEs from Celestrak as a list if possible
#     try:
#         current_raw_tles = download_current_tles_as_list()
#     except requests.exceptions.RequestException:
#         print("An error was encountered when trying to download the TLEs")
#         return None

#     # Ensure each TLE has 3 lines and exit if there's not enough lines
#     if len(current_raw_tles) % 3 != 0:
#         print("The number of lines in the TLEs downloaded is unbalanced")
#         return None

#     (
#         satellites,
#         longitudes,
#         latitudes,
#         altitudes,
#         not_implemented_errors,
#         crashed_errors,
#     ) = tle_lines_to_lists(current_raw_tles, filter_dtc=filter_dtc)

#     data = {
#         "Satellite": satellites,
#         "Longitude": longitudes,
#         "Latitude": latitudes,
#         "Altitude": altitudes,
#     }

#     # If there were any satellites with NotImplementedError, print them
#     if not_implemented_errors:
#         print(
#             f"NotImplementedErrors encountered for satellites: {', '.join(not_implemented_errors)}"
#         )

#     # If there were any satellites that are calculated to have crashed, print them
#     if crashed_errors:
#         print(
#             f"The following satellites were calculated to have crashed: {', '.join(crashed_errors)}"
#         )

#     return pd.DataFrame(data)