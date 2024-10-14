import requests
import pandas as pd
from datetime import datetime, timezone
from tle_processing import tle_lines_to_lists


def download_tle_data(url):
    """
    Downloads TLE data from the given URL and saves it to a text file.
    The filename will automatically include the current date and time.

    Parameters:
    - url: The URL to download the TLE data from.
    """
    try:
        # Send a GET request to the URL
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors

        # Get the current date and time to format the filename
        current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"starlink_tle_{current_time}.txt"

        # Save the content to a text file
        with open(filename, "w") as file:
            file.write(response.text)

        print(f"TLE data successfully downloaded and saved to {filename}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")


def download_current_tles_as_list(
    url="https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle",
):
    """Downloads TLE data from the internet and returns them as a list

    Args:
        url (str, optional): The url to download TLE data from.
            Defaults to "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle".

    Returns:
        list | None: Returns the TLE data as a list of strings (each TLE is a set of 3 strings) or None if the GET request failed
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


def download_current_tles_to_dataframe():
    """Downloads TLE data from the internet and returns them as a pandas DataFrame

    Returns:
        pandas.DataFrame | None: Returns the TLE data as pandas DataFrame object or None if the GET request failed.
            The columns of the DataFrame object are "Satellite" , "Longitude", "Latitude", and "Altitude".
            The latitude and longitude are in decimal degrees, the altitude is in kilometers
    """

    # Get the TLEs from Celestrak as a list if possible
    try:
        current_raw_tles = download_current_tles_as_list()
    except requests.exceptions.RequestException:
        print("An error was encountered when trying to download the TLEs")
        return None

    # Ensure each TLE has 3 lines and exit if there's not enough lines
    if len(current_raw_tles) % 3 != 0:
        print("The number of lines in the TLEs downloaded is unbalanced")
        return None

    (
        satellites,
        longitudes,
        latitudes,
        altitudes,
        not_implemented_errors,
        crashed_errors,
    ) = tle_lines_to_lists(current_raw_tles)

    data = {
        "Satellite": satellites,
        "Longitude": longitudes,
        "Latitude": latitudes,
        "Altitude": altitudes,
    }

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

    return pd.DataFrame(data)


# commented this part out because this will generate a text file in any directory this module is used
# URL for the satellite positions (Starlink TLE data)
# url = "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle"

# Call the function to download and save the TLE data with auto-named file
# download_tle_data(url)
