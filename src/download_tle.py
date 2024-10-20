import requests
import pandas as pd
from datetime import datetime, timezone


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



# URL for the satellite positions (Starlink TLE data)
url = "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle"

# Call the function to download and save the TLE data with auto-named file
download_tle_data(url)
