import requests
from datetime import datetime
import pandas as pd

def download_tle_data_to_dataframe(url):
    """
    Downloads TLE data from the given URL and converts it into a pandas DataFrame.
    
    Parameters:
    - url: The URL to download the TLE data from.
    
    Returns:
    - df: A pandas DataFrame containing the TLE data.
    """
    try:
        # Send a GET request to the URL
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors

        # Split the text data into lines
        tle_lines = response.text.splitlines()

        # Prepare lists to store the TLE data
        satellite_names = []
        tle_line1 = []
        tle_line2 = []

        # Process the TLE data (every set of 3 lines corresponds to a satellite)
        for i in range(0, len(tle_lines), 3):
            satellite_names.append(tle_lines[i].strip())
            tle_line1.append(tle_lines[i + 1].strip())
            tle_line2.append(tle_lines[i + 2].strip())

        # Create a pandas DataFrame from the lists
        df = pd.DataFrame({
            'Satellite Name': satellite_names,
            'TLE Line 1': tle_line1,
            'TLE Line 2': tle_line2
        })

        print("TLE data successfully downloaded and converted to a DataFrame")
        return df

    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None

# URL for the satellite positions (Starlink TLE data)
url = "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle"

# # Call the function to download and convert the TLE data to a DataFrame
# df = download_tle_data_to_dataframe(url)

# # Display the DataFrame (if successfully created)
# if df is not None:
#     print(df.head())
