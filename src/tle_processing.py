from pyorbital.orbital import Orbital
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

# # Function to calculate lat, long, and altitude for a given TLE
# def process_tle(tle_name: str, tle_line1: str, tle_line2: str):
    
#     # define an orbital object with a given name
#     satellite = Orbital(tle_name, line1=tle_line1, line2=tle_line2)
    
#     # TODO: need to change this to take the time the file was created into account
#     # maybe pass a time into the function
#    current_time = datetime.now(timezone.utc)
    
#     # return the longitude, latitude, and altitude of the satellite
#     return satellite.get_lonlatalt(current_time)

# Function to process TLEs into a DataFrame object
def tles_to_dataframe(input_file_name: str, time: datetime):
    
    # Define data directory and input file
    dirpath = Path(__file__).parent.resolve() / ".." / "data"
    input_file_path = dirpath / input_file_name
    
    # Extract the TLEs from the file and strip whitespace
    with open(input_file_path, 'r') as input_file:
        tle_lines = [line.strip() for line in input_file.readlines()]
    
    # Ensure each TLE has 3 lines
    if len(tle_lines) % 3 != 0:
        print("The TLE file is unbalanced and does not have the appropriate number of lines")
        return None

    # Initialize lists for storing satellite data
    satellites = []
    longitudes = []
    latitudes = []
    altitudes = []
    errors = []

    # across all TLE entries
    for i in range(0, len(tle_lines), 3):
        
        name = tle_lines[i]
        tle_line_1 = tle_lines[i+1]
        tle_line_2 = tle_lines[i+2]
        
        if 'DTC' in name:
            continue
        
        # create an Orbital object
        orbital_object = Orbital(satellite=name,
                                 line1=tle_line_1,
                                 line2=tle_line_2)
        
        # I kept getting NotImplementedError:
        # 'Mode "Near-space, simplified equations" not implemented', so
        # I needed this try/except block to filter out error-throwing
        # satellites
        try:
            # extract the longitude, latitude, and altitude from the object
            current_lon, current_lat, current_alt = orbital_object.get_lonlatalt(time)
        except NotImplementedError: 
            continue
        
        # append the values to the appropriate lists
        satellites.append(name)
        longitudes.append(current_lon)
        latitudes.append(current_lat)
        altitudes.append(current_alt)
        
    output = pd.DataFrame({"Satellite": satellites,
                           "Longitude": longitudes,
                           "Latitude": latitudes,
                           "Altitude": altitudes})
    
    return output

test = tles_to_dataframe("starlink_tle_06OCT2024_21_22.txt", datetime.now(timezone.utc))
print(test.head())
print(len(test))