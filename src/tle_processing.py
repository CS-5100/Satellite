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
    
    # create slices for each type of line
    
    # the first element of every three elements is the name
    tle_names = tle_lines[::3]
    
    # the second element of every three elements is the first TLE line
    tle_lines_1 = tle_lines[1::3]
    
    # the last element of every three elements is the second TLE line
    tle_lines_2 = tle_lines[2::3]
    
    # Check that all the TLEs have the appropriate number of lines
    if len(tle_names) != len(tle_lines_1) or len(tle_lines_1) != len(tle_lines_2):
        print("The TLE file is unbalanced and does not have the appropriate number of lines")
        return None
    
    # initialization of empty lists
    satellite = []
    longitude = []
    latitude = []
    altitude = []

    # across all TLE entries
    for i, name in enumerate(tle_names):
        
        if 'DTC' in name:
            continue
        
        # create an Orbital object
        orbital_object = Orbital(satellite=name,
                                 line1=tle_lines_1[i],
                                 line2=tle_lines_2[i])
        
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
        satellite.append(name)
        longitude.append(current_lon)
        latitude.append(current_lat)
        altitude.append(current_alt)
        
    output = pd.DataFrame({"Satellite": satellite,
                           "Longitude": longitude,
                           "Latitude": latitude,
                           "Altitude": altitude})
    
    return output

test = tles_to_dataframe("starlink_tle_06OCT2024_21_22.txt", datetime.now(timezone.utc))
print(test.head())
print(len(test))