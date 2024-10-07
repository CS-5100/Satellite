from pyorbital.orbital import Orbital
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import numpy as np
import pyproj as pyj
import geopandas as gpd

# Function to process TLEs into a DataFrame object
def tles_to_dataframe(input_file_name: str, time: datetime, add_radius = False, add_area = False, angle = 45.0):
    
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

    if add_area and not add_radius:
        print("Radii must be calculated if coverage areas are requested, change both parameters to True")
        return None

    # Initialize lists for satellite data
    satellites, longitudes, latitudes, altitudes = [], [], [], []
    not_implemented_errors, crashed_errors = [], []
    
    radii = [] if add_radius else None
    areas = [] if add_area else None

    # across all TLE entries
    for i in range(0, len(tle_lines), 3):
        
        name = tle_lines[i]
        tle_line_1 = tle_lines[i+1]
        tle_line_2 = tle_lines[i+2]
        
        if 'DTC' in name:
            continue
        
        # I kept getting NotImplementedError:
        # 'Mode "Near-space, simplified equations" not implemented', so
        # I needed this try/except block to filter out error-throwing
        # satellites
        
        # I also got some satellites throwing errors that they have been calculated to
        # have crashed, so I put a branch handling that eventuality
        try:
            # create an Orbital object
            orbital_object = Orbital(satellite=name,
                                 line1=tle_line_1,
                                 line2=tle_line_2)
            # extract the longitude, latitude, and altitude from the object
            current_lon, current_lat, current_alt = orbital_object.get_lonlatalt(time)
            
            # append the values to the appropriate lists
            satellites.append(name)
            longitudes.append(current_lon)
            latitudes.append(current_lat)
            altitudes.append(current_alt)
            
            # if a radius needs to be added, add a calculated radius with the given
            # field of view angle
            if add_radius:
                radius = current_alt * np.tan(angle / 2.0)
                radii.append(radius)
                
                if add_area:
                    area = np.pi * radius**2
                    areas.append(area)
                
        except NotImplementedError: 
            not_implemented_errors.append(name)  # Collect errors for logging if necessary
        except Exception as e:
            if "crash" in str(e):
                crashed_errors.append(name)
                
    
    # create the output DataFrame Object
    data = {
        "Satellite": satellites,
        "Longitude": longitudes,
        "Latitude": latitudes,
        "Altitude": altitudes
    }
    
    if add_radius:
        data["Radius"] = radii
        if add_area:
            data["Area"] = areas
    
    output = pd.DataFrame(data)
    
    # If there were any satellites with NotImplementedError, print them
    if not_implemented_errors:
        print(f"NotImplementedErrors encountered for satellites: {', '.join(not_implemented_errors)}")
    
    # If there were any satellites that are calculated to have crashed, print them
    if crashed_errors:
        print(f"The following satellites were calculated to have crashed: {', '.join(crashed_errors)}")
    
    return output

def tles_to_geodataframe(input_file_name: str, time: datetime, crs_string: str, add_radius = False, add_area=False, angle = 45.0):
    
    # create the DataFrame object from the TLE file
    df = tles_to_dataframe(input_file_name=input_file_name,
                           time=time,
                           add_radius=add_radius,
                           add_area=add_area,
                           angle=angle)
    
    # adding the Shapely object geometries from the data
    df["geometry"] = gpd.points_from_xy(x=df["Longitude"],
                                        y=df["Latitude"])
    
    # generate the GeoDataFrame object
    gdf = gpd.GeoDataFrame(df)
    
    # set the Coordinate Reference System from the user input
    gdf = gdf.set_crs(pyj.CRS.from_user_input(crs_string))
    
    return gdf

# simple impromptu test code
# test = tles_to_dataframe("example_starlink_tle.txt", datetime.now(timezone.utc))
# test = tles_to_dataframe("starlink_tle_06OCT2024_21_22.txt", datetime.now(timezone.utc))
# test = tles_to_geodataframe("example_starlink_tle.txt", datetime.now(timezone.utc),
#                             "+proj=eck4 +datum=WGS84", add_radius=True, add_area=True)
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