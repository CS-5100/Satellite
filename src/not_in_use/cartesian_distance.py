import numpy as np
from pyorbital.orbital import Orbital
from datetime import datetime

# Constants for Earth (if converting to Earth-centered coordinates)
a = 6378.137  # semi-major axis (km)
b = 6356.752  # semi-minor axis (km)

# Satellite names to extract TLE from the file
satellite1 = "STARLINK-6317"
satellite2 = "STARLINK-30090"

# File path of the most recent TLE
file_path = "/data/starlink_tle_2024-10-11_17-39-18.txt"

def read_tle_from_file(file_path, satellite):
    """
    Reads a text file containing the specific TLE data and returns it as a dictionary with satellite names as keys.
    
    Parameters:
    - file_path: The path to the file containing the TLE data.
    
    Returns:
    - tle_dict: A dictionary with satellite names as keys and their TLE data as values.
    """
    tle_dict = {}

    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 3):
            # Extract satellite name, TLE line 1 and TLE line 2
            satellite_name = lines[i].strip()
            if satellite_name == satellite:
                tle_line1 = lines[i + 1].strip()
                tle_line2 = lines[i + 2].strip()

                # Store in dictionary
                tle_dict[satellite_name] = [satellite_name, tle_line1, tle_line2]

            else: 
                continue
            
    return tle_dict

def tle_to_lat_long_alt(tle):
    satellite1 = Orbital(tle[0], line1=tle[1], line2=tle[2])

    # Get current time
    now = datetime.utcnow()

    # Get longitude, latitude, and altitude for the satellite
    lon, lat, alt = satellite1.get_lonlatalt(now)

    # Return the extracted longitude, latitude and altitude
    return lon, lat, alt

def geodetic_to_ecef(lat, lon, alt):
    # Convert degrees to radians
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    
    # Compute eccentricity squared
    ecc = 1 - (b ** 2 / a ** 2)
    
    # Radius of curvature in the prime vertical
    N_lat = a / np.sqrt(1 - (ecc * (np.sin(lat_r) ** 2)))
    
    # Calculate ECEF coordinates
    X = (N_lat + alt) * np.cos(lat_r) * np.cos(lon_r)
    Y = (N_lat + alt) * np.cos(lat_r) * np.sin(lon_r)
    Z = (((1 - ecc) * N_lat) + alt) * np.sin(lat_r)
    
    return X, Y, Z

# TLE from file
tle1 = read_tle_from_file(file_path, satellite1).get(satellite1)
tle2 = read_tle_from_file(file_path, satellite2).get(satellite2)


# Example TLE extracted lat, lon, alt (in degrees and km)
lon1, lat1, alt1 = tle_to_lat_long_alt(tle1)
lon2, lat2, alt2 = tle_to_lat_long_alt(tle2)

# Convert to ECEF
X1, Y1, Z1 = geodetic_to_ecef(lat1, lon1, alt1)
X2, Y2, Z2 = geodetic_to_ecef(lat2, lon2, alt2)

# Calculate Euclidean distance between the two satellites in the ECEF coordinate system
distance = np.sqrt(((X2 - X1)**2) + ((Y2 - Y1)**2) + ((Z2 - Z1)**2))

# Print the result
print(f"Distance = {distance:.2f} km")
