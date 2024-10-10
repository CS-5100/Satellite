import numpy as np
from pyorbital.orbital import Orbital
from datetime import datetime

# Constants for Earth (if converting to Earth-centered coordinates)
a = 6378.137  # semi-major axis (km)
b = 6356.752  # semi-minor axis (km)

# Example TLE
tle1 = ["STARLINK-1282", 
        "1 45409C 20019BB  24277.72826389  .00028003  00000+0  18758-2 0  2773", 
        "2 45409  53.0537 319.7978 0001531  80.1330 356.1358 15.06379228    11"]

tle2 = ["STARLINK-32206", 
        "1 60273C 24131U   24277.42965278 -.00009268  00000+0 -27528-3 0  2776", 
        "2 60273  53.1574 322.0487 0001076  86.9170 331.2866 15.34232211    10"]

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

# Example TLE extracted lat, lon, alt (in degrees and km)
lon1, lat1, alt1 = tle_to_lat_long_alt(tle1)
lon2, lat2, alt2 = tle_to_lat_long_alt(tle2)

# Convert to ECEF
X1, Y1, Z1 = geodetic_to_ecef(lat1, lon1, alt1)
X2, Y2, Z2 = geodetic_to_ecef(lat2, lon2, alt2)

# Calculate Euclidean distance between the two satellites in the ECEF coordinate system
distance = np.sqrt(((X2 - X1)**2) + ((Y2 - Y1)**2) + ((Z2 - Z1)**2))

# Print the result
print(f"Distance between STARLINK-1282 and STARLINK-32206: {distance:.2f} km")