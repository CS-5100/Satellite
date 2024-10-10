from pyorbital.orbital import Orbital
from datetime import datetime
import numpy as np

def calculate_satellite_distance(tle1, tle2):
    """
    Calculate the Euclidean distance between two satellites given their TLE data.
    
    Parameters:
        tle1: List of two-line element (TLE) data for the first satellite.
        tle2: List of two-line element (TLE) data for the second satellite.
    
    Returns:
        Distance in kilometers between the two satellites.
    """
    # Create Orbital objects for both satellites
    satellite1 = Orbital(tle1[0], line1=tle1[1], line2=tle1[2])
    satellite2 = Orbital(tle2[0], line1=tle2[1], line2=tle2[2])

    # Get current time
    now = datetime.utcnow()

    # Get lon, lat, and altitude for both satellites
    lon1, lat1, alt1 = satellite1.get_lonlatalt(now)
    lon2, lat2, alt2 = satellite2.get_lonlatalt(now)

    # Calculate Euclidean distance using the (lon, lat, alt) coordinates
    # This distance is not strictly geographical but just a 3D distance in the lat-lon-alt space
    distance = np.sqrt((lon2 - lon1)**2 + (lat2 - lat1)**2 + (alt2 - alt1)**2)

    return distance

# Example usage
tle1 = ["STARLINK-1282", 
        "1 45409C 20019BB  24277.72826389  .00028003  00000+0  18758-2 0  2773", 
        "2 45409  53.0537 319.7978 0001531  80.1330 356.1358 15.06379228    11"]

tle2 = ["STARLINK-32206", 
        "1 60273C 24131U   24277.42965278 -.00009268  00000+0 -27528-3 0  2776", 
        "2 60273  53.1574 322.0487 0001076  86.9170 331.2866 15.34232211    10"]

# Calculate distance between the two satellites
distance = calculate_satellite_distance(tle1, tle2)
print(f"Distance between STARLINK-1282 and STARLINK-32206: {distance:.2f} km")
