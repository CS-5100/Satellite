import numpy as np

# Constants for Earth (if converting to Earth-centered coordinates)
a = 6378.137  # semi-major axis (km)
b = 6356.752  # semi-minor axis (km)

def geodetic_to_ecef(lat, lon, alt):
    """
    Converts geodetic coordinates (latitude, longitude, altitude) to ECEF coordinates.
    
    Parameters:
    - lat: Latitude in degrees.
    - lon: Longitude in degrees.
    - alt: Altitude in km.
    
    Returns:
    - ECEF coordinates (X, Y, Z) in km.
    """
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

def calculate_distance(position1, position2):
    """
    Calculates the Euclidean distance between two satellites.
    
    Parameters:
    - position1: Tuple representing (latitude, longitude, altitude) for satellite 1.
    - position2: Tuple representing (latitude, longitude, altitude) for satellite 2.
    
    Returns:
    - Distance between the two satellites in km.
    """
    lat1, lon1, alt1 = position1
    lat2, lon2, alt2 = position2
    
    # Convert each position to ECEF coordinates
    X1, Y1, Z1 = geodetic_to_ecef(lat1, lon1, alt1)
    X2, Y2, Z2 = geodetic_to_ecef(lat2, lon2, alt2)
    
    # Calculate Euclidean distance in the ECEF coordinate system
    distance = np.sqrt((X2 - X1)**2 + (Y2 - Y1)**2 + (Z2 - Z1)**2)
    
    return distance

# Example positions
position1 = (550, 53, 0)    # Satellite 1's position
position2 = (540, 50, 20)   # Satellite 2's position

# Calculate the distance between two satellites
distance = calculate_distance(position1, position2)
print(f"Distance = {distance:.2f} km")