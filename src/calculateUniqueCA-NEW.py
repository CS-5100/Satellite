import pandas as pd
import numpy as np
import networkx as nx
from math import sqrt, pi
from sklearn.neighbors import BallTree

def calculate_unique_coverage_area(satellite_data):
    """
    Calculate the total unique coverage area for a set of satellites using BallTree for efficiency.

    Args:
        satellite_data (pd.DataFrame): DataFrame containing satellite information
                                       (Latitude, Longitude, Coverage Area, etc.)

    Returns:
        float: Total unique coverage area in km²
    """
    # Ensure Latitude and Longitude are numeric and drop rows with NaN values
    satellite_data['Latitude'] = pd.to_numeric(satellite_data['Latitude'], errors='coerce')
    satellite_data['Longitude'] = pd.to_numeric(satellite_data['Longitude'], errors='coerce')
    satellite_data = satellite_data.dropna(subset=['Latitude', 'Longitude'])

    # Calculate the radius of coverage for each satellite based on its coverage area
    satellite_data['Coverage Radius (km)'] = np.sqrt(satellite_data['Coverage Area (km^2)'] / np.pi)
    
    print(f"Processing {len(satellite_data)} satellites...")

    # Create a graph for overlapping satellites
    G = nx.Graph()

    # Convert latitude and longitude to radians for BallTree
    coords = np.radians(satellite_data[['Latitude', 'Longitude']].values)
    
    # Build a BallTree for fast neighbor lookup
    tree = BallTree(coords, metric='haversine')

    # Add nodes (satellites)
    for i, sat in satellite_data.iterrows():
        G.add_node(sat['Satellite Name'], radius=sat['Coverage Radius (km)'], 
                   lat=sat['Latitude'], lon=sat['Longitude'])

    print("Nodes added. Finding overlapping satellites using BallTree...")

    # Function to calculate the Haversine distance between two latitude/longitude points
    def haversine(lat1, lon1, lat2, lon2):
        from math import radians, cos, sin, sqrt, atan2
        r = 6371  # Earth's radius in kilometers
        dlat = radians(lat2 - lat1)
        dlon = radians(lon2 - lon1)
        a = sin(dlat / 2) ** 2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2) ** 2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        return r * c

    # Find overlapping satellites using BallTree
    for i, sat in satellite_data.iterrows():
        r1 = sat['Coverage Radius (km)'] / 6371.0  # Convert radius to radians (Earth's radius in km)
        
        # Query the BallTree for neighbors within the coverage radius
        neighbors = tree.query_radius([coords[i]], r=r1)[0]
        
        for neighbor_index in neighbors:
            if neighbor_index != i:  # Ignore self
                sat2 = satellite_data.iloc[neighbor_index]
                
                # Calculate Haversine distance between the two satellites
                distance = haversine(sat['Latitude'], sat['Longitude'], sat2['Latitude'], sat2['Longitude'])

                # Only add an edge if the satellites actually overlap
                if distance < (sat['Coverage Radius (km)'] + sat2['Coverage Radius (km)']):
                    G.add_edge(sat['Satellite Name'], sat2['Satellite Name'])

    print("Edges added. Finding connected components...")

    # Use NetworkX to find connected components
    connected_components = list(nx.connected_components(G))

    print(f"Found {len(connected_components)} connected components. Calculating unique coverage area...")

    # Function to calculate the area of a circle
    def calculate_circle_area(radius):
        return pi * radius ** 2

    # Calculate the total unique coverage area
    total_unique_area = 0
    for idx, component in enumerate(connected_components):
        component_satellites = satellite_data[satellite_data['Satellite Name'].isin(component)]
        total_area = sum(calculate_circle_area(sat['Coverage Radius (km)']) for _, sat in component_satellites.iterrows())
        
        # Progress tracking for each component
        print(f"Processed component {idx+1}/{len(connected_components)}")

        total_unique_area += total_area

    print(f"Total unique coverage area: {total_unique_area} km²")
    
    return total_unique_area
