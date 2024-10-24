import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
from collections import defaultdict
from math import sqrt, acos, pi

# ----------------------------------------------------------
# Encapsulated Function: calculate_unique_coverage_area
# ----------------------------------------------------------

def calculate_unique_coverage_area(satellite_data):
    """
    Calculate the total unique coverage area for a set of satellites.

    Args:
        satellite_data (pd.DataFrame): DataFrame containing satellite information
                                       (Latitude, Longitude, Coverage Area, etc.)

    Returns:
        float: Total unique coverage area in km²
        float: Percentage of Earth's surface covered by the satellites
    """

    # ----------------------------------------------------------
    # Step 1: Process the satellite data
    # ----------------------------------------------------------

    # Calculate the radius of coverage for each satellite based on its coverage area
    # Formula: Radius = sqrt(Coverage Area / pi)
    satellite_data['Coverage Radius (km)'] = np.sqrt(satellite_data['Coverage Area (km^2)'] / np.pi)

    # ----------------------------------------------------------
    # Step 2: Build a BallTree and find neighboring satellites
    # ----------------------------------------------------------

    # Function to calculate the Haversine distance between two latitude/longitude points
    def haversine(lat1, lon1, lat2, lon2):
        """Calculate the great-circle distance between two points on the Earth's surface."""
        from math import radians, cos, sin, sqrt, atan2
        r = 6371  # Earth's radius in kilometers
        dlat = radians(lat2 - lat1)
        dlon = radians(lon2 - lon1)
        a = sin(dlat / 2) ** 2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2) ** 2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        return r * c

    # Build a BallTree using the satellite coordinates (latitude and longitude in radians)
    coords = np.radians(satellite_data[['Latitude', 'Longitude']].values)
    tree = BallTree(coords, metric='haversine')

    # Find neighbors for a satellite based on coverage radius
    def find_neighbors(tree, index, radius):
        """Find neighbors within the given radius using BallTree."""
        query_point = coords[index].reshape(1, -1)
        indices = tree.query_radius(query_point, r=radius / 6371.0)  # Convert radius to radians
        return indices[0]  # Return the indices of the neighbors

    # ----------------------------------------------------------
    # Step 3: Build a graph of overlapping satellites
    # ----------------------------------------------------------

    satellite_graph = defaultdict(list)

    # Populate the graph by connecting satellites whose coverage areas overlap
    for i, sat1 in satellite_data.iterrows():
        r1 = sat1['Coverage Radius (km)']
        neighbors = find_neighbors(tree, i, r1)

        for neighbor in neighbors:
            if neighbor != i:  # Avoid connecting the satellite to itself
                sat2 = satellite_data.iloc[neighbor]
                r2 = sat2['Coverage Radius (km)']
                # Calculate the distance between the two satellites using the Haversine formula
                distance = haversine(sat1['Latitude'], sat1['Longitude'], sat2['Latitude'], sat2['Longitude'])
                if distance < (r1 + r2):  # Connect satellites only if they overlap
                    satellite_graph[sat1['Satellite']].append(sat2['Satellite'])
                    satellite_graph[sat2['Satellite']].append(sat1['Satellite'])

    # ----------------------------------------------------------
    # Step 4: Find connected components of overlapping satellites
    # ----------------------------------------------------------

    def find_connected_components(graph):
        """Find connected components in the graph of overlapping satellites."""
        visited = set()
        components = []

        def dfs(node, component):
            """Depth-First Search to explore and find all satellites in the same connected component."""
            visited.add(node)
            component.append(node)
            for neighbor in graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, component)

        for satellite in graph:
            if satellite not in visited:
                component = []
                dfs(satellite, component)
                components.append(component)

        return components

    connected_components = find_connected_components(satellite_graph)

    # ----------------------------------------------------------
    # Step 5: Calculate the union of coverage areas for each group
    # ----------------------------------------------------------

    def circle_overlap_area(r1, r2, d):
        """Calculate the area of overlap between two circles with radii r1 and r2 and distance d."""
        if d >= r1 + r2:
            return 0  # No overlap
        elif d <= abs(r1 - r2):
            return pi * min(r1, r2) ** 2  # One circle is completely inside the other
        else:
            # Area of overlap between two circles using geometry
            part1 = r1 ** 2 * acos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
            part2 = r2 ** 2 * acos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
            part3 = 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
            return part1 + part2 - part3

    def calculate_circle_area(radius):
        """Calculate the area of a circle given its radius."""
        return pi * radius ** 2

    def calculate_group_union_area(group, satellite_data):
        """Calculate the total unique area covered by a group of overlapping satellites."""
        total_area = 0
        overlap_area = 0
        counted_pairs = set()  # To track satellite pairs already considered for overlap

        # Step 1: Calculate the total area for all satellites in the group
        for satellite in group:
            sat_data = satellite_data[satellite_data['Satellite'] == satellite].iloc[0]
            r = sat_data['Coverage Radius (km)']
            area = calculate_circle_area(r)
            total_area += area

        # Step 2: Calculate the overlapping area within the group
        for i, sat1 in enumerate(group):
            for j, sat2 in enumerate(group):
                if i < j:
                    sat1_data = satellite_data[satellite_data['Satellite'] == sat1].iloc[0]
                    sat2_data = satellite_data[satellite_data['Satellite'] == sat2].iloc[0]

                    # Only calculate overlap if this pair has not been processed
                    if (sat1, sat2) not in counted_pairs:
                        distance = haversine(sat1_data['Latitude'], sat1_data['Longitude'],
                                             sat2_data['Latitude'], sat2_data['Longitude'])
                        r1 = sat1_data['Coverage Radius (km)']
                        r2 = sat2_data['Coverage Radius (km)']

                        if distance < (r1 + r2):
                            overlap = circle_overlap_area(r1, r2, distance)
                            overlap = min(overlap, 0.25 * min(calculate_circle_area(r1),
                                                              calculate_circle_area(r2)))
                            overlap_area += overlap

                        # Mark this pair as processed
                        counted_pairs.add((sat1, sat2))
                        counted_pairs.add((sat2, sat1))

        return total_area - overlap_area

    # ----------------------------------------------------------
    # Step 6: Calculate the total unique coverage area
    # ----------------------------------------------------------

    earth_surface_area_km2 = 510072000  # Earth's total surface area in km²

    # Calculate unique area for each connected component (group of satellites)
    total_unique_area = 0
    for group in connected_components:
        group_unique_area = calculate_group_union_area(group, satellite_data)
        total_unique_area += group_unique_area

    # Calculate percentage of Earth's surface covered
    coverage_percentage = (total_unique_area / earth_surface_area_km2) * 100

    return total_unique_area, coverage_percentage


# ----------------------------------------------------------
# Example Usage
# ----------------------------------------------------------

file_path = r'C:\Users\danny\PycharmProjects\Satellite\data\satellite_coverage.csv'
satellite_data = pd.read_csv(file_path)

total_unique_area, coverage_percentage = calculate_unique_coverage_area(satellite_data)

print(f"Total Unique Coverage Area: {total_unique_area} km²")
print(f"Percentage of Earth's Surface Covered: {coverage_percentage:.2f}%")
