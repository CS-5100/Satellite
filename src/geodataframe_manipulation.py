# select import statements
from pathlib import Path
from datetime import datetime, timezone
from pyorbital.orbital import Orbital


# full package import statements
import pandas as pd
import numpy as np
import pyproj as pyj
import geopandas as gpd
import shapely as sp

# custom package import statements
import tle_processing as tlp


def satellite_and_map_data_intersection_area(satellite_data: gpd.GeoDataFrame, map_data: gpd.GeoDataFrame):
    
    # flatten the geometries contained in each GeoDataFrame into one multi-shape
    satellite_union_shape = satellite_data['geometry'].union_all()
    map_union_shape = map_data['geometry'].union_all()
    
    # get the intersection
    satellite_map_intersection_shape = map_union_shape.intersection(satellite_union_shape)
    
    return satellite_map_intersection_shape.area

def geodataframe_total_unique_area(satellite_data: gpd.GeoDataFrame, deep_copy = True, angle = 45.0):
    
    # seeing if operating on a temporary deep copy fixes the bug
    satellites = satellite_data.copy(deep=deep_copy)
    
    # defining a couple of epsg codes
    equal_distance_epsg = 4087
    equal_area_epsg = 6933
    
    # re-project the data into a Cartesian coordinate system that preserves distances
    satellite_distance = satellites.to_crs(epsg=equal_distance_epsg)
    
    # generate radii from the input GeoDataFrame
    # selected Cartesian coordinate system is in meters,
    # so altitude needs to be converted from kilometers to meters here
    satellite_distance['Radius'] = (satellite_distance['Altitude'] * 1000) * np.tan((np.deg2rad(angle)) / 2.0)
    
    # create internal function here for buffering the first argument with the second argument
    def buffer_point(point: sp.Point, radius: float):
        return point.buffer(radius)
    
    # expand geometry here from points to circles
    satellite_distance['geometry'] = satellite_distance['geometry'].combine(satellite_distance['Radius'],
                                                                            buffer_point)
    
    # project the GeoDataFrame into an equal-area coordinate system
    satellite_area =  satellite_distance.to_crs(epsg=equal_area_epsg)
    
    # filtering out invalid geometries
    satellite_area_valid = satellite_area['geometry'].is_valid
    satellite_area = satellite_area[satellite_area_valid]
    
    # flatten all the areas into a MultiPolygon
    satellite_area_union = satellite_area.union_all()
    
    # get the union area in square kilometers
    satellite_unique_area = (satellite_area_union.area / (1000**2))
    
    return satellite_unique_area


def geodataframe_total_area(satellite_data: gpd.GeoDataFrame, angle = 45.0):
    
    # defining a couple of epsg codes
    equal_distance_epsg = 4087
    equal_area_epsg = 6933
    
    # re-project the data into a Cartesian coordinate system that preserves distances
    satellite_distance = satellite_data.to_crs(epsg=equal_distance_epsg)
    
    # generate radii from the input GeoDataFrame
    # selected Cartesian coordinate system is in meters,
    # so altitude needs to be converted from kilometers to meters here
    satellite_distance['Radius'] = (satellite_distance['Altitude'] * 1000) * np.tan((np.deg2rad(angle)) / 2.0)
    
    # create internal function here for buffering the first argument with the second argument
    def buffer_point(point: sp.Point, radius: float):
        return point.buffer(radius)
    
    # expand geometry here from points to circles
    satellite_distance['geometry'] = satellite_distance['geometry'].combine(satellite_distance['Radius'],
                                                                            buffer_point)
    
    # filtering out invalid geometries
    
    # project the GeoDataFrame into an equal-area coordinate system
    satellite_area =  satellite_distance.to_crs(epsg=equal_area_epsg)
    satellite_area_valid = satellite_area['geometry'].is_valid
    satellite_area = satellite_area[satellite_area_valid]
    
    # create internal function here for getting the area from a point
    def get_area(circle: sp.Polygon):
        return circle.area
    
    # store the area of each circle
    satellite_area['Area'] = satellite_area['geometry'].apply(get_area)
    
    # get the sum of all the areas
    total_possible_area = satellite_area['Area'].sum()
    total_possible_area_sq_km = total_possible_area / (1000**2)
    
    return total_possible_area_sq_km