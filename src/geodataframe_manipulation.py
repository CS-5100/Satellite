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

# function to get population covered by satellites
# (dependent on getting a reliable population density map)
def geodataframe_get_population_coverage():
    return None