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

# function to add the area to each point, to be used in internal functions
def geodataframe_add_buffer():
    return None

# function to get total area covered by satellites
def geodataframe_get_total_area():
    return None

# function to get land area covered by satellites
def geodataframe_get_land_area():
    return None

# function to get population covered by satellites
# (dependent on getting a reliable population density map)
def geodataframe_get_population_coverage():
    return None