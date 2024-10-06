# trying to use and work through techniques outlined in Python for Geographic Data Analysis
# (https://python-gis-book.readthedocs.io/en/develop/index.html)

import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pyproj as pyj
import pandas as pd
import geopandas as gpd

# constructing the relative data directory path
dirpath = Path(__file__).parent.resolve() / ".." / "data"

# getting the file path of the current starlink positions
position_filepath = dirpath / "starlink_positions.csv"

# defining a common projection for data
common_crs = pyj.CRS.from_proj4("+proj=eck4")

# importing the positions as a pandas DataFrame object
starlink_raw = pd.read_csv(filepath_or_buffer=position_filepath, header=0)

# renaming the columns because there's whitespace in the column names
# due to how the data was converted to a csv (for now, may encapsulate and modify functions for data conversion)
starlink_raw = starlink_raw.rename(columns={" Latitude" : "latitude",
                             " Longitude" : "longitude",
                             " Altitude" : "altitude"})

# adding the radii to the dataframe
starlink_raw["radius"] = np.pi * (starlink_raw["altitude"] * np.tan(45.0 / 2.0))**2

# adding the coordinates as a Shapely geometry object to the DataFrame
starlink_raw["geometry"] = gpd.points_from_xy(x=starlink_raw["longitude"],
                                              y=starlink_raw["latitude"])

# create a GeoDataFrame object
starlink = gpd.GeoDataFrame(starlink_raw)

# projecting the data to the Eckert IV equal area projection
starlink = starlink.set_crs(pyj.CRS.from_epsg("4326"))
starlink = starlink.to_crs(common_crs)

# create a map projection that can be reasoned about
# and project it to the same map projection
map_filepath = dirpath / "ne_10m_land.zip"
world_map = gpd.read_file(filename=map_filepath)
world_map = world_map.to_crs(common_crs)

# showing the plot
fig, ax = plt.subplots()
world_map.plot(ax=ax)
starlink.plot(ax=ax, color="red", markersize = 1)
plt.show()