# trying to use and work through techniques outlined in Python for Geographic Data Analysis
# (https://python-gis-book.readthedocs.io/en/develop/index.html)

import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pyproj as pyj
import pandas as pd
import geopandas as gpd
import tle_processing as tlp

# # constructing the relative data directory path
dirpath = Path(__file__).parent.resolve() / ".." / "data"

# defining a common projection for data
common_crs = pyj.CRS.from_proj4("+proj=eck4")

# create a GeoDataFrame object
starlink = tlp.tles_to_geodataframe(input_file_name=dirpath/"starlink_tle_06OCT2024_21_22.txt",
                                      time=datetime.now(timezone.utc),
                                      crs_string="4326",
                                      add_radius=True,
                                      add_area=True,
                                      set_crs=True)

# projecting the data to the Eckert IV equal area projection
starlink = starlink.to_crs(common_crs)

# create a map projection that can be reasoned about
# and project it to the same map projection
map_filepath = dirpath / "ne_10m_land.zip"
world_map = gpd.read_file(filename=map_filepath)
world_map = world_map.to_crs(common_crs)

# Tested that this map file does actually label the areas that are land
# so a spatial join will probably be able to get us land area
# print(world_map.head())

# showing the plot
fig, ax = plt.subplots()
world_map.plot(ax=ax)
starlink.plot(ax=ax, color="red", markersize = 1)
plt.show()