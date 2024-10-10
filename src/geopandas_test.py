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
                                      buffer_points=True,
                                      angle=20.0)

# projecting the data to the Eckert IV equal area projection
starlink = starlink.set_crs(epsg="4326")
starlink = starlink.to_crs(common_crs)

# create a map projection that can be reasoned about
# and project it to the same map projection
land_filepath = dirpath / "ne_10m_land_scale_rank.zip"
ocean_filepath = dirpath / "ne_10m_ocean_scale_rank.zip"
land_map = gpd.read_file(filename=land_filepath)
ocean_map = gpd.read_file(filename=ocean_filepath)
land_map = land_map.to_crs(common_crs)
ocean_map = ocean_map.to_crs(common_crs)

# Tested that this map file does actually label the areas that are land
# so a spatial join will probably be able to get us land area
# print(world_map.head())
# print(starlink.head())


# showing the plot
fig, ax = plt.subplots()
land_map.plot(ax=ax, color="#228B22")
ocean_map.plot(ax=ax, color="#246BCE")
starlink.plot(ax=ax, color="#C51E3A", markersize = 1)
plt.show()