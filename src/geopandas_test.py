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
import shapely as sp

# # constructing the relative data directory path
dirpath = Path(__file__).parent.resolve() / ".." / "data"

# defining a common projection for data
common_crs = pyj.CRS.from_proj4("+proj=eck4")

# create a GeoDataFrame object
starlink = tlp.tles_to_geodataframe(input_file_name=dirpath/"starlink_tle_06OCT2024_21_22.txt",
                                      time=datetime.now(timezone.utc),
                                      buffer_points=False,
                                      angle=20.0)

# projecting the data to the Eckert IV equal area projection
starlink = starlink.set_crs(epsg="4326")
starlink = starlink.to_crs(common_crs)

# selecting a subset of satellites for testing
starlink_sample = starlink.head(20).copy(deep=True)

# buffering the points outside of the generation function
# starlink_sample_buffered = starlink_sample.copy(deep=True)
# starlink_sample_buffered["geometry"] = starlink_sample_buffered["geometry"].combine(starlink_sample_buffered["Radius"],
#                                                                                     lambda x, y: x.buffer(y))
# print(starlink_sample_buffered.head())


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
# starlink_sample.plot(ax=ax, color="#C51E3A", markersize = 1)
# starlink_sample_buffered.plot(ax=ax, color="#C51E3A", markersize = 1)
plt.show()