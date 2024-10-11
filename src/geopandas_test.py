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

# # constructing the relative data directory path and the example tle file path
dirpath = Path(__file__).parent.resolve() / ".." / "data"
example_tle_txt_file = dirpath / "tle_text_files" / "starlink_tle_06OCT2024_21_22.txt"

# defining a common projection for data
common_crs = pyj.CRS.from_proj4("+proj=eck4")

# getting the current time
current_time = datetime.now(timezone.utc)

# create a GeoDataFrame object
starlink = tlp.tles_to_geodataframe(
    input_file_name=example_tle_txt_file,
    time=datetime.now(timezone.utc),
    buffer_points=False,
    angle=20.0,
)

# projecting the data to the Eckert IV equal area projection
starlink = starlink.set_crs(epsg="4326")
starlink = starlink.to_crs(common_crs)

# selecting a subset of satellites for testing
# starlink_sample = starlink.head(20).copy(deep=True)

# creating just a dataframe of starlink data with one points
starlink_minimal = tlp.tles_to_dataframe(
    input_file_name=example_tle_txt_file, time=current_time
).head(1)

starlink_minimal['geometry'] = gpd.points_from_xy(x=starlink_minimal["Longitude"],
                                                  y=starlink_minimal["Latitude"])


starlink_minimal["geometry"] = starlink_minimal["geometry"].combine(starlink_minimal["Radius"], lambda x,y: x.buffer(y))

starlink_minimal_gpd = gpd.GeoDataFrame(starlink_minimal)
starlink_minimal_gpd = starlink_minimal_gpd.set_crs(epsg="4326")
starlink_minimal_gpd = starlink_minimal_gpd.to_crs(common_crs)


print(starlink_minimal)

# buffering the points outside of the generation function
# starlink_sample_buffered = starlink_sample.copy(deep=True)
# starlink_sample_buffered["geometry"] = starlink_sample_buffered["geometry"].combine(starlink_sample_buffered["Radius"],
#                                                                                     lambda x, y: x.buffer(y))
# print(starlink_sample_buffered.head())


# create a map projection that can be reasoned about
# and project it to the same map projection
land_filepath = dirpath / "map_data" / "ne_10m_land_scale_rank.zip"
ocean_filepath = dirpath / "map_data" / "ne_10m_ocean_scale_rank.zip"
land_map = gpd.read_file(filename=land_filepath)
ocean_map = gpd.read_file(filename=ocean_filepath)
land_map = land_map.to_crs(common_crs)
ocean_map = ocean_map.to_crs(common_crs)

# Tested that this map file does actually label the areas that are land
# so a spatial join will probably be able to get us land area
# print(world_map.head())
# print(starlink.head())


# plotting the map
fig, ax = plt.subplots()
# land_map.plot(ax=ax, color="#228B22")
# ocean_map.plot(ax=ax, color="#246BCE")

# plotting starlink data
# starlink.plot(ax=ax, color="#C51E3A", markersize = 1)
# starlink_sample.plot(ax=ax, color="#C51E3A", markersize = 1)
# starlink_sample_buffered.plot(ax=ax, color="#C51E3A", markersize = 1)
starlink_minimal_gpd.plot(ax=ax, color="#C51E3A")

# showing the plot
plt.show()
