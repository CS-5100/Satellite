# trying to use and work through techniques outlined in Python for Geographic Data Analysis
# (https://python-gis-book.readthedocs.io/en/develop/index.html)

import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timezone
import pyproj as pyj
import geopandas as gpd
import download_tle as dtl

# defining the data directory as a relative path
dirpath = Path(__file__).parent.resolve() / ".." / "data"

# getting the current time
current_time = datetime.now(timezone.utc)

# creating a DataFrame object
starlink_current = dtl.download_current_tles_to_dataframe()
starlink_gpd = gpd.GeoDataFrame(data=starlink_current,
                                geometry=gpd.points_from_xy(x=starlink_current["Longitude"],
                                                            y=starlink_current["Latitude"]),
                                crs="EPSG:4326")

# importing map data
land_filepath = dirpath / "map_data" / "ne_10m_land_scale_rank.zip"
ocean_filepath = dirpath / "map_data" / "ne_10m_ocean_scale_rank.zip"
land_map = gpd.read_file(filename=land_filepath)
ocean_map = gpd.read_file(filename=ocean_filepath)

# plotting the map
fig, ax = plt.subplots()
land_map.plot(ax=ax, color="#228B22")
ocean_map.plot(ax=ax, color="#246BCE")
starlink_gpd.plot(ax=ax, color="#C51E3A", markersize = 1)
plt.show()

# TODO: figure out projection to cartesian coordinate system then buffer
