# trying to use and work through techniques outlined in Python for Geographic Data Analysis
# (https://python-gis-book.readthedocs.io/en/develop/index.html)

import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import geopandas as gpd
import download_tle as dtl

# defining a few constants
EARTH_SURFACE_AREA_SQ_KM = 509600000
EARTH_LAND_AREA_SQ_KM = 148326000

# defining the data directory as a relative path
dirpath = Path(__file__).parent.resolve() / ".." / "data"

# getting the current time
current_time = datetime.now(timezone.utc)

# creating a GeoDataFrame object
starlink_current = dtl.download_current_tles_to_dataframe()
starlink_gpd = gpd.GeoDataFrame(data=starlink_current,
                                geometry=gpd.points_from_xy(x=starlink_current["Longitude"],
                                                            y=starlink_current["Latitude"]),
                                crs="EPSG:4326")
starlink_gpd_points = starlink_gpd.copy(deep=True)

# importing map data
land_filepath = dirpath / "map_data" / "ne_10m_land_scale_rank.zip"
ocean_filepath = dirpath / "map_data" / "ne_10m_ocean_scale_rank.zip"
land_map = gpd.read_file(filename=land_filepath)
ocean_map = gpd.read_file(filename=ocean_filepath)

# defining a couple of epsg codes
equal_distance_epsg = 4087
equal_area_epsg = 6933

# re-projecting the starlink gpd data onto the equal distance
# projection to get all the points within a particular distance
# of each satellite point
starlink_gpd = starlink_gpd.to_crs(epsg=equal_distance_epsg)

# units of distance are in meters, 24.13 km diameter = 12.065 km radius, from satnav website
starlink_gpd['geometry'] = starlink_gpd['geometry'].buffer(12065)

# Re-projection onto equal area map
land_map = land_map.to_crs(epsg=equal_area_epsg)
ocean_map = ocean_map.to_crs(epsg=equal_area_epsg)
starlink_gpd = starlink_gpd.to_crs(epsg=equal_area_epsg)
starlink_gpd_points = starlink_gpd_points.to_crs(epsg=equal_area_epsg)

# filtering out any invalid geometries
starlink_valid = starlink_gpd[starlink_gpd['geometry'].is_valid]

# plotting the map
fig, ax = plt.subplots()
land_map.plot(ax=ax, color="#228B22")
ocean_map.plot(ax=ax, color="#246BCE")
starlink_valid.plot(ax=ax, color="#C51E3A")
# starlink_gpd_points.plot(ax=ax, color="#C51E3A", markersize = 1)
plt.show()

# getting the total area covered by all satellites in square meters
starlink_all_areas = starlink_valid['geometry'].union_all()
land_all = land_map['geometry'].union_all()
ocean_all = ocean_map['geometry'].union_all()

starlink_all_areas_square_km = starlink_all_areas.area / (1000**2)
land_all_areas_square_km = land_all.area / (1000**2)
ocean_all_areas_square_km = ocean_all.area / (1000**2)
earth_area = land_all_areas_square_km + ocean_all_areas_square_km

print(starlink_all_areas_square_km)
print(land_all_areas_square_km)
print(EARTH_LAND_AREA_SQ_KM)
print(earth_area)
print(EARTH_SURFACE_AREA_SQ_KM)

# land area comparison doesn't make sense given we need intersection logic
# can't do this for land area yet
print(starlink_all_areas_square_km / earth_area) 