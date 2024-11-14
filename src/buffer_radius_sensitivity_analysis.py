import datetime
import geodataframe_processing as gdfp

# EPSG Codes
EQUAL_DISTANCE_EPSG = 4087
EQUAL_AREA_EPSG = 6933

# Global Parameters
BUFFER_RADIUS = 12065 # 121065
PERTURB_DISTANCE_KM = 50

# Time the TLE data was collected
TLE_TIME = datetime.datetime(2024, 11, 14, 16, 12, 32, 202983)


# Importing a static set of TLEs for sensitivity analysis
base_gdf = gdfp.load_satellites_from_file(EPSG=EQUAL_AREA_EPSG,
                                          input_filename="parameter_sensitivity_analysis_tles.txt",
                                          time=TLE_TIME)

