from pyorbital.orbital import Orbital
from datetime import datetime

# Define the TLE
tle_name = "STARLINK-31725"
tle_line1 = "1 59538C 24074E   24277.42965278  .00011426  00000+0  27260-3 0  2774"
tle_line2 = "2 59538  43.0026  60.9100 0001512 263.4205 201.9955 15.40689613    19"

# Create an Orbital object
satellite = Orbital(tle_name, line1=tle_line1, line2=tle_line2)

# Get current time in UTC
current_time = datetime.utcnow()

# Compute position and velocity
position = satellite.get_lonlatalt(current_time)

# Print the position and velocity with more details
print(f"Position (km) at {current_time}: {position}")

# NOTE: We should use a standard time otherwise the output will be different