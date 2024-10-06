from pyorbital.orbital import Orbital
from datetime import datetime

def process_tles(file_name, output_file_name):
    # Open the TLE file for reading
    with open(file_name, 'r') as tle_file:
        lines = tle_file.readlines()
    
    # Prepare the output file to write results
    with open(output_file_name, 'w') as output_file:
        output_file.write("Satellite, Latitude, Longitude, Altitude\n")
        
        # Process the TLEs (each TLE consists of 3 lines: name, line1, line2)
        for i in range(0, len(lines), 3):
            tle_name = lines[i].strip()          # Satellite name
            tle_line1 = lines[i+1].strip()       # Line 1 of the TLE
            tle_line2 = lines[i+2].strip()       # Line 2 of the TLE
            
            # Create an Orbital object
            satellite = Orbital(tle_name, line1=tle_line1, line2=tle_line2)
            
            # Get the current time in UTC
            current_time = datetime.utcnow()
            
            # Calculate the latitude, longitude, and altitude
            lon, lat, alt = satellite.get_lonlatalt(current_time)
            
            # Write the results to the output file
            output_file.write(f"{tle_name}, {lat:.6f}, {lon:.6f}, {alt:.2f}\n")
    
    print(f"Lat, Long, Alt data saved to {output_file_name}")

# Specify the file containing TLEs and the output file
# NOTE: NEED TO UPDATE PATH BASED ON WHERE YOU SAVE THE FILES
process_tles(r'C:\Users\chaud\Desktop\cs5100\starlinkTLEs.txt', 'starlink_positions.csv')
