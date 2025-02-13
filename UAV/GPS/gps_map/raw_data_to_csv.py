#13/2/25 Victor Meimaris: Created script that turns raw gps data from log and turns it into CSV for testing purposes.

import csv # CSV lib to write on the CSV with actual data
import re # Lib for matching
import os # Lib for os.path.expand_user for easier path access.

import sys # Lib for paths
libs_path = os.path.expanduser("~/Documents/GitHub/UAV/useful_libs")
sys.path.append(libs_path)

from file_namer import FileNamer # Custom lib for renaming if named file already exists

# NMEA conversion to decimal degrees
def convert(nmea_coord, dir): #NMEA coordinates and direction
    deg = int(nmea_coord / 100)
    minutes = nmea_coord - deg * 100
    dec_deg = deg + minutes / 60
    if dir in ['S', 'W']:
        dec_deg = -dec_deg
    return dec_deg # return decimal degree after conversion.

# Input log file where gps data is stored after test.                            
input_file = os.path.expanduser("~/Documents/GitHub/UAV/GPS/gps_map/data/raw_logs/log_gps.txt")  # Replace with your log file path if script isn't on documents--> github.
output_file = os.path.expanduser("~/Documents/GitHub/UAV/GPS/gps_map/data/CSV/gps_raw.csv")  # Output CSV file to desired file path, change if you want.

file_namer = FileNamer(output_file)
output_file = file_namer.get_filename() # If output name already exists, changes it from gps_raw to gps_raw_0 for example

# Regex pattern to extract lattitude and longitude from raw log file
pattern = re.compile(r"^\+CGNSSINFO: (?:[^,]*,){4}\s*([\d.]+)\s*,\s*([NS])\s*,\s*([\d.]+)\s*,\s*([EW])")

with open(input_file, 'r') as log_file, open(output_file, 'w', newline='') as csvfile: # Open log file and create and write on gps_raw_0.csv.
    writer = csv.writer(csvfile)
    writer.writerow(['lat', 'lon'])  # Write CSV header.

     # !!! Change if format in regex if needed (if output format of gps is changed) !!!
    for line in log_file: # For each line in the log file we extract data using regex with the format below. 
        match = pattern.search(line) # Search for the pattern in each line
        if match:
            lat_nmea, lat_dir, lon_nmea, lon_dir = match.groups()
            lat_nmea, lon_nmea = float(lat_nmea), float(lon_nmea) # Storing nmea lattitude and longitude into lat_nmea, lon_nmea.

            # Converting coordinates to deccimal degrees
            lat = convert(lat_nmea, lat_dir)
            lon = convert(lon_nmea, lon_dir)

            # Write to CSV
            writer.writerow([lat, lon])