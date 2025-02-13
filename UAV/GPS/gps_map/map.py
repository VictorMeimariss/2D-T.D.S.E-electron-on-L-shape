""" 
12/2/25 Victor Meimaris creating map that shows actual and gps coordinates with cubic splines (if uncommented) on a map to visualise error.
"""
import numpy as np #numbers library for array and function use.
import folium # Library for map creation.
import pandas as pd # Pandas library for reading gps data from data CSV
from scipy.interpolate import CubicSpline # Import cubic spline to interpolate function from actual data on the map.
import os # Lib for os.path.expand_user for easier path access.
import sys # Lib for paths
libs_path = os.path.expanduser("~/Documents/GitHub/UAV/useful_libs")
sys.path.append(libs_path)

from file_namer import FileNamer # Custom lib for renaming if named file already exists

columns = ['lat', 'lon'] # Columns of data in data CSV.
file_path = os.path.expanduser("~/Documents/GitHub/UAV/GPS/gps_map/data/CSV/gps_raw.csv") # Change if github folder not on documents.
data = pd.read_csv(file_path, usecols=columns) # Read data from csv file

lat = data['lat'].values # Latitude
lon = data['lon'].values # Longitude

# Cubic spline interpolation method for actual data which is known, creating functions with CupicSpline
lat_func = CubicSpline(np.arange(len(lat)),lat) 
lon_func = CubicSpline(np.arange(len(lon)),lon)
"""
# Least squares script, too much loss of accuracy for real data, reccomended only for faulty raw gps data.
# Acquire function coefficients for the least squares polynomial "fitting to minimise difference betwen curve and actual points".
degree = 2 # 2 --> Least squares polynomial, for more accuracy could use least cubes degree= 3 if gps loses accuracy on curves.
lat_coeffs = np.polyfit(np.arange(len(lat)), lat, degree)  # Fit latitude --> Creating function coefficients with np.polyfit .
lon_coeffs = np.polyfit(np.arange(len(lon)), lon, degree)  # Fit longitude             np.arange returns every value spaced out

# np.poly1d creates functions out of the coefficients                                  polyfit(num,array,degree)
lat_func = np.poly1d(lat_coeffs)
lon_func = np.poly1d(lon_coeffs)
"""
""" # !!!!!!!! Uncomment if you want the spline, reccomended for actual data points for example from map_datactual.csv which is from google maps.
# Points along the fitted curve
num_points = 500  # Number of points to generate
t = np.linspace(0, len(lat) - 1, num_points)  # Curve parameter linspace like in matlab np.linspace generates 100 points from 0 (len(lat)-1) to plot the function on them.
lat_fitted = lat_func(t)  # Fitted latitude points
lon_fitted = lon_func(t)  # Fitted longitude points

fitted_points = list(zip(lat_fitted, lon_fitted)) # fitted point list of (lan,lon) couples 
# I am using zip because folium expects the wanted points as tuples within a list.
# Also using a list because tuples cannot be fit in np.array because they aren't numbers.'''
"""
# Map creation with folium
m = folium.Map( 
    location=[lat[0], lon[0]], # Center of map
    zoom_start=14, # Starting zoom
    tiles='https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',  # Hybrid tiles for both sattelite photo as well as street names
    attr='Google'  # Attribution
)
""" # !!!!!!!! Uncomment if you want the spline
# Polyline on top of the map
folium.PolyLine(
    locations= fitted_points,  # List of (lat, lon) points
    color='blue',  # Color of the line connecting these points
    weight=5,  # Line thickness 
    opacity=0.8  # Line transparency
).add_to(m) # Addition of the Polyline item to the map.
"""
# Add markers for the original GPS points.
for point in zip(lat, lon): # For each actual point in zip(lan,lon) (again using zip because of folium) add a marker.
    folium.Marker(point).add_to(m)

output_file = os.path.expanduser("~/Documents/GitHub/UAV/GPS/gps_map/data/maps/map.html")
file_namer = FileNamer(output_file) 
output_file = file_namer.get_filename() # Unique map name 
m.save(output_file) # Saving map to an HTML file on data folder.