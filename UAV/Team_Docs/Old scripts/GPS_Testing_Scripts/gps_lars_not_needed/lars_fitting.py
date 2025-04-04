# 10/2/25 Victor Meimaris: 
# Fitting Least Angle Regression (LARS) for optimised accuracy, this is the training file so we won't have to train it everytime
# Needs work!!!!!
import numpy as np
from joblib import dump # To dump training into a file
from sklearn.linear_model import Lars
import pandas as pd # Pandas library for reading gps data from map_data.csv

# Setting up falty and actual coordinates in NMEA DMM format: 4108.445,N,02453.67190,E ->4108.445 & 2453.67190 lat and long
# Actual GPS data from google maps (path taken)
columns = ['lat_a', 'lon_a'] # Columns of data in map_data.csv
file_path = 'C:/Users/Victor Meimaris/Documents/GitHub/UAV/GPS_Module/map_data.csv' # Change if on other laptop/pc
data = pd.read_csv(file_path, usecols=columns) 
lat = data['lat'].values 
lon = data['lon'].values

X_train = np.array([[latf1, lotf1],[latf2, lonf2],[latf3, lonf3]]) # Falty latitude and longitude coordinates taken from gps module
Y_train = np.array([[lata1, lota1],[lata2, lona2],[lata3, lona3]]) # Actual latitude and longitude coordinates

lars = Lars().fit(X_train, Y_train) # Train model based on collected data 
dump(lars, 'lars.joblib') # Dump trained model into a file
print("Model trained and saved to 'lars_model.joblib'.")