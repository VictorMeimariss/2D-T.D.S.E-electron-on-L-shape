# 10/2/25 Victor Meimaris: 
# Created Least Angle Regression (LARS) model for optimised accuracy

import numpy as np
import re # Re lib for searching
from joblib import load # Lib to import trained model
from sklearn.linear_model import Lars

lars = load('lars.joblib') # Load trained model after executing lars_training.py
class Lars_Predict:
    def __init__(self, buff): # Initialise lars_predict
        self.buff = buff
    def corrected_coordinates(self):
        match = re.search(r"\+CGPSINFO:\s*(\d+\.\d+),([NS]),(\d+\.\d+),([EW])", self.buff) # Extract data from buffer with re.search
        if match:
            lat_nmea, lat_dir, lon_nmea, lon_dir = match.groups()
            lat_nmea, lon_nmea = float(lat_nmea), float(lon_nmea)
        gps_data = np.array([[lat_nmea, lon_nmea]]) # Extracted lattitude and longitude in DMM format
        corrected_data = lars.predict(gps_data)[0] # Correcting data with LARS
        corr_lat_nmea= corrected_data[0]
        corr_lon_nmea= corrected_data[1]
        return f"Enhanced GPS Data: Latitude:{corr_lat_nmea}, Longitude:{corr_lon_nmea}" # Prediction of actual coordinates based on old data 