#!/usr/bin/python
# -*- coding:utf-8 -*-
import RPi.GPIO as GPIO # Used for GPIO pins
import serial # Serial port communication
import time # Time lib
import os
from coord_convert import Convert # Coordinate conversion class

import sys # Lib for paths
libs_path = os.path.expanduser("~/Documents/GitHub/UAV/Raspberry_Pi/libs")
sys.path.append(libs_path)
from post_data import Post # Post data lib 

'''
Victor Meimaris 25/2/25: Created boot file for GPS module to run on boot and start sending coordinates to API and Ardupilot, gps now retains its last known location
even if signal is lost, this might change later because ardupilot is doing this automatically and this might produce errors. Timestamp 
will also be the same for the fast determination of time when the signal was lost, this might cahnge if ardupilot can't work with outdated timestamps.
'''

gps_on = 'AT+CGPS=1' # Start GPS session, =1,(1,2,3) for standalone, based, assisted mode, default is standalone
gps_off = 'AT+CGPS=0'
gps_info = 'AT+CGPSINFO'
gps_info_alt = 'AT+CGNSSINFO' # GNSS info 

url = 'http://ip:8000/gps'# GPS API url
'''
Initialize serial port and clear contents
Two-way communication with the GPS module is achieved through serial commands
'''
ser = serial.Serial('/dev/ttyS0',115200) # Serial port connection with the gps module and Baudrate of 115200.
ser.flushInput() # Clear serial input buffer of any leftover data.
# Begin GPS session.
ser.write((gps_on+'\r\n').encode())
time.sleep(5)
print("Begining GPS Session")
'''
Cold Start will begin after first reboot of module, ~1.5-2 minutes empirically required for satellite fix.
Restarting GPS session without resetting the module retains satellite fix data and jumps to hot start, ~30 second required for sat fix 
Cold start can be skipped(shortened?) by implementing Assisted-GNSS (AGNSS), which requires established server connection through the module in order
to download tracking data
'''
buff_raw = '' # Initialize buffer.Data sent from the module through the serial port will be stored here
API_obj = Post(url, buff_raw) # Create instance of class Post
temp = 'No GPS data yet'; # Print temp at boot until data arrives.
while(True):
    buff_raw = '' # Clear contents of buffer.
    ser.write((gps_info_alt+'\r\n').encode()) # Request data received through the GPS (location/time)
    time.sleep(1)
    if ser.inWaiting(): # Check if data has arrived in the serial buffer. inWaiting() returns number of bytes waiting
        buff_raw = ser.read(ser.inWaiting()).decode('utf-8', errors='ignore') # Copy serial output to local buffer
        buff_conv = Convert(buff_raw)# Create instance of class Convert, to convert to GPAA for Ardupilot
        buff = buff_conv.parse_cgnssinfo(buff_raw) # Converting to wanted format
        if ',,,,,,,,,,,,,,,' in buff_raw or buff_raw is None: # Validate this !!!!!!!!
            #print(temp)
            # Here is the line where we would send data to the ardupilot serially telling it "No fix" since gps has lost signal and buff = None
            if 'No GPS data' in temp: # The reason why we dont keep a temp for ardupilot is because it keeps the last location automatically and enters failsafe mode.
                API = API_obj.post_API(url, temp)
            else:
                API = API_obj.post_GPS(url, temp) # Posts to the API, function needs raw buffer input to work.
        else:
            #print(buff)
            # Here is the line where we would send gps data to the ardupilot serially, through the buff variable.
            API = API_obj.post_GPS(url, buff_raw)
            temp = buff_raw # As soon as buffer gets its first reading, whenever the signal is lost the gps will be displaying the last known location.
    else:
        print("GPS not functioning correctly...") # Error message shown if module is turned off or other error occurs, for example 