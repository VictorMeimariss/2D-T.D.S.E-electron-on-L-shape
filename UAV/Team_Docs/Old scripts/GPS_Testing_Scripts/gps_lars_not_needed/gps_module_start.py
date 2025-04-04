#!/usr/bin/python 
# -*- coding:utf-8 -*-
# First line Path to /python to be run easily from pi and second is encoding declaration.
import RPi.GPIO as GPIO # Used for GPIO pins
import serial # Serial port communication
import time # Time lib
from lars import Lars_Predict # LARS 
'''
5/12/24 Orestis Keisoglu making:
PROTOTYPE of a GPS session launch, run and termination for the sole purpose of rapid TESTING AT commands and GPS accuracy performance.
Below code is as of now extremely rudimentary and does not include error correction or any means of failsafe.

6/12/24 Victor Meimaris:
Added some comments, added .decode('utf-8', errors='ignore') in line 48, GNSS support command, potential commands that might be helpful at 17.2.12, 17.2.20
first is to delete gps information, second is used to respond to incoming positioning request. Lastly 17.2.8 to configure NMEA sentence type, we might be getting 
information, only from certain sattelites.

10/2/25 Victor Meimaris:
Added Least Angle Regression (LARS) for optimised accuracy in while loop after cold start.
'''
gps_on = 'AT+CGPS=1' # Start GPS session, =1,(1,2,3) for standalone, based, assisted mode, default is standalone
gps_off = 'AT+CGPS=0' 
gps_info = 'AT+CGPSINFO'
gps_set_acc = 'AT+CGPSHOR=5' # Accuracy of 0-1800000 meters, now set to 5
gps_sat = 'AT+CGNSSMODE=15,1' # Bit 0 - GLONASS Bit 1 - BEIDOU Bit 2 - GALILEO Bit 3 - QZSS, 1-enable 0-disable 15,1 ==> 1111 = 15 in Binary, using all sattelites 15,1 enable all sattelites

'''
Initialize serial port and clear contents
Two-way communication with the GPS module is achieved through serial commands
'''
ser = serial.Serial('/dev/ttyS0',115200) # Baud rate of 115200
ser.flushInput() # Clear serial input buffer of any leftover data 

ser.write((gps_off+'\r\n').encode()) # Shutdown GPS session, documentation recommends before cold start
print("Resetting GPS")
time.sleep(2) # Wait for 2 seconds empirically, might change later

# Config settingss
ser.write((gps_set_acc+'\r\n').encode()) # Set desired GPS accuracy
print("Setting Accuracy")
time.sleep(2)
ser.write((gps_sat+'\r\n').encode()) # Enabling GNSS Support
print("Enabling GNSS support from selected sattelites")
time.sleep(2)

# Begin GPS session
ser.write((gps_on+'\r\n').encode())
time.sleep(2)
print("Beginning GPS Session")
'''
Cold Start will begin after first reboot of module, ~10 minutes empirically required for satellite fix.
Restarting GPS session without resetting the module retains satellite fix data and jumps to hot start, ~30 second required for sat fix 
Cold start can be skipped(shortened?) by implementing Assisted-GNSS (AGNSS), which requires established server connection through the module in order
to download tracking data
'''
try:
    while(True):
        buff = ''# Initialize and clear contents of buffer.Data sent from the module through the serial port will be stored here
        ser.write((gps_info+'\r\n').encode()) # Request data received through the GPS (location/time)
        time.sleep(0.2)
        if ser.inWaiting(): # Check if data has arrived in the serial buffer. inWaiting() returns number of bytes waiting
            buff = ser.read(ser.inWaiting()).decode('utf-8', errors='ignore') # Copy serial output to local buffer
            Lars_Predict(buff) # Using LARS to achieve optimal accuracy
        else:
            print("GPS not ready yet...")
except KeyboardInterrupt: # Catch Ctrl+C in order to terminate program
    print("GPS Session Terminated")
    ser.write((gps_off+'\r\n').encode())# Shutdown GPS session, redundant to line 21
    time.sleep(2)
    exit()