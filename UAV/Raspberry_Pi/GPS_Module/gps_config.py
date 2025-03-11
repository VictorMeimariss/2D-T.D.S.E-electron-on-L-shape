#!/usr/bin/python
# -*- coding:utf-8 -*-
import RPi.GPIO as GPIO # Used for GPIO pins
import serial # Serial port communication
import time # Time lib

'''Victor Meimaris 25/2/25:Created seperate file for configuration settings to decreace boot time'''

gps_on = 'AT+CGPS=1' # Start GPS session, =1,(1,2,3) for standalone, based, assisted mode, default is standalone
gps_off = 'AT+CGPS=0'
gps_set_sentence = 'AT+CGPSNMEA=200191'
gps_set_rate = 'AT+CGPSNMEARATE=1'# 1 = 10Hz and 0 = 1Hz
gps_set_acc = 'AT+CGPSHOR=5' # Accuracy of 0-1800000 meters, now set to 5
gps_sat = 'AT+CGNSSMODE=15,1' # Bit 0 - GLONASS Bit 1 - BEIDOU Bit 2 - GALILEO Bit 3 - QZSS, 1-enable 0-disable 15,1 ==> 1111 = 15 in Binary, using all sattelites 15,1 enable all sattelites

'''
Initialize serial port and clear contents
Two-way communication with the GPS module is achieved through serial commands
'''
ser = serial.Serial('/dev/ttyS0',115200) # Baud rate of 115200 may change to lower baud rate like 9600
ser.flushInput() # Clear serial input buffer of any leftover data

ser.write((gps_off+'\r\n').encode()) # Shutdown GPS session, documentation recommends before cold start
print("Resetting GPS")
time.sleep(2) # Wait for 2 seconds empirically, might change later

# Config settingss
ser.write((gps_set_acc+'\r\n').encode()) # Set desired GPS accuracy
print("Setting Accuracy")
time.sleep(2)
ser.write((gps_set_rate+'\r\n').encode())
print("Setting rate")
time.sleep(2)
ser.write((gps_sat+'\r\n').encode()) # Enabling GNSS Support
print("Enabling GNSS support from selected sattelites")
time.sleep(2)
ser.write((gps_set_sentence+'\r\n').encode())
time.sleep(2)

# Begin GPS session
ser.write((gps_on+'\r\n').encode())
time.sleep(5)
print("New GPS Session has begun")