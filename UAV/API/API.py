"""
29/11/24 Victor Emmanuel Meimaris : Creating Falcon API for Pi to be able to send and receive information,for example gps location,
a.i module results, updating code etc etc. We are using Falcon because of it's resourcefullness, to use as little processing power
as possible. In order to run the program you will have to install the requirements first.

26/02/25 Ioannis Zikos: Now we create instances of the YOLOv9_Detection() and ImageRequest() classes and we can use curl to give 
                        any image url we want to the model. New requirements are needed (I will post my current venv packages and versions 
                        although i do not have the neccessary packages for the gps)

"""

#!/usr/bin/python # Path to /python to be run easily from pi
# -*- coding:utf-8 -*- # Encoding declaration

import falcon,json # API libs
#from AI_API import fire_detection # AI lib
from AI_API import ImageRequest, YOLOv9_Detection # For requesting an image with c-url and runs the YOLO obj detection model
from GPS_API import gps_location # GPS lib
from wsgiref.simple_server import make_server # Running server directly

# Creating Falcon App, gps, a.i etc "modes" and routes

app = falcon.App() 

# Modes

#detection = fire_detection()
detection = YOLOv9_Detection()
image_req = ImageRequest()
gps = gps_location()

# Routes

app.add_route('/gps', gps) 
app.add_route('/detection', detection)
app.add_route('/image_req', image_req)

# Run the server directly by running program, for debugging

if __name__ == '__main__':
    print("Server running on http://localhost:8000")
    with make_server('localhost', 8000, app) as httpd:
        httpd.serve_forever()
