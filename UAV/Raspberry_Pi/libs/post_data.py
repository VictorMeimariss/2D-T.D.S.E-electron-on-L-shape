import json
import requests
import re

# Victor Meimaris 26/2/25: Creating lib for posting data to the API or Ardupilot
class Post:
    def __init__(self, url, payload):
        self.url = url
        self.payload = payload

    # Function to post data directly to the API
    def post_API(self, url, payload):
        response = requests.post(url, json=payload) # json= is like using json.dumps(payload), so there will be no need for conversion.
        #print(response.status_code, response.json())
        return response
    
    # NMEA conversion to decimal degrees
    @staticmethod # Convert does not use instance parameters 
    def convert(nmea_coord, dir):  # NMEA coordinates and direction
        deg = int(nmea_coord / 100)
        minutes = nmea_coord - deg * 100
        dec_deg = deg + minutes / 60
        if dir in ['S', 'W']:
            dec_deg = -dec_deg
        return dec_deg  # return decimal degree after conversion.
    
    # Function to post GPS data to the API
    def post_GPS(self, url, payload):  
        pattern = re.compile(r"^\+CGNSSINFO: (?:[^,]*,){4}\s*([\d.]+)\s*,\s*([NS])\s*,\s*([\d.]+)\s*,\s*([EW])\s*,\s*([\d]+)\s*,\s*([\d.]+)")
        match = pattern.search(payload)
        if match:
            lat_nmea, lat_dir, lon_nmea, lon_dir, date, timestamp = match.groups()  # Storing data into variables
            lat_nmea, lon_nmea, timestamp = float(lat_nmea), float(lon_nmea), float(timestamp)

            # Converting coordinates to decimal degrees
            lat = Post.convert(lat_nmea, lat_dir) # Static methods need Post.
            lon = Post.convert(lon_nmea, lon_dir)

            # Storing data in a JSON
            data = {
                'Latitude': lat,
                'Longitude': lon,
                'Timestamp': timestamp,
                'Date': date
            }

            response = self.post_API(url, data)
            return response
        else:
            return None
