#29/11/24 Victor Emmanuel Meimaris : Creating gps library that handles incoming POST requests.

import falcon,json # API libs

class gps_location: # Handles incoming POST requests with GPS data
    def on_post(self, req, resp):
        # Get JSON data from request body and convert to JSON obj automatically with req.media
        gps_data = req.media 
        if not gps_data:
            resp.status = falcon.HTTP_400
            resp.text = "gps_data is missing."
            return
        # Print GPS data to terminal
        print(f"Received GPS data: {gps_data}")

        # Response if communication was successful

        resp.status = falcon.HTTP_200
        resp.content_type = falcon.MEDIA_JSON
        resp.text = json.dumps({
            "status": "success",
            "message": "GPS data received successfully",
            "data": gps_data
        })