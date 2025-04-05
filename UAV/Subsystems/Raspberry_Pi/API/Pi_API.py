"""
27/03/2025 Zikos Ioannis: Created a basic API for the Raspberry Pi that takes curl requests and sends and saves their data at the Server_API host server
                          (We will need to add a reverse ssh tunnel in the future for remote accessing. I am currently working on it)
"""

import falcon, requests, json
from ImageReqHeader import ImgReq 

class ImageHandler:
    def on_post(self, req, resp):
        try:
            data = req.media # Get JSON data from the request

            if "image_url" not in data:
                resp.status =  falcon.HTTP_400
                resp.media = {"error": "No image url provided!"}
                return
        
            ImgReq.set_url(data['image_url']) # Store URL
            resp.status = falcon.HTTP_200
            resp.media = {"error": "Image received successfully."}

        except Exception as e:
            resp.status = falcon.HTTP_500
            resp.media = {"error": f"Internal server error: {str(e)}"}

class DataSender:
    def on_get(self, req, resp):
        data = ImgReq.get_url() # Fetch stored image URL

        # Send data to FastAPI server
        ServerURL = "http://Servers_IP:8000/data"
        headers = {"Content-Type": "application/json"}
        response = requests.post(ServerURL, json=data, headers=headers)

        resp.media = {
            "message": "Sent image url to FastAPI",
            "server_response": response.json()
        }


app = falcon.App()
app.add_route('/image_req', ImageHandler()) # For storing image_url
app.add_route('/send_to_server', DataSender()) # For sending to FastAPI (For now its manual, will be automated)

# Run with: gunicorn -b Pis_IP:5000 Pi_API:app --workers 1
# Or for Windows testing run: python -m waitress --listen=0.0.0.0:5000 Pi_API:app
# Check: http://Pis_IP:5000/send_to_server

"""
The curl request should look like:
    curl -X POST "http://192.168.1.7:5000/image_req" -H "Content-Type: application/json" -d '{"image_url": "https://example.com/image.jpg"}'
"""
