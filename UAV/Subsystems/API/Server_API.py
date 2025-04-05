"""
27/03/2025 Zikos Ioannis: Created a basic API for the server host that is going to receive all the data (right now just some test URL curl requests)
                          from the Raspberry Pi API (Pi_API.py) and saves the URL at a txt file (logs.txt) for later use with the AI detection model
                          (the URLs are going to be replaced by photos from the camera)
"""

from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI()

# Define a data model for incoming requests
class SensorData(BaseModel):
    image_url: str

@app.get("/")
def read_root():
    return {"message": "FastAPI server is running!"}

@app.post("/data")
def receive_data(data: SensorData):
    print(f"Received data: {data}")  # Print received data in terminal
    logs = open("logs.txt", "a") # Open logs.txt and append to the end of the file
    logs.write(data.image_url + "\n") # Write line to log
    logs.close 
    return {"status": "success", "received": data}
    
# Run with: uvicorn Server_API:app --host 0.0.0.0 --port 8000 
#   With 0.0.0.0 it automatically scans the port for activity in the network, otherwise give global for WAN
# Check if running: curl -X GET "http://Servers_IP:8000/"

