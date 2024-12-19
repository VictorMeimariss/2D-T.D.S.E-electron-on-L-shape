#29/11/24 Victor Emmanuel Meimaris : Creating A.I model library for different modes (fire detection,rescue etc etc)

# Environmet setup
import os
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0" #This disables oneDNN otherwise terminal sends errors
HF_HUB_OFFLINE = 1 # No attempt to download from the HUB
HF_DATASETS_OFFLINE = 1 # Only local datasets

# Main AI Libraries
import falcon, json, torch # API libs and torch
from PIL import Image # Image lib
from transformers import pipeline, AutoImageProcessor, AutoModelForImageClassification  # This will be deleted as soon as custom pipeline is created
#from pipeline_placeholder import pipeline # Custom pipeline lib for ONNX and tensorflow lite that works without pi-torch for raspberry pi

#Initialise Processor locally
LocalProcessor = AutoImageProcessor.from_pretrained("C:/path/to/processor_file",local_files_only=True)

""" 
Fire detection A.I mode, Stage #1 -> Using most robust model tested by us and huggingface, if indication higher than certain number,
moving to Stage #2 -> Using all models and getting the average.
"""
class fire_detection:
    def __init__(self): # Setting up models

        self.model_path_stack = [] # Model stack

        # Defining path to models 
        self.model_path_stack.append("C:/path/to folder of model_0 not file because it won't work") 

        # Initialising models
        self.models = [AutoModelForImageClassification.from_pretrained(path, local_files_only=True) 
                       for path in self.model_path_stack]
        # Loading Models
        self.pipelines = [pipeline("image-classification", model = model, image_processor = LocalProcessor, framework="pt")
                          for model in self.models]
        # First Stage ready
        self.first_stage_method = self.pipelines[0]
        # Second stage method
        self.second_stage_method = self.pipelines

    def on_post(self, req, resp): # Handles POST requsts with image_directory data, opens image and returns results

        # Get JSON data from request body and convert to JSON obj automatically with req.media, then handles image data
        data = req.media 
        image_directory = data['image_directory']

        # Error Handling if image_directory isn't in the requests
        if not image_directory:
            resp.status = falcon.HTTP_400
            resp.text = "Image directory is missing."
            return

        # Opens image from directory and stores data in "image" if there is no error unless there is an exception
        try:
            image = Image.open(image_directory)
        except Exception as e:
            resp.status = falcon.HTTP_400
            resp.text = f"Error opening image: {str(e)}"
            return

        # Getting Results by using first stage image-classification (with error handling)
        try:
            results = self.first_stage_method(image)
            # Uses function next() to look for item 'label'=="Fire" and get the score's value, if not fire_results = 0
            fire_results = next((item['score'] for item in results if item['label'] == "Fire"), 0) #!!!!We need to fine tune our models to have Fire labels
        except Exception as e:
            resp.status = falcon.HTTP_500
            resp.text = f"Error on first stage_one detection: {str(e)}"
            return
        if(fire_results >= 0.3): # If fire_level >= 0.3 we go to stage_two, "Number will be changed later for optimal performance"
            sum_fire_levels = 0 # Initialising variables for calculation of average
            count = 0
            try:
                for pipeline in self.second_stage_method:# For loop runs second stage method
                    results = pipeline(image)
                    fire_level = next((item['score'] for item in results if item['label'] == "Fire"), 0)
                    sum_fire_levels += fire_level
                    count += 1
                    if count == 0:
                        raise ValueError("No valid fire levels detected from any model.")
                average_fire_level = sum_fire_levels / count # Calculating average of fire_levels
                resp.text = json.dumps({"Average Fire level": average_fire_level}) # Returns result average of all models
            except Exception as e: # Error Handling
                resp.status = falcon.HTTP_500
                resp.text = f"Error on first stage_two detection: {str(e)}"
                return
        else:
            resp.text = json.dumps({"Fire level": fire_results}) # Returns results from model_0