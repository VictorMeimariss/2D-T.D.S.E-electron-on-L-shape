import os
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0" #This disables oneDNN otherwise terminal sends messages about disabling oneDNN.

import requests
from PIL import Image
from transformers import pipeline


class fire_detect:
    def __init__(self, model_name="EdBianchi/vit-fire-detection"):
        self.pipe = pipeline("image-classification", model=model_name, framework="pt")#Initialise model
    def detect(self,image_url):
        image = Image.open(requests.get(image_url, stream=True).raw)
        results=self.pipe(image)
        return results
