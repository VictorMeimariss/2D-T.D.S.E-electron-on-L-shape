from transformers import AutoModel, AutoImageProcessor

# Specify the model name
model_name = "EdBianchi/vit-fire-detection"

# Load the model and the image processor
model = AutoModel.from_pretrained(model_name, cache_dir="./models")
image_processor = AutoImageProcessor.from_pretrained(model_name, cache_dir="./models")

print("Model and image processor loaded successfully!")

