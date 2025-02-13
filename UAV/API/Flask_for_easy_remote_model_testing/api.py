from flask import Flask, request, jsonify
from flask_restful import Api  # Import restful api libs 
from ai_model import fire_detect  # Ensure fire_detect is properly defined in ai_model

app = Flask(__name__)
api = Api(app)
detector = fire_detect()

@app.route('/detection', methods=['POST'])
def detection():
    data = request.get_json()
    if 'image_url' not in data:
        return jsonify({"error": "No image url was provided"}), 400
    
    image_url = data['image_url']
    
    
    try:
        # Call the fire detection function with the image URL
        results = detector.detect(image_url)
        return jsonify(results), 200 
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)
