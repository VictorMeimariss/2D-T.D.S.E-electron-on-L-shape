Easy use API for model testing-->

1)Use open cmd within the folder and type"pip install -r requirements.txt".
2)To open server type"python api.py".
3)Then open curl with git bash and paste that command along with your url


curl -X POST http://127.0.0.1:5000/detection -H "Content-Type: application/json" -d '{"image_url": "---TYPE--URL--HERE---"}'
curl -X POST http://127.0.0.1:5000/detection -H "Content-Type: application/json" -d "{\"image_url\": \"C:/Users/victo/Downloads/fire.jpg\"}"