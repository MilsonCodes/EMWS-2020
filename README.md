# Electromagnetic Wave Scattering
### 2020 Update
Developers: [Noah Templet](https://github.com/w3aseL); [Michael Sheppard III](https://github.com/MilsonCodes)\
Research Supervisor: [Stephen Shipman](https://www.math.lsu.edu/~shipman)

[Backend Deployment Instructions](https://docs.google.com/document/d/1Hm3J-an80XWGBYOH-XYcGJwbcKlsgqoGKKTIXb3xzgI/edit?usp=sharing)

### Starting the App (Frontend)
Best option is to use Python. Creates a HTTP server in the frontend directory.

Simply just run this:
```
cd frontend && python -m http.server
```

### Starting the App (Backend)
Python is required here.

Install the following dependencies using pip:
 - flask
 - flask_cors
 - numpy
 - scipy

Simply just run this:
```
cd backend && python -m api.app
```

