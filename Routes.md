# Routes #

A list of all routes to run calculations from the backend.

*<backend_link>* = *http://some-link-to-backend.tbd/api*

### Modes ###
Route *<backend_link>/structure/modes*\
Method Type: **POST**\
POST Body:
```JSON
{
  "omega": 0.368,
  "k1": 0.5,
  "k2": 0.22,
  "layers": [
    {
      "name": "Ambient Left",
      "length": 10,
      "epsilon": [],
      "mu": []
    }
  ]
}
```
Response Status Code(s): 201 Created, 400 Bad Request, 500 Server Error\
Response Body:
```JSON
{
  "maxwell_matrices": [
    []
  ],
  "eigenvalues": [
    []
  ],
  "eigenvectors": [
    []
  ]
}
```

