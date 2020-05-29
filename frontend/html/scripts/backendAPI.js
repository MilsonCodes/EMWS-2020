/** 
 * backendAPI.js
 * -------------------------
 * The file for accessing the Python backend server to run calculations
 * 
*/
var backendAPI = backendAPI || {};

const API_HOST = "http://localhost:5000/";      //TBD

const request = async (route, body, request_type="POST", content_type="application/json") => {
  var numReqSent = 0, maxReqAttempt = 5, res = null

  while(numReqSent < maxReqAttempt) {
    try {
      res = await fetch(API_HOST + route, {
        method: request_type,
        body: JSON.stringify(body),
        headers: {
          'Content-Type': content_type
        }
      })
    } catch (e) {
      res = e.response
      numReqSent++
    }
  }

  return res.json();
}

const convertComplexMatrixToPythonParsableMatrix = matrix => {
  var newMat = matrix

  for(var i = 0; i < matrix.length; i++) {
    for(var j = 0; j < matrix[i].length; i++) {
      var val = matrix[i][j]

      if(val.re || val.im)
        newMat[i][j] = { re: val.re, im: val.im }
      else
        newMat[i][j] = val
    }
  }

  return newMat
}

const convertJSLayersToPythonLayers = layers => {
  var newLayers = []

  for(var i = 0; i < layers.length; i++) {
    var curLayer = layers[i]

    newLayers[i] = {
      name: curLayer.layerName,
      length: curLayer.length,
      epsilon: convertComplexMatrixToPythonParsableMatrix(curLayer.epsilonA),
      mu: convertComplexMatrixToPythonParsableMatrix(curLayer.muA)
    }
  }

  return newLayers
}

backendAPI.buildStructure = async (omega, k1, k2, layers) => {
  var data = {
    omega, k1, k2, layers: convertJSLayersToPythonLayers(layers)
  }

  var res = await request("structure/modes", data);

  return res
}