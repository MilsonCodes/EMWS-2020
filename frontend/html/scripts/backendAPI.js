/** 
 * backendAPI.js
 * -------------------------
 * The file for accessing the Python backend server to run calculations
 * 
*/
var backendAPI = backendAPI || {};

const hostname = window && window.location && window.location.hostname;

const API_HOST = hostname === "math.lsu.edu" ? "https://emws.pythonanywhere.com/" : "http://localhost:5000/";      //TBD

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

      console.log(res)
      break;
    } catch (e) {
      res = e.response
      numReqSent++
    }
  }

  return res.json();
}

const convertComplexMatrixToPythonParsableMatrix = matrix => {
  var newMat = new Array(matrix.length)

  for(var i = 0; i < matrix.length; i++) {
    newMat[i] = new Array(matrix[i].length)

    for(var j = 0; j < matrix[i].length; j++) {
      var val = matrix[i][j]

      if(val.re != 0 || val.im != 0) {
        if(val.re != 0 && val.im == 0)
          newMat[i][j] = val.re
        else
          newMat[i][j] = { re: val.re, im: val.im }
      } else {
        newMat[i][j] = 0
      }
    }
  }

  return newMat
}

const convertJSLayersToPythonLayers = layers => {
  var newLayers = []

  console.log(layers)

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

//const convert

backendAPI.createStructureObject = (omega, k1, k2, layers) => {
  return new Structure(omega, k1, k2, layers)
}

class Structure {
  constructor(omega, k1, k2, layers=[]) {
    this.omega = omega
    this.k1 = k1
    this.k2 = k2
    this.layers = layers
  }

  getOmega() {
    return this.omega
  }

  getk1() {
    return this.k1
  }

  getk2() {
    return this.k2
  }

  getLayers() {
    return this.layers
  }

  async buildStructure() {
    var data = { 
      omega: this.omega,
      k1: this.k1,
      k2: this.k2,
      layers: convertJSLayersToPythonLayers(this.layers)
    }

    var res = await request("structure/modes", data, "POST");

    this.maxwell_matrices = res.maxwell_matrices
    this.eigenvalues = res.eigenvalues
    this.eigenvectors = res.eigenvectors
  }

  async getConstantVector(incoming=[1, 0, 0, 0]) {
    var data = {
      omega: this.omega,
      k1: this.k1,
      k2: this.k2,
      layers: convertJSLayersToPythonLayers(this.layers)
    }
  }
}