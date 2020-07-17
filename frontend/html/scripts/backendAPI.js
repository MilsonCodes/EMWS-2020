/** 
 * backendAPI.js
 * -------------------------
 * The file for accessing the Python backend server to run calculations
 * 
*/
var backendAPI = backendAPI || {};

const hostname = window && window.location && window.location.hostname;

//const API_HOST = hostname === "math.lsu.edu" ? "https://emws.pythonanywhere.com/" : "http://localhost:5000/";      //TBD
const API_HOST = "https://emws.pythonanywhere.com/"

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
      console.log(e.response)
      res = e.response
      numReqSent++
    }
  }

  return res.json();
}

const simplifyJSComplex = val => {
  if(val.re != 0 || val.im != 0) {
    if(val.re != 0 && val.im == 0)
      return val.re
    else
      return { re: val.re, im: val.im }
  } else {
    return 0
  }
}

const convertEigenvaluesToPythonParsable = eigVals => {
  var newArr = new Array(eigVals.length)

  console.log(newArr)

  for(var i = 0; i < eigVals.length; i++) {
    newArr[i] = new Array(eigVals[i].length)

    for(var j = 0; j < eigVals[i].length; j++) {
      newArr[i][j] = simplifyJSComplex(eigVals[i][j])
    }
  }

  return newArr
}

const convertComplexMatrixToPythonParsableMatrix = matrix => {
  var newMat = new Array(matrix.length)

  for(var i = 0; i < matrix.length; i++) {
    newMat[i] = new Array(matrix[i].length)

    for(var j = 0; j < matrix[i].length; j++) {
      newMat[i][j] = simplifyJSComplex(matrix[i][j])
    }
  }

  return newMat
}

const convertMatrixArrayToPythonParsable = matrixArr => {
  var newArr = new Array(matrixArr.length)

  for(var i = 0; i < matrixArr.length; i++) {
    newArr[i] = convertComplexMatrixToPythonParsableMatrix(matrixArr[i])
  }

  return newArr
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

const parseStringVal = str => {
  var val = 1

  try {
    val = math.complex(str)
  } catch(e) {
    console.log("Could not parse string inputted -- reverting to 1")
  }

  return val
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
    this.incoming = [1,0,0,0]
  }

  updateValues(omega, k1, k2, layers) {
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

  getMaxwell() {
    return this.maxwell_matrices ? this.maxwell_matrices : null
  }

  getEigenvalues() {
    return this.eigenvalues ? this.eigenvalues : null
  }

  getEigenvectors() {
    return this.eigenvectors ? this.eigenvectors : null
  }

  async getConstantVector() {
    var data = {
      omega: this.omega,
      k1: this.k1,
      k2: this.k2,
      layers: convertJSLayersToPythonLayers(this.layers),
      incoming: this.incoming
    }

    if(this.eigenvalues && this.eigenvectors && this.maxwell_matrices) {
      data = {
        ...data,
        maxwell_matrices: convertMatrixArrayToPythonParsable(this.maxwell_matrices),
        eigenvalues: convertEigenvaluesToPythonParsable(this.eigenvalues),
        eigenvectors: convertMatrixArrayToPythonParsable(this.eigenvectors)
      }
    }

    var res = await request("structure/constants", data, "POST")

    if(!this.eigenvalues && !this.eigenvectors && !this.maxwell_matrices) {
      this.eigenvalues = res.eigenvalues
      this.eigenvectors = res.eigenvectors
      this.maxwell_matrices = res.maxwell_matrices
    }

    this.scattering = res.scattering
    this.constants = res.constants
  }

  setIncoming(incoming) {
    if(incoming) {
      var parsedIncoming = new Array(incoming.length)

      for(var i = 0; i < incoming.length; i++) {
        parsedIncoming[i] = simplifyJSComplex(parseStringVal(incoming[i]))
      }

      this.incoming = parsedIncoming
    }
  }

  getIncoming() {
    return this.incoming
  }

  getScatteringMatrix() {
    return this.scattering ? this.scattering : null
  }

  getConstantsVector() {
    return this.constants ? this.constants : null
  }

  async determineField() {
    var data = {
      omega: this.omega,
      k1: this.k1,
      k2: this.k2,
      layers: convertJSLayersToPythonLayers(this.layers),
      incoming: this.incoming
    }

    if(this.eigenvalues && this.eigenvectors && this.maxwell_matrices) {
      data = {
        ...data,
        maxwell_matrices: convertMatrixArrayToPythonParsable(this.maxwell_matrices),
        eigenvalues: convertEigenvaluesToPythonParsable(this.eigenvalues),
        eigenvectors: convertMatrixArrayToPythonParsable(this.eigenvectors)
      }
    }

    console.log(data)

    var res = await request("structure/field", data, "POST")

    console.log(res)

    if(!this.eigenvalues && !this.eigenvectors && !this.maxwell_matrices) {
      this.eigenvalues = res.eigenvalues
      this.eigenvectors = res.eigenvectors
      this.maxwell_matrices = res.maxwell_matrices
    }

    this.scattering = res.scattering
    this.constants = res.constants
    this.field = res.field
  }

  getField() {
    return this.field ? this.field : null
  }
}