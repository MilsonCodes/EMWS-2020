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

  //console.log(newArr)

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

  //console.log(layers)

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

const orderEigenvalues = (eVals, orderLeft, orderRight) => {
  var newEigVals = new Array(eVals.length)

  for(let i = 0; i < newEigVals.length; i++) {
    newEigVals[i] = new Array(eVals[i].length)

    for(let j = 0; j < newEigVals[i].length; j++) {
      if(i == 0) {
        var indexToGrab = orderLeft[j], val = eVals[i][indexToGrab]

        newEigVals[i][j] = val
      } else if(i == newEigVals.length - 1) {
        var indexToGrab = orderRight[j], val = eVals[i][indexToGrab]

        newEigVals[i][j] = val
      } else {
        newEigVals[i][j] = eVals[i][j]
      }
    }
  }

  console.log(newEigVals)

  return newEigVals
}

const orderEigenvectors = (eVecs, orderLeft, orderRight) => {
  var newEigVecs = new Array(eVecs.length)

  for(let i = 0; i < newEigVecs.length; i++) {
    newEigVecs[i] = new Array(eVecs[i].length)

    for(let j = 0; j < newEigVecs[i].length; j++) {
      if(i == 0) {
        var indexToGrab = orderLeft[j], val = eVecs[i][indexToGrab]

        newEigVecs[i][j] = val
      } else if(i == newEigVecs.length - 1) {
        var indexToGrab = orderRight[j], val = eVecs[i][indexToGrab]

        newEigVecs[i][j] = val
      } else {
        newEigVecs[i][j] = eVecs[i][j]
      }
    }
  }

  console.log(newEigVecs)
  
  return newEigVecs
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
    this.eigenOrderLeft = [0,1,2,3]
    this.eigenOrderRight = [0,1,2,3] 
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
    }, res

    try {
      res = await request("structure/modes", data, "POST");
    } catch(e) {
      console.error(e)
      console.log("Failed to get modes!")
      return
    }

    console.log(res.maxwell_matrices)

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
    }, res

    if(this.eigenvalues && this.eigenvectors && this.maxwell_matrices) {
      data = {
        ...data,
        maxwell_matrices: convertMatrixArrayToPythonParsable(this.maxwell_matrices),
        eigenvalues: convertEigenvaluesToPythonParsable(this.eigenvalues),
        eigenvectors: convertMatrixArrayToPythonParsable(this.eigenvectors)
      }
    }

    try {
      res = await request("structure/constants", data, "POST")
    } catch(e) {
      console.error(e)
      console.log("Failed to get constants!")
      return
    }

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

  permuteOrder(oldIndex, newIndex, toRight) {
    // bruh: https://stackoverflow.com/questions/5306680/move-an-array-element-from-one-array-position-to-another
    const array_move = (arr, old_index, new_index) => {
      if (new_index >= arr.length) {
          var k = new_index - arr.length + 1;
          while (k--) {
              arr.push(undefined);
          }
      }

      arr.splice(new_index, 0, arr.splice(old_index, 1)[0]);
      return arr; // for testing
    };

    if(toRight)
      this.eigenOrderRight = array_move(this.eigenOrderRight, oldIndex, newIndex)
    else
      this.eigenOrderLeft = array_move(this.eigenOrderLeft, oldIndex, newIndex)
  }

  resetPermuteOrder(toRight) {
    if(toRight)
      this.eigenOrderRight = [0,1,2,3]
    else
      this.eigenOrderLeft = [0,1,2,3]
  }

  getPermuteOrder(toRight) {
    if(toRight)
      return this.eigenOrderRight
    else
      return this.eigenOrderLeft
  }

  getPositionInPermuteOrder(id, toRight) {
    var arr = toRight ? this.eigenOrderRight : this.eigenOrderLeft

    for(let i = 0; i < arr.length; i++) {
      if(arr[i] === id) return i
    }
  }

  reorderPartOfOrder(startingPos, toRight) {
    var arr = toRight ? this.eigenOrderRight : this.eigenOrderLeft

    var preArr = arr.slice(0, startingPos), sortedArr = arr.slice(startingPos).sort(), newArr = []

    newArr.push.apply(newArr, preArr.concat(sortedArr))

    if(toRight)
      this.eigenOrderRight = newArr
    else
      this.eigenOrderLeft = newArr
  }

  placePositionBackInOrder(id, startingPos, toRight) {
    var arr = toRight ? this.eigenOrderRight : this.eigenOrderLeft

    var arrayContainsElem = false

    for(let i = startingPos; i < arr.length; i++) {
      if(arr[i] === id) arrayContainsElem = true
    }

    if(!arrayContainsElem) {
      var curPos = this.getPositionInPermuteOrder(id, toRight)

      this.permuteOrder(curPos, startingPos, toRight)
    }

    this.reorderPartOfOrder(startingPos, toRight)
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
    }, res

    if(this.eigenvalues && this.eigenvectors && this.maxwell_matrices) {
      data = {
        ...data,
        maxwell_matrices: convertMatrixArrayToPythonParsable(this.maxwell_matrices),
        eigenvalues: orderEigenvalues(convertEigenvaluesToPythonParsable(this.eigenvalues), this.eigenOrderLeft, this.eigenOrderRight),
        eigenvectors: orderEigenvectors(convertMatrixArrayToPythonParsable(this.eigenvectors), this.eigenOrderLeft, this.eigenOrderRight)
      }
    }

    console.log(data)

    try {
      res = await request("structure/field", data, "POST")
    } catch(e) {
      console.error(e)
      console.log("Failed to get field!")
      return
    }

    if(!this.eigenvalues && !this.eigenvectors && !this.maxwell_matrices) {
      console.log(res.maxwell_matrices)

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