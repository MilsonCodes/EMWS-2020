var EigenCalc = {};

/**
 * GetEigenvaluesAndEigenvectors
 * ----------------------
 * Takes a real matrix and returns the real Eigenvalues and Eigenvectors
 * using one of the many QR Algorithm methods.
 * 
 * http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
 * 
 * @param Matrix a 
 */
EigenCalc.getEigenvaluesAndEigenvectors = function(a) {
    var eigenvalues = new Array();
    var eigenvectors = new Array();

    console.time("eigTest")

    var characteristicPol = characteristicPolynomial(a)

    console.log( { A: a, characteristicPolynomial: characteristicPol })

    var roots = rootFinding(characteristicPol)

    console.log( { A: a, eigvals: roots })

    var eigVecs = findEigenvectors(a, roots)

    console.log( { A: a, eigenvalues: roots, eigenvectors: eigVecs } )

    var eigRes = organizePairOrdering(roots, eigVecs)

    console.timeEnd("eigTest")

    return eigRes
}

function printVectorOfComplexNums(vec) {
    var str = "[ ", arr = vec.toArray();
    
    for(var i = 0; i < arr.length; i++){
        str += "[" + math.complex(arr[i][0]).toString() + "]";

        if(i != arr.length - 1) str += ", "
    }

    str += " ]";

    console.log(str);
}

function isZeroBlock(matrix) {
    var hasNumber = false, matArr = matrix.toArray(), n = matrix.size()[0];

    for(let i = 0; i < n; i++) {
        for(let j = 0; j < n; j++) {
            if(matArr[i][j] != 0) hasNumber = true;
        }
    }

    return !hasNumber;
}

/**Mulitply Vector Matrices
 * ------------------------
 * Takes a regular vector and a transposed vector
 * and multiplies it into a square matrix.
 * 
 * @param Matrix v - Vector
 * @param Matrix vT - Transposed Vector
 */
function multiplyVectorMatrices(v, vT) {
    var resMat = math.zeros(v.size()[0], vT.size()[1]);

    var matArr = resMat.toArray();

    for(var i = 0; i < v.size()[0]; i++) {
        for(var j = 0; j < vT.size()[1]; j++){
            matArr[i][j] = math.multiply(math.subset(v, math.index(i, 0)), math.subset(vT, math.index(0, j)));
        }
    }

    resMat = math.matrix(matArr);

    //console.log(resMat);

    return resMat;
}

/** 
 * EIGENVALUE-FINDING SECTION
 */

/**
 * Characteristic Polynomial Algorithm
 * -----------------------------------
 * https://math.stackexchange.com/questions/405822/what-is-the-fastest-way-to-find-the-characteristic-polynomial-of-a-matrix
 * 
 * @param Matrix a
 */
function characteristicPolynomial(A) {
    var m = A.size()[0], n = A.size()[1], coeff = new Array();

    if (m != n) throw "Not a square matrix!";

    for(var i = 0; i <= n; i++) coeff[i] = 1

    var C = A
    for(var k = 1; k <= n; k++) {
        if(k > 1) C = math.multiply(A, math.add(C, math.multiply(coeff[n-k+1], math.identity(n))))

        coeff[n-k] = math.multiply(math.unaryMinus(math.divide(1,k)), math.trace(C))
    }

    return coeff
}

/**
 * Root Finding using Durand-Kerner method
 * ---------------------------------------
 * https://en.wikipedia.org/wiki/Durand%E2%80%93Kerner_method
 * 
 * @param Array coefficients 
 */
function rootFinding(coeff) {
    var n = coeff.length - 1, precision = 1e-14, prevVals = new Array(), curVals = new Array(), allUnderPrecision = false;

    for(var i = 0; i < n; i++) {
        curVals[i] = math.complex(0.4, 0.9);
        prevVals[i] = 0;
    }

    //console.log({ curVals: curVals, prevVals: prevVals })

    var iterationCount = 0

    while(!allUnderPrecision) {
        allUnderPrecision = true

        //console.log({ iterationCount: iterationCount, curVals: curVals })

        for(var i = 0; i < n; i++) {
            var curVal = curVals[i]

            var denominator = math.complex(1, 0)

            for(var j = 0; j < n; j++){
                if (j != i){
                    denominator = math.multiply(denominator, math.subtract(curVal, (j > i ? prevVals[j] : curVals[j])))
                    //console.log({ j: j, denominator: denominator})
                }
            }

            var funcRes = polynomialFunction(coeff, curVal)

            prevVals[i] = curVal
            curVals[i] = math.subtract(curVal, math.divide(funcRes, denominator))

            allUnderPrecision = checkPrecision(curVals[i], prevVals[i], precision)
        }

        iterationCount++

        if (iterationCount > 1000) break
    }

    //console.log({ iterationCount: iterationCount })

    return curVals
}

function polynomialFunction(coeff, x) {
    var val = 0;

    for(var i = 0; i < coeff.length; i++) {
        val = math.add(val, math.multiply(coeff[i], math.pow(x, i)))
    }

    return val;
}

function testPolynomial(coeff, roots) {
    var res = new Array()

    for(var i = 0; i < roots.length; i++)
        res[i] = polynomialFunction(coeff, roots[i])

    return res
}

function checkPrecision(val1, val2, precision) {
    var finVal = math.subtract(val1, val2)

    return Math.abs(finVal.re) < precision && Math.abs(finVal.im) < precision;
}

/**
 * EIGENVECTOR-FINDING SECTION
 */

 /**
  * Find Eigenvectors
  * ------------------------
  * Solves for the eigenvectors of matrix A using the 
  * 
  * @param {*} Matrix A
  * @param {*} Array eigenvectors 
  */
 function findEigenvectors(A, eigenvalues) {
    var eigenvectors = new Array()

     for(let i = 0; i < eigenvalues.length; i++) {
         var eigVal = eigenvalues[i]

         //A-lambdaI
         var Aeig = math.subtract(A, math.multiply(eigVal, math.identity(A.size()[0], A.size()[1])))

        var vec = solveVectorNullspace(Aeig)        

        eigenvectors[i] = math.transpose(vec).toArray()[0]
     }

     return math.matrix(eigenvectors)
 }

 /**
  * Solve Vector Nullspace
  * -----------------------
  * Algorithm that solves for the vector that completes the nullspace (other than 0).
  * 
  * Ax = 0
  * 
  * @param {*} Matrix a
  */
 function solveVectorNullspace(A) {
    var Agauss = performGaussianElimination(A)

    //console.log({ A: A, gaussElim: Agauss })

    var vector = solveSystem(Agauss)

    //console.log({ A: A, vector: vector, gauss: Agauss, check: matrixMultiplication(A, vector) })

    return vector 
 }

 /**
  * PerformGaussianElimination
  * ---------------------------
  * Performs Gaussian Elimination on the (A-lambdaI) matrix
  * 
  * @param {Matrix} A 
  */
 function performGaussianElimination(A) {
    var m = A.size()[0], n = A.size()[1], h = 0, k = 0;

    while (h < m && k < n){
       var pivotLoc = h, pivotVal = 0;

        for(let i = h; i < m; i++) {
            var val = A.get([i, k]), comp = math.abs(math.add(math.pow(val.re, 2), math.pow(val.im, 2)))

            if(pivotVal < comp) {
               pivotVal = comp;
               pivotLoc = i;
            }
        }

        var checkVal = A.get([pivotLoc, k])

        if (isZero(checkVal, 1e-9))
            k++
        else {
            A = swapMatrixRows(A, h, pivotLoc)

            for(let i = h+1; i < m; i++) {
                var f = math.divide(A.get([i, k]), A.get([h,k]))

                A.set([i,k], 0)

                for(let j = k+1; j < n; j++)
                    A.set([i,j], math.subtract(A.get([i,j]), math.multiply(A.get([h,j]), f)))
            }

            h++
            k++
        }
    }


    return A;
 }

 /**
  * Solve System
  * ----------------
  * Solves the system of the matrix that has gone through Gaussian elimination
  * @param {Matrix} A 
  */
 function solveSystem(A) {
    var m = A.size()[0], n = A.size()[1];

    var vector = math.zeros(m, 1)

    for(let i = m-1; i >= 0; i--) {
        var isRowZero = true

        for(let j = 0; j < n; j++)
            if(!isZero(A.get([i,j]), 1e-6)) isRowZero = false

        if(isRowZero)
            vector.set([i, 0], math.complex(1))
        else if(i == m-1 && !isRowZero) {
            console.log(A._data[i])
            throw "Last row is not a zero row!"
        } else {
            vector.set([i, 0], solveEquation(A._data[i], i, vector))
        }
    }

    return math.divide(vector, math.sqrt(complexVectorDotProduct(vector, vector)))
 }

 function swapMatrixRows(A, row1, row2) {
    var temp = A._data[row2]
    A._data[row2] = A._data[row1]
    A._data[row1] = temp

    return A
 }

 function isZero(val, precision) {
     return ((val.re < precision && val.re > -precision) && (val.im < precision && val.im > -precision)) || val == 0
 }

 function solveEquation(coeff, curRow, vector) {
    var val = 0

    for(let i = curRow+1; i < coeff.length; i++) 
        val = math.add(val, math.multiply(math.unaryMinus(coeff[i]), vector.get([i, 0])))

    return math.divide(val, coeff[curRow])
 }

/**
 * ORDERING VALUES
 */

 function organizePairOrdering(eigenvalues, eigenvectors) {
    const PRECISION = 1e-7

    for(var i = 0; i < eigenvalues.length; i++) {
        var eigVal = eigenvalues[i]

        if (!isNumZero(eigVal.im, PRECISION) && isNumZero(eigVal.re, PRECISION)) {
            //console.log(`Eigenvalue ${eigVal} is imaginary`)

            if(eigVal.im > 0 && i != 0) {
                //console.log(`Eigenvalue ${eigVal} is positive imaginary. Swapping to front!`)
                eigenvalues = swapArrayIndexes(eigenvalues, i, 0)
                eigenvectors = swapMatrixRows(eigenvectors, i, 0)
            }
        } else if(isNumZero(eigVal.im, PRECISION) && !isNumZero(eigVal.re, PRECISION)) {
            //console.log(`Eigenvalue ${eigVal} is real`)
            
            if(eigVal.re < 0 && i != 1) {
                //console.log(`Eigenvalue ${eigVal} is negative real. Swapping to second index!`)
                eigenvalues = swapArrayIndexes(eigenvalues, i, 1)
                eigenvectors = swapMatrixRows(eigenvectors, i, 1)
            }
        }
    }

    return { eigenvalues: eigenvalues, eigenvectors: eigenvectors }
 }

 function isNumZero(num, precision) {
    return (num < precision && num > -precision) || num == 0;
 }

 function swapArrayIndexes(array, val1, val2) {
    var temp = array[val2]
    array[val2] = array[val1]
    array[val1] = temp

    return array
 }

/**
 * Complex Vector Dot Product
 * -------------------------
 * Computes the dot product of two complex vectors
 * 
 * @param {Vector} a 
 * @param {Vector} b 
 */
function complexVectorDotProduct(a, b) {
    var m = a.size()[0], n = a.size()[1], p = b.size()[0], q = b.size()[1], sum = 0.0;

    if (m != p || n != q) throw "The provided vectors are not the same size!"

    if(m == 1 && n >= 1) {
        for(var i = 0; i < n; i++) {
            var num = math.multiply(math.conj(a.get([0,i])), b.get([0,i]))

            sum += num.re
        }
    } else if(n == 1 && m >= 1) {
        for(var i = 0; i < m; i++) {
            var num = math.multiply(math.conj(a.get([i,0])), b.get([i,0]))

            sum += num.re
        }
    }

    return sum
}