// Global namespace
var emScattering2 = emScattering2 || {};


// Intrinsic impedance of free space in Ohms
emScattering2.ETTA_0 = 376.73031;

// Speed of light in free space in m/s -- CURRENTLY IN UNITS WHERE THE NUMERICAL VALUE IS 1
// emScattering2.C = 299792458.0; // m/s
emScattering2.C = 1;

/*
Layer Object
------------------------------------------------------------------------------------
*/

/**
Represents a single material layer in a 1D photonic crystal. 'epsilon' and 'mu' 
represent the relative electric permittivity and the relative magnetic permeability, 
respectively, of the material. Currently, only isotropic materals are suported, so 'epsilon' 
and 'mu' are scalars. 'length' represents the length of the layer.
*/
emScattering2.Layer = function(eMat, mMat, length) {
    this.epsilon = eMat;
    this.mu = mMat;
    this.length = length;
};

/* 

GENERAL FUNCTIONS 
------------------------------------------------------------------------------------
*/

/**
 Expects array of epsilon values, array of mu values, and the scalar number of layers,
 Parses user input from web into a usable standardized format
 Outputs Array containing arrays of real and imaginary parts for epsilon and mu
 */
emScattering2.Parse = function(eArray,mArray,numLayers){
    var i,j,k, tmpe, tmpm, eParse = [], mParse = [];            
    for( i = 0; i < numLayers; i++){
        tmpm = mArray[i];
        tmpe = eArray[i];
        for(j = 0; j < 3; j++){
            for(k = 0; k < 3; k++){
               
               tmpm[j][k] = math.complex(tmpm[j][k]);
               tmpe[j][k] = math.complex(tmpe[j][k]);
               
            }
        }
        mParse[i] = math.matrix([[tmpm[0][0],tmpm[0][1],tmpm[0][2]], [tmpm[1][0],tmpm[1][1],tmpm[1][2]], [tmpm[2][0],tmpm[2][1],tmpm[2][2]]]);
        eParse[i] = math.matrix([[tmpe[0][0],tmpe[0][1],tmpe[0][2]], [tmpe[1][0],tmpe[1][1],tmpe[1][2]], [tmpe[2][0],tmpe[2][1],tmpe[2][2]]]);
            
    }
    
    return [eParse,mParse];
    
};


/**
Generates an array of layer objects. Expects array of Epsilon Matrices,
Mu Matrices, and lengths.
*/
emScattering2.createLayers = function(eMat,mMat, lengths) {
    var layers = [];
    for (var i = 0; i < eMat.length; i++) {
        var layer = new emScattering2.Layer(eMat[i], mMat[i], lengths[i]);
        layers.push(layer);
    }
    return layers;
};

/**
 Takes an epislon matrix, mu matrix, scalar kx, and scalar ky
 Generates the anisotropic maxwell matrix for the given values.
 This matrix is found on page 9, slide 17 from the following link
    http://emlab.utep.edu/ee5390cem/Lecture%204%20--%20Transfer%20Matrix%20Method.pdf

 The actual matrix in the Maxwell ODE system is the output of this function times omega.
 */
emScattering2.Maxwell = function(eMat, mMat, kx, ky){
    var i = math.complex('i'), exx = eMat[0][0], exy = eMat[0][1], exz = eMat[0][2], eyx = eMat[1][0], eyy = eMat[1][1],
        eyz = eMat[1][2], ezx = eMat[2][0], ezy = eMat[2][1], ezz = eMat[2][2],mxx = mMat[0][0], mxy = mMat[0][1], 
        mxz = mMat[0][2], myx = mMat[1][0], myy = mMat[1][1], myz = mMat[1][2], mzx = mMat[2][0], mzy = mMat[2][1], mzz = mMat[2][2],
        m = math.matrix();
                                                                        
    m.set([0,0],math.add(math.multiply(ky,math.divide(myz,mzz)),math.multiply(math.divide(ezx,ezz),kx)));
    m.set([0,1],math.multiply(kx,math.subtract(math.divide(myz,mzz),math.divide(ezy,ezz))));
    m.set([0,2],math.subtract(math.add(math.divide(math.multiply(kx,ky),ezz),myx),math.divide(math.multiply(myz,mzx),mzz)));
    m.set([0,3],math.subtract(math.add(math.unaryMinus(math.divide(math.pow(kx,2),ezz)),myy),math.divide(math.multiply(myz,mzy),mzz)));
    m.set([1,0],math.multiply(ky,math.subtract(math.divide(mxz,mzz),math.divide(ezx,ezz))));
    m.set([1,1],math.add(math.multiply(kx,math.divide(mxz,mzz)),math.multiply(math.divide(ezy,ezz),ky)));
    m.set([1,2],math.add(math.subtract(math.divide(math.pow(ky,2),ezz),mxx),math.divide(math.multiply(mxz,mzx),mzz)));
    m.set([1,3],math.add(math.subtract(math.unaryMinus(math.divide(math.multiply(kx,ky),ezz)),mxy),math.divide(math.multiply(mxz,mzy),mzz)));

    //m.set([2,0],math.subtract(math.add(math.divide(math.multiply(kx,ky),mzz),eyx),math.divide(math.multiply(eyz,ezx),ezz)));              //Old
    //m.set([2,1],math.subtract(math.add(math.unaryMinus(math.divide(math.pow(kx,2),mzz)),eyy),math.divide(math.multiply(eyz,ezy),ezz)));
    m.set([2,0],math.add(math.subtract(math.unaryMinus(math.divide(math.multiply(kx,ky),mzz)),eyx),math.divide(math.multiply(eyz,ezx),ezz)));//New
    m.set([2,1],math.add(math.subtract(math.divide(math.pow(kx,2),mzz),eyy),math.divide(math.multiply(eyz,ezy),ezz)));
    m.set([2,2],math.add(math.multiply(ky,math.divide(eyz,ezz)),math.multiply(math.divide(mzx,mzz),kx)));
    m.set([2,3],math.multiply(kx,math.subtract(math.divide(eyz,ezz),math.divide(mzy,mzz))));
    //m.set([3,0],math.add(math.subtract(math.divide(math.pow(ky,2),mzz),exx),math.divide(math.multiply(exz,ezx),ezz)));                      //Old
    //m.set([3,1],math.add(math.subtract(math.unaryMinus(math.divide(math.multiply(kx,ky),mzz)),exy),math.divide(math.multiply(exz,ezy),ezz)));
    m.set([3,0],math.subtract(math.add(math.unaryMinus(math.divide(math.pow(ky,2),mzz)),exx),math.divide(math.multiply(exz,ezx),ezz)));   //New
    m.set([3,1],math.subtract(math.add(math.divide(math.multiply(kx,ky),mzz),exy),math.divide(math.multiply(exz,ezy),ezz)));
    m.set([3,2],math.multiply(ky,math.subtract(math.divide(exz,ezz),math.divide(mzx,mzz))));
    m.set([3,3],math.add(math.multiply(kx,math.divide(exz,ezz)),math.multiply(math.divide(mzy,mzz),ky)));
    
    
    return math.multiply(m, i);
};

/**
 Calculates the four eigenvalues for the 4x4 matrix of a single layer using the
 C++ wrapped function complex_eigenvalues. Expects a maxwell matrix,
 and the wrapped complex_eigenvalues function. 
 Outputs array containing the complex eigenvalues.
*/
emScattering2.calcEigsVa = function(maxwell, complex_eigenvalues){
    var i,j, buf, tmp, foo = [], arr = [], ret = [];
    
    for( i = 0; i < maxwell._size[0]; i++){
        for( j = 0; j < maxwell._size[1]; j++){
            foo.push(maxwell._data[i][j].re);
            foo.push(maxwell._data[i][j].im);
        }
    }
   
    buf = Module._malloc(foo.length*8);
    Module.HEAPF64.set(new Float64Array(foo),buf/8);
    
    for (i = 0; i < 8; i++){
        tmp = Module.getValue(complex_eigenvalues(buf) + (8*i),'double');
        arr.push(tmp);
    }

    
    for (i = 0, j = 0; i < 8; i = i+2, j++){
        foo = math.complex(arr[i],arr[i+1]);
        ret[j] = foo;
    }
    
    Module._free(buf);
    return ret;
};

/**
 Calculates the eigenvectors of the complex 4x4 maxwell matrix using the C++ function 
 complex_eigenvectors.  Expects a maxwell matrix and the 
 wrapped complex_eigenvectors function. Outputs matrix of eigenvectors.
 */
emScattering2.calcEigsVe = function (maxwell,complex_eigenvectors){
    var i,j, buf, tmp, ret, arr = [], foo = [];
    ret = math.matrix();
    for( i = 0; i < maxwell._size[0]; i++){
        for( j = 0; j < maxwell._size[1]; j++){
            foo.push(maxwell._data[i][j].re);
            foo.push(maxwell._data[i][j].im);
        }
    }

    buf = Module._malloc(foo.length*8);
    Module.HEAPF64.set(new Float64Array(foo),buf/8);
    
    for (i = 0; i < 32; i++){
        tmp = Module.getValue(complex_eigenvectors(buf) + (8*i),'double');
        arr.push(tmp);
    }

    tmp = 0;
    for( i = 0; i < 4; i++){
        for( j = 0; j < 4; j++, tmp = tmp + 2){
            ret.set([i,j], math.complex(arr[tmp],arr[tmp+1]));
        }
    }

    Module._free(buf);
    return ret;
};

/**
 Expects a diagonal matrix of eigenvalues, and a znorm
 Returns a diagonal matrix where each element on the diagonal is of the form: exp[eigenvalue * znorm]
 */
emScattering2.expEigenvaluesDiag = function(eigenvalues,znorm){
   var tmp = [eigenvalues[0],eigenvalues[1],eigenvalues[2],eigenvalues[3]];
   return math.diag(math.exp(math.multiply(tmp,znorm)));
};

/**Calculates the Constants for the Scattering problem*/
emScattering2.calculateConstants = function(S,Modes,T0) {
    var a = math.zeros(S._data.length),tmp,b0, bN, unknown_c, c = math.zeros(S._data.length + 3);

    c._data[0] = Modes[0];
    c._data[1] = Modes[1];
    c._data[c._data.length - 1] = Modes[2];
    c._data[c._data.length] = Modes[3];
    
    tmp = math.matrix([
                       [T0._data[0][0],T0._data[0][1]],
                       [T0._data[1][0],T0._data[1][1]],
                       [T0._data[2][0],T0._data[2][1]],
                       [T0._data[3][0],T0._data[3][1]]
                      ]);
    
    b0 = math.multiply(tmp,math.matrix([[Modes[0]],[Modes[1]]]));
    for(var i = 0; i < b0._data.length; i++ )
        a._data[i] = b0._data[i];
    
    tmp = math.matrix([
                       [0,0],
                       [0,0],
                       [-1,0],
                       [0,-1]
                     ]);
    
    bN = math.multiply(tmp,math.matrix([[Modes[2]],[Modes[3]]]));
    for(var i = 0; i < bN._data.length; i++)
        a._data[a._data.length - bN._data.length + i] = bN._data[i];

    unknown_c = math.lusolve(S,math.unaryMinus(a._data));


    
    for(var i = 0; i < unknown_c._data.length; i++)
        c._data[i + 2] = unknown_c._data[i][0];

    return c;
};



//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================


/**
Eigenpair Object and Methods
------------------------------------------------------------------------------------
An Eigenpair is a pairing between an eigenvalue and eigenvector that has either rightward or leftward propagation
*/
emScattering2.Eigenpair = function(eigenvalue, eigenvector) {
    this.eigenvalue = eigenvalue;
    this.eigenvector = eigenvector;
    this.isRightward = false;      //if 1, rightward, if 0, leftward
};

emScattering2.Eigenpair.prototype.setOrientation = function(orientation){
    if(orientation === true || orientation === false)
        this.isRightward = orientation;
    else
        throw "Only true and false are allowed for orientation";
    
};


//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================


/**
Struct Object and Methods
------------------------------------------------------------------------------------
A Struct is a collection of layers and values corresponding to each layer. Methods included are used to solve the scattering problem.
*/
emScattering2.Struct = function(layers, lengths, kx, ky, omega, Modes) {
    this.layers = layers;
    this.lengths = lengths;
    this.numLayers = layers.length;
    this.maxwellMatrices = new Array(this.numLayers);
    this.Eigensystems = new Array(this.numLayers);
    this.eigenvectors = new Array(this.numLayers);
    this.eigenvalues = new Array(this.numLayers); 
    this.transferMatrices = new Array(this.numLayers - 1);
    this.kx = kx;
    this.ky = ky;
    this.omega = omega;
    this.Modes = Modes;
    this.scatteringMatrix;
};


/**
 Constructor for structure, expects matrix of epsilon matrices, matrix of mu matrices, kx, ky and omega scalars
 returns object of type Struct
 */
emScattering2.makeStructure = function(eMat,mMat,length, kx, ky ,omega,Modes) {
    var Structure = new emScattering2.Struct(emScattering2.createLayers(eMat, mMat, length),length,kx ,ky ,omega, Modes);
    return Structure;   

};

/**
 Creates the anisotropic maxwell matrix for each layer and stores in Struct object.
 The output of emScattering2.Maxwell needs to be multiplied by this.omega to obtain the matrix for the Maxwell ODE system.
 */
emScattering2.Struct.prototype.makeMaxwell = function(){
    var tmp, kx_ = this.kx/this.omega, ky_ = this.ky/this.omega;
    for( var i = 0; i < this.numLayers; i++){
        this.maxwellMatrices[i] = math.multiply(emScattering2.Maxwell(this.layers[i].epsilon._data, this.layers[i].mu._data, kx_, ky_), this.omega);
    }
    
};

/**
Calculates Eigenvalues and Eigenvectors for each layer and returns them in an array.
*/
emScattering2.Struct.prototype.calcEigs = function(){
    var i, maxwell, complex_eigenvalues, complex_eigenvectors;
    var ret_values = new Array(this.numLayers), ret_vectors = new Array(this.numLayers);
    complex_eigenvalues = Module.cwrap('complex_eigenvalues','number',['number']);
    complex_eigenvectors = Module.cwrap('complex_eigenvectors','number',['number']);
    for( i = 0; i<this.numLayers; i++){
        maxwell = this.maxwellMatrices[i];
        //ret_values[i] = emScattering2.calcEigsVa(maxwell, complex_eigenvalues);
        //ret_vectors[i] = emScattering2.calcEigsVe(maxwell, complex_eigenvectors);

        var eigResults = EigenCalc.getEigenvaluesAndEigenvectors(maxwell), omega = this.omega;

        console.log({omega, maxwell, eigResults});

        ret_values[i] = eigResults.eigenvalues;
        ret_vectors[i] = eigResults.eigenvectors;

        //console.log(EigenCalc.getEigenvaluesAndEigenvectors(maxwell));                     //Calculates eigenvectors using returned eigenvalues - currently printed to console to check for correction
    }
    return [ret_values, ret_vectors];
};

/**
 For each layer sorts eigenvalues according to the parity of their real and imaginary parts into rightward and leftward pointing modes
 If re(eigenvlaue) > 0 the eigenvalue is considered leftward, if re(eigenvalue) < 0 then it is considered rightward
 If |re(eigenvalue)| < 10^-14 then we consider the imaginary part of the eigenvalue. If im(eigenvalue) > 0 then right, otherwise left.
 Each layer has exactly 2 right and 2 left. Stores in the structure data struct
 */
emScattering2.Struct.prototype.organizeEigenvalues = function(){
    var rightward, leftward;
    for(var i = 0; i < this.numLayers; i++){
        rightward = 0;
        leftward = 0;
        for (var j = 0; j < 4; j++) {
            if(rightward === 2){
                ++leftward;
                this.Eigensystems[i][j].setOrientation(false);
            }
            else if(leftward === 2){
                ++rightward;
                this.Eigensystems[i][j].setOrientation(true);
            }
            else {
                if(math.abs(math.re(this.Eigensystems[i][j].eigenvalue)) < math.pow(10,-14) ){
                   if(math.im(this.Eigensystems[i][j].eigenvalue) > 0){
                       ++rightward;
                       this.Eigensystems[i][j].setOrientation(true);
                   }
                   else{
                       ++leftward;
                       this.Eigensystems[i][j].setOrientation(false);
                   }
                }
                else if(math.re(this.Eigensystems[i][j].eigenvalue) > 0){
                    ++leftward;
                    this.Eigensystems[i][j].setOrientation(false);
                }
                else if(math.re(this.Eigensystems[i][j].eigenvalue) < 0){
                    ++rightward;
                    this.Eigensystems[i][j].setOrientation(true);
                }
            }
        }
    }
};

/**Calculates the eigensystem for each layer where an Eigensystem is a collection of Eigenpairs 
* Stores the completed system in the Struct object.
*/
emScattering2.Struct.prototype.calcEigensystems = function(){
    var tmp, esystem = new Array(this.numLayers), pairs = new Array(4), eva, eve;
    tmp = this.calcEigs(math);
    eva = tmp[0], eve = tmp[1];
    for(var j = 0; j < this.numLayers; j++){
        for(var i = 0; i < 4; i++){
            pairs[i] = new emScattering2.Eigenpair(eva[j][i],eve[j]._data[i]);
        }
        esystem[j] = pairs;
    }
    this.eigenvectors = eve;
    this.eigenvalues = eva;
    this.Eigensystems = esystem;
    this.organizeEigenvalues();
};

/**
  Calculates the transfer matrices for the interface between each layer
  Stores in Struct object
 */

emScattering2.Struct.prototype.calcTransfer = function(){
    //Matrices per interface
    if (this.numLayers > 1){
        for(var i = 0, N = this.transferMatrices.length; i < N; i ++){
            var w, wNext, tmp, zNorm, expDiagonal;
            wNext = math.inv(this.eigenvectors[i+1]);
            w = this.eigenvectors[i];
            if (i === 0)
                zNorm = 0;
            else
                zNorm = math.multiply (this.omega,this.layers[i].length);
            expDiagonal = emScattering2.expEigenvaluesDiag(this.eigenvalues[i], zNorm);
            tmp = math.multiply(w,expDiagonal);
            this.transferMatrices[i] = math.multiply(wNext,tmp);
        }
    }

};

/**
Constructs the matrix used to solve the scattering problem (i.e. finding the unknown constants governing propagation
given the incoming constants and the calculated transfer matrices) and stored in Struct 
 */

emScattering2.Struct.prototype.calcScattering = function(){
    var tm = this.transferMatrices, S = math.zeros(tm.length*4,tm.length*4);

    //Set -1's
    for(var i = 0; i < tm.length*4 - 2; i++){
        S._data[i][i+2] = -1;
    }
    //Set top 2x2 transfer matrix block
    for( var i = 0; i < 4;i++){
        for( var j = 0; j < 2; j++){
        S._data[i][j] = tm[0]._data[i][2 + j];
        }
    }
    //Set transfer matrices
    for(var k = 1; k < tm.length; k++){
        for( var i = 0; i < 4;i++){
            for( var j = 0; j < 4; j++){
                if( k === 1)
                    S._data[4*k+i]  [2*k+j] = tm[k]._data[i][j];
                else
                    S._data[4*k+i][4*k+j-2] = tm[k]._data[i][j];
            }
        }
    }
    
   
    this.scatteringMatrix = S;

};








//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================



/**
 Photonic Crystal Object and Methods
------------------------------------------------------------------------------------
A Photonic Crystal is a wrapper object for a structure. Methods included are used to complete and display the experiment.
*/
emScattering2.PhotonicCrystal = function(Struct){
    this.Struct = Struct;
};

emScattering2.createPhotonicCrystal = function(Struct){
    var crystal = new emScattering2.PhotonicCrystal(Struct);
    return crystal;
};

/**
Returns the field values in all the layers given the coefficients of the incoming modes.
Returned object has properties z for the coordinates used, Ex, Ey, Hx, 
and Hy. There is a one-to-one correspondence between an element in z and the other 4 Arrays.
*/
emScattering2.PhotonicCrystal.prototype.determineField = function() {
    var numLayers = this.Struct.numLayers, numPoints = 100, layerNormZ, W, lambda, c, current_c,
    expDiag, result, currentLeftZ = 0, currentRightZ = this.Struct.layers[0].length,
    _Ex = new Array(), _Ey = new Array(), _Hx = new Array(), _Hy = new Array(), _z = new Array();
    
    c = emScattering2.calculateConstants(this.Struct.scatteringMatrix,this.Struct.Modes,this.Struct.transferMatrices[0]);           //Creates a constant vector using the scattering matrix, coefficients, and transfer matrix

    console.log(c);

    current_c = c._data.slice(0,4);                                                                         //Slices the constant vector to a smaller vector
    for(var i = 0; i < numLayers; i++){
        if(i === 0){
            layerNormZ = numeric.linspace(-this.Struct.omega*this.Struct.layers[i].length,0, 
                                 Math.floor(this.Struct.layers[i].length)*numPoints);                       //Creates an array of Z values (with omega taken into account) with the size of the length and number of points per Z value (Ex. 10 size of layer * 100 points)
        }
        else{
            layerNormZ = numeric.linspace(0, this.Struct.omega*this.Struct.layers[i].length,
                                 Math.floor(this.Struct.layers[i].length)*numPoints);
        }
        
        W = this.Struct.eigenvectors[i];                                                                //Sets W to the current layer eigenvector matrix
        lambda = this.Struct.eigenvalues[i];                                                            //Sets lambda to the current layer eigenvalues
        for(var j = 0; j < layerNormZ.length; j++){
            expDiag = emScattering2.expEigenvaluesDiag(lambda, layerNormZ[j]);                          //Creates a diagonal matrix with the exponential of each eigenvalue times the current z value
            result = math.multiply(W,(math.multiply(expDiag,current_c)));                               //Multiplies the exponential diagonal times the eigenvector matrix times the constant vector
            _Ex.push(result._data[0].re);
            _Ey.push(result._data[1].re);
            _Hx.push(result._data[2].re);
            _Hy.push(result._data[3].re);
        }    
         _z = _z.concat(numeric.linspace(currentLeftZ, currentRightZ, 
                        numPoints*Math.floor(this.Struct.layers[i].length)));                           //Adds Z values to an array (Z values are DIFFERENT from calculated above)
         if (i+1 < numLayers) {
            currentLeftZ += this.Struct.layers[i].length;                                               //Extends lower layer limit
            currentRightZ += this.Struct.layers[i+1].length;                                            //Extends upper layer limit
            current_c = c._data.slice(4+4*i,8+4*i);                                                     //Shifts constant vector
        } 
    }
    
    return {z: _z, Ex: _Ex, Ey: _Ey, Hx: _Hx, Hy: _Hy};
};

/**DetermineFieldAtZPoint
 * ----------
 * Takes a point of z and calculates the value of the fields (Ex, Ey, Hx, and Hy) at
 * that point. Extreme WIP - This code will need to be rewritten and, once working,
 * will allow for much more efficient Transmission calculations
 */
emScattering2.PhotonicCrystal.prototype.determineFieldAtZPoint = function(zPoint){
    var W, lambda, c, current_c, expDiag, result, currentLayer, zIntoLayer = zPoint, normZ, interfaces = this.materialInterfaces();

    for(var i = 0; i < this.Struct.numLayers; i++){
        if(zPoint == 0 && i == 0) currentLayer = i;

        if(interfaces[i] < zPoint && interfaces[i+1] >= zPoint) currentLayer = i;
    }

    //console.log("Layer: " + currentLayer);

    for(var i = 0; i < currentLayer; i++){
        zIntoLayer -= this.Struct.layers[i].length;
    }

    //console.log("Z-Point in Layer: " + zIntoLayer);

    if (this.Struct.omega < 0){
        if(currentLayer === 0) normZ = (-this.Struct.omega * this.Struct.layers[currentLayer].length) - math.abs(zIntoLayer * this.Struct.omega);
        else normZ = (this.Struct.omega * this.Struct.layers[currentLayer].length) + math.abs((zIntoLayer - this.Struct.layers[currentLayer].length) * this.Struct.omega);
    } else if(this.Struct.omega > 0) {
        if(currentLayer === 0) normZ = (-this.Struct.omega * this.Struct.layers[currentLayer].length) + math.abs(zIntoLayer * this.Struct.omega);
        else normZ = (this.Struct.omega * this.Struct.layers[currentLayer].length) - math.abs((zIntoLayer - this.Struct.layers[currentLayer].length) * this.Struct.omega);
    }

    //console.log("Value used to calculate: " + normZ);
    
    c = emScattering2.calculateConstants(this.Struct.scatteringMatrix,this.Struct.Modes,this.Struct.transferMatrices[0]);
    
    current_c = c._data.slice(0,4);

    for(var i = 0; (i < currentLayer) && (i < this.Struct.numLayers - 1); i++){
        current_c = c._data.slice(4+4*i,8+4*i);
    }

    W = this.Struct.eigenvectors[currentLayer];
    lambda = this.Struct.eigenvalues[currentLayer];
    
    expDiag = emScattering2.expEigenvaluesDiag(lambda, normZ);
    result =  math.multiply(W,(math.multiply(expDiag,current_c)));

    return {Ex: result._data[0].re, Ey: result._data[1].re, Hx: result._data[2].re, Hy: result._data[3].re};
}

/**
Setup for Electric Field for Mathbox that is run each time CreateAnim() is run with the runsetup flag set
Outputs Ex and Ey in complex polar form, where each part is a separate array output 
 */
emScattering2.PhotonicCrystal.prototype.mathboxSetupEf = function() {
    var numLayers = this.Struct.numLayers,W, lambda, c, current_c,
    tmpRX = [],tmpPhiX = [],Ex = [],Ey = [],ExR = [],EyR = [],ExPhi = [],EyPhi = [],
    tmpRY = [],tmpPhiY = [];
    
    
    c = emScattering2.calculateConstants(
            this.Struct.scatteringMatrix,this.Struct.Modes,this.Struct.transferMatrices[0]);
    current_c = c._data.slice(0,4);
    
    for(var i = 0; i < numLayers; i++){
        lambda = this.Struct.eigenvalues[i];
        W = this.Struct.eigenvectors[i];
        Ex = math.dotMultiply(current_c,W._data[0]);
        Ey = math.dotMultiply(current_c,W._data[1]);
        for(var k = 0; k < Ex.length; k++) {
            Ex[k] = Ex[k].toPolar();
            Ey[k] = Ey[k].toPolar();
        }

        for(var k = 0; k < Ex.length; k++){
            tmpRX[k] = Ex[k].r;
            tmpPhiX[k]= Ex[k].phi;
            tmpRY[k] = Ey[k].r;
            tmpPhiY[k] = Ey[k].phi;            
        }
   
        ExR.push(tmpRX);
        ExPhi.push(tmpPhiX);
        EyR.push(tmpRY);
        EyPhi.push(tmpPhiY);
        
        
        tmpRX = [],tmpPhiX = [];
        tmpRY = [],tmpPhiY = [];
        Ex = [],Ey = [];
        
        
        if (i+1 < numLayers) {
            current_c = c._data.slice(4+4*i,8+4*i);
            
        } 
    }
    
    return {ExR: ExR, ExPhi: ExPhi, EyR: EyR, EyPhi: EyPhi};
};

/**
Calculates the Ex and Ey bits for Mathbox graphing.
Expects the array of layer lengths, the current timestamp t, the position z, 
ExR, ExPhi, EyR, EyPhi.
 */
emScattering2.PhotonicCrystal.prototype.mathboxEf = function(lengths,t,z,ExR,ExPhi,EyR,EyPhi) {
    var layerNum = 0, inf = 0, sup = lengths[0], lambdas = this.Struct.eigenvalues, vecs = this.Struct.eigenvectors, omega = this.Struct.omega,
            Ex,Ex1,Ex2,Ex3,Ex4,Ey,Ey1,Ey2,Ey3,Ey4;
    for(var i = 0; i < lengths.length; i++){
        if(inf <= z && z <= sup){
            layerNum = i;
            break;
        }
        else{
            inf = sup;
            sup = sup + lengths[i + 1];
        }          
    }
    
    // Ex1 = ExR[layerNum][0] * math.exp(math.re(lambdas[layerNum][0]) * z) * math.cos(ExPhi[layerNum][0] + math.im(lambdas[layerNum][0]) * z - omega * t);
    // Ex2 = ExR[layerNum][1] * math.exp(math.re(lambdas[layerNum][1]) * z) * math.cos(ExPhi[layerNum][1] + math.im(lambdas[layerNum][1]) * z - omega * t);
    // Ex3 = ExR[layerNum][2] * math.exp(math.re(lambdas[layerNum][2]) * z) * math.cos(ExPhi[layerNum][2] + math.im(lambdas[layerNum][2]) * z - omega * t);
    // Ex4 = ExR[layerNum][3] * math.exp(math.re(lambdas[layerNum][3]) * z) * math.cos(ExPhi[layerNum][3] + math.im(lambdas[layerNum][3]) * z - omega * t);
    // Ex = Ex1 + Ex2 + Ex3 + Ex4;
    
    // Ey1 = EyR[layerNum][0] * math.exp(math.re(lambdas[layerNum][0]) * z) * math.cos(EyPhi[layerNum][0] + math.im(lambdas[layerNum][0]) * z - omega * t);
    // Ey2 = EyR[layerNum][1] * math.exp(math.re(lambdas[layerNum][1]) * z) * math.cos(EyPhi[layerNum][1] + math.im(lambdas[layerNum][1]) * z - omega * t);
    // Ey3 = EyR[layerNum][2] * math.exp(math.re(lambdas[layerNum][2]) * z) * math.cos(EyPhi[layerNum][2] + math.im(lambdas[layerNum][2]) * z - omega * t);
    // Ey4 = EyR[layerNum][3] * math.exp(math.re(lambdas[layerNum][3]) * z) * math.cos(EyPhi[layerNum][3] + math.im(lambdas[layerNum][3]) * z - omega * t);
    // Ey = Ey1 + Ey2 + Ey3 + Ey4;
    // var tmpX = Ex;
    // var tmpY = Ey;
    Ex = 0;
    Ey = 0;
    for (i = 0; i < 4; i++){
        Ex += ExR[layerNum][i] * math.exp(math.re(lambdas[layerNum][i]) * z) * math.cos(ExPhi[layerNum][i] + math.im(lambdas[layerNum][i]) * z - omega * t);
        Ey += EyR[layerNum][i] * math.exp(math.re(lambdas[layerNum][i]) * z) * math.cos(EyPhi[layerNum][i] + math.im(lambdas[layerNum][i]) * z - omega * t);
    }
    // console.log(Ex);
    // console.log(Ey);
    return {Ex: Ex, Ey: Ey};
};

/**
Setup for Magnetic Field for Mathbox that is run each time CreateAnim() is run with the runsetup flag set
Outputs Hx and Hy in complex polar form, where each part is a separate array output 
 */
emScattering2.PhotonicCrystal.prototype.mathboxSetupHf = function() {
    var numLayers = this.Struct.numLayers,W, lambda, c, current_c,
    tmpRX = [],tmpPhiX = [],Hx = [],Hy = [],HxR = [],HyR = [],HxPhi = [],HyPhi = [],
    tmpRY = [],tmpPhiY = [];
    
    c = emScattering2.calculateConstants(this.Struct.scatteringMatrix,this.Struct.Modes,this.Struct.transferMatrices[0]);
    current_c = c._data.slice(0,4);
    
    for(var i = 0; i < numLayers; i++){
        lambda = this.Struct.eigenvalues[i];
        W = this.Struct.eigenvectors[i];
        Hx = math.dotMultiply(current_c,W._data[2]);
        Hy = math.dotMultiply(current_c,W._data[3]);
        for(var k = 0; k < Hx.length; k++) {
            Hx[k] = Hx[k].toPolar();
            Hy[k] = Hy[k].toPolar();
        }
        for(var k = 0; k < Hx.length; k++){
            tmpRX[k] = Hx[k].r;
            tmpPhiX[k]= Hx[k].phi;
            tmpRY[k] = Hy[k].r;
            tmpPhiY[k] = Hy[k].phi;            
        }
   
        HxR.push(tmpRX);
        HxPhi.push(tmpPhiX);
        HyR.push(tmpRY);
        HyPhi.push(tmpPhiY);
        
        
        tmpRX = [],tmpPhiX = [];
        tmpRY = [],tmpPhiY = [];
        Hx = [],Hy = [];
        
        
        if (i+1 < numLayers) {
            current_c = c._data.slice(4+4*i,8+4*i);
            
        } 
    }
    
    return {HxR: HxR, HxPhi: HxPhi, HyR: HyR, HyPhi: HyPhi};
};

/**
Calculates the Hx and Hy bits for Mathbox graphing.
Expects the array of layer lengths, the current timestamp t, the position z, 
HxR, HxPhi, HyR, HyPhi.
 */
emScattering2.PhotonicCrystal.prototype.mathboxHf = function(lengths,t,z,HxR,HxPhi,HyR,HyPhi) {
    var layerNum = 0, inf = 0, sup = lengths[0],Hx,Hx1,Hx2,Hx3,Hx4,Hy,Hy1,Hy2,Hy3,Hy4,
            o = this.Struct.omega, lambdas = this.Struct.eigenvalues;
    for(var i = 0; i < lengths.length; i++){
        if(inf <= z && z <= sup){
            layerNum = i;
            break;
        }
        else{
            inf = sup;
            sup = sup + lengths[i + 1];
        }          
    }

    Hx1 = HxR[layerNum][0] * math.exp(math.re(lambdas[layerNum][0]) * z) * math.cos(HxPhi[layerNum][0] + math.im(lambdas[layerNum][0])*z - o*t);
    Hx2 = HxR[layerNum][1] * math.exp(math.re(lambdas[layerNum][1]) * z) * math.cos(HxPhi[layerNum][1] + math.im(lambdas[layerNum][1])*z - o*t);
    Hx3 = HxR[layerNum][2] * math.exp(math.re(lambdas[layerNum][2]) * z) * math.cos(HxPhi[layerNum][2] + math.im(lambdas[layerNum][2])*z - o*t);
    Hx4 = HxR[layerNum][3] * math.exp(math.re(lambdas[layerNum][3]) * z) * math.cos(HxPhi[layerNum][3] + math.im(lambdas[layerNum][3])*z - o*t);
    Hx = Hx1+Hx2+Hx3+Hx4;
   
    Hy1 = HyR[layerNum][0] * math.exp(math.re(lambdas[layerNum][0]) * z) * math.cos(HyPhi[layerNum][0] + math.im(lambdas[layerNum][0])*z - o*t);
    Hy2 = HyR[layerNum][1] * math.exp(math.re(lambdas[layerNum][1]) * z) * math.cos(HyPhi[layerNum][1] + math.im(lambdas[layerNum][1])*z - o*t);
    Hy3 = HyR[layerNum][2] * math.exp(math.re(lambdas[layerNum][2]) * z) * math.cos(HyPhi[layerNum][2] + math.im(lambdas[layerNum][2])*z - o*t);
    Hy4 = HyR[layerNum][3] * math.exp(math.re(lambdas[layerNum][3]) * z) * math.cos(HyPhi[layerNum][3] + math.im(lambdas[layerNum][3])*z - o*t);
    Hy = Hy1+Hy2+Hy3+Hy4;
    return {Hx: Hx, Hy: Hy};
};


/**
Returns the z-coordiantes of the edges of the interfaces. The edges of the ambient mediums 
are determined by the user entered length of the ambeint mediums even though they do extend 
to infinity. The hard limit is more for visual purposes and not because of any underlying 
feature of the structure. The first element is the leftmost extent of the ambient material 
on the left, and the last element is the rightmost extent of the ambient material on the 
right. The second to last element is the rightmost extent of the last layer before the 
ambient material on the right of the strucutre. The ith layer's leftmost extent is in 
the ith position of the returned array. The returned coordinates correspond the the same 
coordinate system used to plot the field values and in other methods of this object. 
*/
emScattering2.PhotonicCrystal.prototype.materialInterfaces = function() {
    var interfaces = new Array();
    interfaces.push(0);
    for (var i = 0; i < this.Struct.numLayers; i++){
        interfaces.push(interfaces[i] + this.Struct.layers[i].length);
    }
    //console.log(interfaces);
    return interfaces;
};



//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================
//=====================================================================================================================================




/** Calls all structure functions in proper order and fills all structure variables, returns filled structure */
emScattering2.computeStructure = function(eArray, mArray, length, numLayers, constants, Modes){
    var values,struct;
    Modes = math.complex(Modes);
    values = emScattering2.Parse(eArray,mArray,numLayers);
    struct = new emScattering2.makeStructure(values[0],values[1],length,constants[0],constants[1],constants[2], Modes);
    struct.makeMaxwell();
    struct.calcEigensystems();
    struct.calcTransfer();
    struct.calcScattering();
    return struct;
};    

//emScattering2.vectorL2Norm = function(vector){
//    var norm = 0 ;
//   
//    if(vector._data[0][0].__proto__.isComplex){
//        for(var i = 0; i < vector._size[0]; i++){
//            norm = norm + math.re(math.multiply(vector._data[i][0],math.conj(vector._data[i][0])));
//        }
//    }
//    else{
//        for(var i = 0; i < vector._size[0]; i++){
//            norm = norm + math.pow(vector._data[i][0],2);         
//        }
//    }
//    return math.sqrt(norm);
//}
//
//emScattering2.unMinor = function(minor,size){
//    var matrix = math.zeros(size,size);
//    var minor_size = minor._size[0];
//    console.log(minor);
//    for(var i = size - minor_size; i < size; i++){
//        for(var j = size - minor_size; j < size; j++)
//            matrix._data[i][j] = minor._data[i-size+minor_size][j-size+minor_size];
//    }
//    console.log(matrix);
//    return matrix;
//}
//
//emScattering2.computeQ = function(vector,e,size){
//    var u,v,transpose,Q;
//    u = math.subtract(vector,e);
//    v = math.divide(u,emScattering2.vectorL2Norm(u));
//    if(v._data[0][0].__proto__.isComplex) {transpose = math.ctranspose(v);}
//    else{transpose = math.transpose(v);}
//    Q = math.subtract(math.identity(size),math.multiply(2,math.multiply(v,transpose)));
//    return Q;
//}
//
//emScattering2.houseHolder = function(matrix){
//    var size = matrix._size[0], index = [], u ,e = [],v, x, transpose,Q,M;
//    
//    for (var i = 0; i < size; i++) {index.push(i); e.push([0]);}
//    
//    x = math.subset(matrix,math.index(index,0));
//    e[0] = [emScattering2.vectorL2Norm(x)];
//    if(x._data[0][0].__proto__.isComplex){
//        e[0][0] = math.multiply(math.unaryMinus(math.exp(math.multiply(math.complex("i") , math.arg(x._data[0][0])))),e[0][0]);
//    }
//    Q = emScattering2.computeQ(x,e,size);
//    M = math.multiply(Q,matrix); 
//   
//    for(var i = 0; i < size; i++){
//        for(var j = 0; j < size; j++){
//            if(math.abs(M._data[i][j]) < math.pow(10,-14)){M._data[i][j] = 0;}
//        }
//    }
//    return [M,Q]; 
//};
//
//emScattering2.QR = function(matrix){
//    var setup = true, size = matrix._size[0],M, Q = [], tmp, minor; 
//    for(var i = 0; i < size - 1; i++){
//        if(setup){
//            tmp = emScattering2.houseHolder(matrix);
//            setup = false;
//        }
//        else{
//            var index_col = [], index_row = [] ;
//            for(var j = 0; j < size - i; j++ ){
//                index_row.push(j + i);
//                index_col.push(j+i);
//            }
//            minor = math.subset(M,math.index(index_row,index_col));
//            tmp = emScattering2.houseHolder(minor);
//            tmp[1] = emScattering2.unMinor(tmp[1],size);
//        }
//        M = tmp[0];
//        Q.push(tmp[1]);
//    }
//    console.log(Q);
//}

/**Performs experiment when run, returns completed crystal*/
emScattering2.Driver = function(eArray, mArray, length, numLayers,constants,Modes){
    var a,b,c ,d = [],struct, crystal;
    struct = emScattering2.computeStructure(eArray, mArray, length, numLayers, constants, Modes);
    crystal = emScattering2.createPhotonicCrystal(struct);
    //crystal.determineField();  
    //console.log(crystal);
//    c = math.complex("1 + i");
//    a = math.matrix([[12,-51,4],
//                    [6,167,-68],
//                    [-4,24,-41]]);
//    b = math.matrix([[c,c,c,c],
//                     [c,c,c,c],
//                     [c,c,c,c],
//                     [c,c,c,c]]);
//    b._datatype = "complex";
//    
//    emScattering2.QR(a);
    
    
    return crystal;    
    
    
    
};

/**CreateTransmissionArrays
 * --------------------
 * Takes a range of Omega values (omegaLow-omegaHigh) and increments them based on n (omegaPoints) times. Calculates the crystal and fields for
 * each calculated value of Omega. Takes the Ex, Ey, Hx, and Hy values at point z and puts them into an array. If Omega is zero, no calculations will be done
 * and the method will continue. WIP - Extremely inefficient and can take minutes based on the value of omegaPoints
 * 
 */
emScattering2.createTransmissionArrays = function(eArray, mArray, length, numLayers, k1, k2, modes, omegaLow, omegaHigh, omegaPoints, zPoint) {
    var _omegas = new Array(), _Ex = new Array(), _Ey = new Array(), _Hx = new Array(), _Hy = new Array();

    omegaHigh = Number(omegaHigh);
    omegaLow = Number(omegaLow);
    omegaPoints = Number(omegaPoints);
    zPoint = Number(zPoint);

    var omegaInterval = (omegaHigh - omegaLow) / omegaPoints;

    //console.log(omegaInterval);

    for(var i = 0; i <= omegaPoints; i++) {
        var tempOmega = omegaLow + (omegaInterval * i);

        //console.log(tempOmega);

        if(tempOmega == 0) continue;

        _omegas.push(tempOmega);

        var constants = [k1, k2, tempOmega];
        var crystal = this.Driver(eArray, mArray, length, numLayers, constants, modes);

        checkBoxesForModes(crystal);                      //Disabled until we can come up with an algoritm to correct modes throughout all crystals

        //DetermineField Method
        /*
        var zIndex = zPoint * 100;

        var interfaces = crystal.materialInterfaces();

        if(zPoint == interfaces[interfaces.length - 1]) zIndex--;

        var fields = crystal.determineField();

        console.log(crystal.determineFieldAtZPoint(zPoint));

        _Ex.push(fields.Ex[zIndex]);
        _Ey.push(fields.Hx[zIndex]);
        _Hx.push(fields.Hx[zIndex]);
        _Hy.push(fields.Hy[zIndex]);
        */

        //DetermineFieldAtZPoint Method
        ///*
        var field = crystal.determineFieldAtZPoint(zPoint);

        _Ex.push(field.Ex);
        _Ey.push(field.Ey);
        _Hx.push(field.Hx);
        _Hy.push(field.Hy);
        //*/
    }

    return {omegas: _omegas, Ex: _Ex, Ey: _Ey, Hx: _Hx, Hy: _Hy};
}
/* TODO (Transmission Arrays):
*   - Make method more efficient for higher precision calculations
*   - Hope numbers are correct when issue with calculations is fixed
*/

/**CheckBoxesForModes
 * ------------
 * For the transmission tab. Takes the crystal and rearranges the modes based on the selected check
 * boxes in the "Incoming Modes in Ambient Medium" section.
 * WIP - UNUSED
 */
function checkBoxesForModes(crystal) {
    var backChecked = 0;
    var forChecked = 0;

    for(let i = 1; i <= 4; i++){
        if(document.getElementById("backModeChkT" + i).checked == true) backChecked++;
        if(document.getElementById("forModeChkT" + i).checked == true) forChecked++;
    }

    var backArr = crystal.Struct.eigenvalues[0];
    var backVecArr = crystal.Struct.eigenvectors[0]._data;
    var forArr = crystal.Struct.eigenvalues[crystal.Struct.eigenvalues.length-1];
    var forVecArr = crystal.Struct.eigenvectors[crystal.Struct.eigenvalues.length-1]._data;

    if(backChecked == 2){
        let j = 0;

        for(let i = 1; i <= 4; i++){
            if(document.getElementById("backModeChkT" + i).checked == true){
                backArr.splice(j, 0, backArr.splice(i-1, 1)[0]);            //Swap array around based on order of checked boxes
                backVecArr.splice(j, 0, backVecArr.splice(i-1, 1)[0]);      //Swap eigenvectors around based on order of checked boxes
            }
        }
    }

    if(forChecked == 2){
        let j = 0;

        for(let i = 1; i <= 4; i++){
            if(document.getElementById("forModeChkT" + i).checked == true) {
                forArr.splice(j, 0, forArr.splice(i-1, 1)[0]);              //Swap array around based on order of checked boxes
                forVecArr.splice(j, 0, forVecArr.splice(i-1, 1)[0]);        //Swap eigenvectors around based on order of checked boxes
            }
        }
    }
}