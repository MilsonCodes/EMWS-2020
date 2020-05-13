// Global namespace
var emScattering = emScattering || {};


// Intrinsic impedance of free space in Ohms
emScattering.ETTA_0 = 376.73031;

// Speed of light in free space in m/s -- CURRENTLY IN UNITS WHERE THE NUMERICAL VALUE IS 1
// emScattering.C = 299792458.0; // m/s
emScattering.C = 1;

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
emScattering.Layer = function(epsilon, mu, length) {
    this.epsilon = epsilon;
    this.mu = mu;
    this.length = length;
};

emScattering.Layer.prototype.toString = function() {
    console.log("Layer: " + this.epsilon + ", " + this.mu + " ," + this.length);
};


/* 

GENERAL FUNCTIONS 
------------------------------------------------------------------------------------
*/

/**
Prints out the numeric array to the console. 'pre' is a string that is printed before 'mat'.
*/
emScattering.print = function(pre, mat) {
    document.write(pre + numeric.prettyPrint(mat) + "<br>");
};

/**
Prints out to the console a large array in an easily copied format.
*/
emScattering.printFields = function(fields) {
    console.log("z Ex Ey Hx Hy");
    for (var i = 0, N = fields.z.length; i < N; i++) {
        console.log(fields.z[i]+" "+fields.Ex[i]+" "+fields.Ey[i]+" "+fields.Hx[i]
            +" "+fields.Hy[i]);
    }
};

/**
Prints out a large array in an easily copied format.
*/
emScattering.printFields2 = function(fields) {
    document.write("<table>");
    document.write("<tr><td>z</td> <td>Ex</td> <td>Ey</td> <td>Hx</td> <td>Hy</td></tr>");
    for (var i = 0, N = fields.z.length; i < N; i++) {
        document.write("<tr>");
        document.write("<td>" + fields.z[i] + "</td>");
        document.write("<td>" + fields.Ex[i] + "</td>");
        document.write("<td>" + fields.Ey[i] + "</td>");
        document.write("<td>" + fields.Hx[i] + "</td>");
        document.write("<td>" + fields.Hy[i] + "</td>");
        document.write("</tr>");
    }
    document.write("</table>");
};

emScattering.printDispersion = function(dispersion) {
    document.write("<table><tr><td>kz</td>");
    for (var i = 0; i < dispersion.layersDispersions.length; i++)
        document.write("<td>Layer " + i + "</td>")
    document.write("</tr>");
    
    for (i = 0; i < dispersion.kz.length; i++) {
        document.write("<tr>");
        document.write("<td>" + dispersion.kz[i] + "</td>");
        for (var j = 0; j < dispersion.layersDispersions.length; j++)
            document.write("<td>" + dispersion.layersDispersions[j][i] + "</td>");
        document.write("</tr>");
    }
};

emScattering.printModes = function(modes) {
    for (var i = 0; i < modes.length; i++) {
        document.write(numeric.prettyPrint(modes[i].eigenvalue) + " ; ");
        document.write(numeric.prettyPrint(modes[i].eigenvector) + " ; ");
        document.write(numeric.prettyPrint(modes[i].forward) + " ; ");
        document.write("<br>");
        document.write("<br>");
    }
};

emScattering.printTransmissionGraph = function(tGraph) {
    document.write("<table><tr><td>omega</td>");
    for (var i = 0; i < tGraph.kzList.length; i++)
        document.write("<td>kz_3=" + tGraph.kzList[i] + "</td>")
    document.write("</tr>");

    for (i = 0; i < tGraph.omegaRange.length; i++) {
        document.write("<tr>");
        document.write("<td>" + tGraph.omegaRange[i] + "</td>");
        for (var j = 0; j < tGraph.kzList.length; j++)
            document.write("<td>" + tGraph.transmissionCoeffArrays[j][i] + "</td>");
        document.write("</tr>");
    }
};

/**
Returns an nxn matrix with every element set to 0.
*/
emScattering.zeroMatrix = function(n) {
    return numeric.T.rep([n, n], new numeric.T(0,0));
};

/**
Converts an array of epsilon, mu, and length values into an array of layers where 
layers[i] has the value of epsilon[i], mu[i], and length[i].
*/
emScattering.createLayers = function(epsilon, mu, lengths) {
    var layers = [];
    for (var i = 0, N = epsilon.length; i < N; i++) {
        var layer = new emScattering.Layer(epsilon[i], mu[i], lengths[i]);
        layers.push(layer);
    }
    return layers;
};

/**
Calculates the eigenvector matrix for a layer of isotropic material. 'eps' is the 
relative permittivity and 'mu' us the relative permeability of the material. The first 
two columns are the forward propogating modes, and the last two columns are the 
backwards propogating modes. 'kx' and 'ky' are the transverse wavenumbers. For normal incidence,
they are 0.
*/
emScattering.eigenvectorsIsotropic = function(epsilon, mu, kx, ky) {
    var W = emScattering.zeroMatrix(4);
    var forwardEig = emScattering.eigenvaluesIsotropic(epsilon, mu, kx, ky)[0];  
    var element1, element2, element3;
    
    // Convenience values for constructing the matrix
    element1 = new numeric.T(kx*ky/mu,0);
    element1 = element1.div(forwardEig);
    
    element2 = new numeric.T((mu*epsilon - Math.pow(kx, 2))/mu,0);
    element2 = element2.div(forwardEig);
    
    element3 = new numeric.T((mu*epsilon - Math.pow(ky, 2))/mu,0);
    element3 = element3.div(forwardEig); 
    
    // Eigenvector 1 - forward prop
    W.set([0, 0], new numeric.T(1, 0));
    W.set([2, 0], element1);
    W.set([3, 0], element3.dot(-1));
    
    // Eigenvector 2 - forward prop
    W.set([1, 1], new numeric.T(1, 0));
    W.set([2, 1], element2);
    W.set([3, 1], element1.dot(-1));
    
    // Eigenvector 3 - backward prop
    W.set([0, 2], new numeric.T(1, 0));
    W.set([2, 2], element1.dot(-1));
    W.set([3, 2], element3);
    
    // Eigenvector 4 - backward prop
    W.set([1, 3], new numeric.T(1, 0));
    W.set([2, 3], element2.dot(-1));    
    W.set([3, 3], element1);
    
    return W;
};

/**
Returns an array containing the eigenvalues for a layer of isotropic material for 
a plane wave at oblique incidence. The eigenvalues correspond to the ordering of the 
eigenvectors returned by emScattering.eigenvectorsIsotropic(). Use emScattering.diag() 
to create a diagonal matrix of the eigenvalues, if necessary. 
*/
emScattering.eigenvaluesIsotropic = function(epsilon, mu, kx, ky) {
    var N = (Math.pow(kx*ky, 2) - (mu*epsilon-Math.pow(kx, 2))*(epsilon*mu-Math.pow(ky, 2)))/(mu*epsilon);
    var forwardEig, backwardEig;
    if (N >= 0) {
        N = Math.sqrt(N);
        forwardEig = new numeric.T(N, 0);
        backwardEig = forwardEig.dot(-1);
    }
    else if (N<0){
        N = Math.sqrt(-1*N);
        forwardEig = new numeric.T(0, N);
        backwardEig = forwardEig.dot(-1);
    }
    return [forwardEig, forwardEig, backwardEig, backwardEig];
};

/**
Returns a diagonal matrix with the eigenvalues for an isotropic material layer.
*/
emScattering.eigenvaluesIsotropicDiag = function(epsilon, mu, kx, ky) {
    var eigs = emScattering.eigenvaluesIsotropic(epsilon, mu, kx, ky);        
    var lambda = emScattering.zeroMatrix(4);
    eigs.forEach(function(item, i, array) {
        lambda.set([i, i], item);
    });
    return lambda;
};

/**
Returns exp(lambda*znorm) where 'lambda' is a diagonal matrix of the eigenvalues for a 
isotropic material layer. 'znorm' is the normalized z coordinate in the layer. 
znorm = 0 is the left side of the layer, and znorm = L is the right side of the layer 
where L is the length of the layer. Lambda is assumed to have the same ordering as that 
returned by emScattering.eigenvaluesIsotropicDiag() and 
emScattering.eigenvaluesIsotropic().
*/
emScattering.expEigenvaluesIsotropicDiag = function(lambda, znorm) {
    var diag = lambda.getDiag();
    var diagX = diag.x, diagY = diag.y, x, y;
    x = numeric.mul(diagX, [znorm, znorm, znorm, znorm]);
    y = numeric.mul(diagY, [znorm, znorm, znorm, znorm]);
    zPrime = new numeric.T(x,y);
    return numeric.T.diag(zPrime.exp());
};

emScattering.EigenVal2 = function(lambda, theta, k0, Zed, Time){
	var diag = lambda.getDiag();
    var diagX = diag.x, diagY = diag.y, x, y;
    var xprime = numeric.mul(diagX, [Zed, Zed, Zed, Zed]);
    var yprime = numeric.mul(diagY, [Zed, Zed, Zed, Zed]);
    var xprime2 = numeric.add(xprime,theta);
    var yprime2 = numeric.add(yprime,theta);
    x = numeric.sub(xprime2,k0*Time);
    y = numeric.sub(yprime2,k0*Time);
    zPrime = new numeric.T(x,y);
    return numeric.T.diag(zPrime.cos());
};

/**
Returns a normalized value z' for a given 'z'. 'z' is a numeric vector.
z' = k_0 * z where k_0 is the freespace wavelength.
*/
emScattering.normalizeZ = function(z, k0) {
    return z.mul(numeric.rep(numeric.dim(z), new numeric.T(k0, 0)));
};


/*
Structure Object
------------------------------------------------------------------------------------
*/

/**
A structure is a 1D photonic crystal. It is composed of several layers each with
their own epsilon, mu, and length. The layers define the scattering properties of the
crystal, and functions are provided to solve the scattering problem and determine the
field within the structure. Includes ambient materials to the left and the right, which
are of infinite extent. The first an last elements in 'layers' are for the ambient medium
on either side of the structure. The 'length' property of each of these layers is ignored. 
*/
emScattering.Structure = function(layers) {
    this.layers = layers;
    this.numLayers = layers.length;
    
    this.transferMatrices = new Array(this.numLayers-1);
    this.eigenvectors = new Array(this.numLayers);
    this.eigenvalues = new Array(this.numLayers);
    
    this.generalTransferMatrix = numeric.T.identity(4);
};

emScattering.Structure.prototype.calcEigenvectorsIsotropic = function(kx, ky) {
    for (var i = 0; i < this.numLayers; i++) {
        var layer = this.layers[i];
        this.eigenvectors[i] = emScattering.eigenvectorsIsotropic(layer.epsilon, layer.mu, kx, ky);
    }
};

emScattering.Structure.prototype.calcEigenvaluesIsotropic = function(kx, ky) {
    for (var i = 0; i < this.numLayers; i++) {
        var layer = this.layers[i];
        this.eigenvalues[i] = emScattering.eigenvaluesIsotropicDiag(layer.epsilon, layer.mu, kx, ky);
    }
};

/**
Creates the transfer matrices for the given free space wavelength k_0. 
*/
emScattering.Structure.prototype.calcTransferMatrices = function(k_0, kx, ky) {
    this.transferMatrices = new Array(this.numLayers-1);
    this.calcEigenvectorsIsotropic(kx/k_0, ky/k_0);
    this.calcEigenvaluesIsotropic(kx/k_0, ky/k_0);
    
    // Calculates the transfer matrices between each layer in the structure
    // Number of Layers - 1 = number of Interfaces
    if (this.numLayers > 1) {
        for (var i = 0, N = this.transferMatrices.length; i < N; i++) {
            var wNext = this.eigenvectors[i+1].inv();
            var w = this.eigenvectors[i];
            var znorm;
            if (i == N-1 || i == 0)
                znorm = 0;
            else
                znorm = this.layers[i].length * k_0;
            var expLambda = emScattering.expEigenvaluesIsotropicDiag(this.eigenvalues[i], znorm);
            this.transferMatrices[i] = wNext.dot(w.dot(expLambda));
        }
    }
    //this.transferMatrices[0] = this.eigenvectors[1].inv().dot(this.eigenvectors[0]);
    
    // General transfer matrix 
    for (i = 0; i < this.numLayers-1; i++) {
        this.generalTransferMatrix = this.transferMatrices[i].dot(this.generalTransferMatrix);
    }
};

/**
Creates the scattering matrix for the given free space wavelength k_0.
If U- and U+ represent the backwards and forewards propogating modes that travel 
away from the structure and J+ and J- represent the forewards and backwards propogating 
modes the travel towards the structure, the [U- U+] = S * [J+ J-]. S is a 4x4 matrix, and 
the Us and Js are 2x1 vectors.
J+ = [J+_1 J+_2] -- the two forward propogating incoming modes on the left of the strucutre
J- = [J+_1 J+_2] -- the two backwards propogating incoming modes on the right of the structure
U- = [U-_1 U-_2] -- the two backwards propogating outgoing modes on the left of the structure
U+ = [U+_1 U+_2] -- the two forward propogating outgoing modes on the right of the structure
*/
emScattering.Structure.prototype.calcScatteringMatrix = function(k_0, kx, ky) {
    this.calcTransferMatrices(k_0, kx, ky);

    // Block matrices composing the general transfer matrix T = [T11 T12; T21 T22]
    var T11 = this.generalTransferMatrix.getBlock([0,0], [1,1]);
    var T12 = this.generalTransferMatrix.getBlock([0,2], [1,3]);
    var T21 = this.generalTransferMatrix.getBlock([2,0], [3,1]);
    var T22 = this.generalTransferMatrix.getBlock([2,2], [3,3]);
    
    var iT22 = T22.inv();
    
    // The scattering matrix's blocks
    var S11 = iT22.dot(T21).mul(numeric.T.rep([2,2], new numeric.T(-1,0)));
    var S12 = iT22;
    var S21 = T11.sub(T12.dot(iT22.dot(T21)));
    var S22 = T12.dot(iT22);
    
    // emScattering.zeroMatrix() isn't used because setBlock() has issues with a complex valued 
    // zero matrix. The final two setBlocks() fail if emScattering.zeroMatrix() (or if numeric
    // is used directly to create the same thing), so a matrix of nonzero complex numbers (both 
    // the imaginary and the real part must be nonzero) is created. The actual complex numbers 
    // don't matter because they will all be overwritten by the setBlock() calls.
    var S = new numeric.T.rep([4,4], new numeric.T(1,1));
    S.setBlock([0,0], [1,1], S11);
    S.setBlock([0,2], [1,3], S12);
    S.setBlock([2,0], [3,1], S21);
    S.setBlock([2,2], [3,3], S22);
    
    return S;
};



/**
Object representing a scattering experiment on a fixed structure. Provides functions to 
calculate the fields and dispersion relationship.
------------------------------------------------------------------------------------
*/
emScattering.PhotonicStructure1D = function(epsilon, mu, length) {
    this.layers = emScattering.createLayers(epsilon, mu, length);
    this.crystal = new emScattering.Structure(this.layers);
    this.N = this.layers.length;
    this.k0 = -1;   // TODO: Use to cache the results
    this.S = [];
};

/**
Solves the scattering problem at the given free space wavelength k0 and the specified 
incoming mode which are an Array J = [J+_1, J+_2, J-_1, J-_2] are the field amplitudes of the 
2 incoming modes on both sides of the structure for a total of 4 incoming modes. Returns 
U = [U-_1, U-_2, U+_1, U+_2]
*/
emScattering.PhotonicStructure1D.prototype.scattering = function(k0, kx, ky, J) {
    this.S = this.crystal.calcScatteringMatrix(k0, kx, ky);
    return this.S.dot(J);
};
// currently only building Ex lists
emScattering.PhotonicStructure1D.prototype.generateLists = function(k0, kx, ky, J) {
	console.log("go1");
	var N = this.crystal.numLayers;
	var r_Ex = new Array(), theta_Ex = new Array();
	var U = this.scattering(k0, kx, ky, J);
	console.log("go2");
	for(i=0; i<N; i++){
		W = this.crystal.eigenvectors[i];

	}
	return{RListEx: r_Ex, ThetaListEx: theta_Ex};

};


emScattering.PhotonicStructure1D.prototype.MathboxStuff = function(k0, kx, ky, J, Zed, Time){
	console.log("ho");
    var N = this.crystal.numLayers;
	var _Ex = new Array(), _Ey = new Array(), _Hx = new Array(), _Hy = new Array();

	var U = this.scattering(k0, kx, ky, J);

	var cx = [J[0], J[1], U.x[0], U.x[1]];
    var cy = [0, 0, U.y[0], U.y[1]];
    var c = new numeric.T(cx, cy);   // Coefficients in the current layer

    var W, lambda, result, i, j, cosDiag;
    console.log("ho1");
    for(i=0; i<N; i++){
	    W = this.crystal.eigenvectors[i];
	    lambda = this.crystal.eigenvalues[i];
	    console.log("ho2");
	    var A = W.dot(c.x);
	    var B = W.dot(c.y);
	    var R = Math.sqrt(Math.pow(A,2)+Math.pow(B,2));
		console.log("ho3");
		var theta = Math.atan2(B,A);
		console.log("ho4");
		cosDiag = emScattering.EigenVal2(lambda, theta, k0, Zed, Time);
		console.log("ho5");

		var ans = R*cosDiag;
		console.log("ho6");
		if (i+1 < N) {
            c = this.crystal.transferMatrices[i].dot(c); //updates coefficients for next layer
        } 
        console.log("ho7");
        _Ex.push(ans[0]);
        _Ey.push(ans[1]);
        _Hx.push(ans[2]);
        _Hy.push(ans[3]);
        console.log("ho8");
	}
    
    return {Ex: _Ex, Ey: _Ey, Hx: _Hx, Hy: _Hy};

};

emScattering.PhotonicStructure1D.prototype.checkZMathBox = function(z){
    var indexOfMatrix = 0;
    var prepLayers = this.prepareLayer();
        for(i=0; i< prepLayers.length - 1; i++){
        	var ind = i;
            if((prepLayers[i] <= z) && (z <= prepLayers[i+1])){
                indexOfMatrix++;
                break;
            }
        }
    return indexOfMatrix;
};
/**returns the layers in a list compatible for checking in checkZMathBox*/
emScattering.PhotonicStructure1D.prototype.prepareLayer = function(){
    var prepLayer = new Array();
    var totalLength = 0;
    prepLayer.push(totalLength);
    for(i=0; i < this.crystal.layers.length; i++){
            totalLength = totalLength + this.crystal.layers[i].length;
            prepLayer.push(totalLength);
    }
    return prepLayer;
}
/**
Returns the field values in all the layers given the coefficients of the incoming modes. The 
fields should be ordered so that J = [J+_1, J+_2, J-_1, J-_2] where F is a 4x1 array or 
numeric marix. Returned object has properties z for the coordinates used, Ex, Ey, Hx, 
and Hy. There is a one-to-one correspondence between an element in z and the other 4 Arrays.
*/
emScattering.PhotonicStructure1D.prototype.determineField = function(k0, kx, ky, J) {
    var N = this.crystal.numLayers;
    var numPoints = 100;    // Number of points calculated per unit length
    var _Ex = new Array(), _Ey = new Array(), _Hx = new Array(), _Hy = new Array(), _z = new Array();
     
    var U = this.scattering(k0, kx, ky, J);
    
    var cx = [J[0], J[1], U.x[0], U.x[1]];
    var cy = [0, 0, U.y[0], U.y[1]];
    var c = new numeric.T(cx, cy);   // Coefficients in the current layer
    
    // Determine the fields in each layer
    var W, lambda, layerNormZ, result, i, j, expDiag;
    var currentLeftZ = 0, currentRightZ = this.crystal.layers[0].length;
    for (i = 0; i < N; i++) {
        layerNormZ = numeric.linspace(0, k0*this.crystal.layers[i].length, 
                                this.crystal.layers[i].length*numPoints);
        W = this.crystal.eigenvectors[i];
        lambda = this.crystal.eigenvalues[i];
        
        for (j = 0, M = layerNormZ.length; j < M; j++) {
            expDiag = emScattering.expEigenvaluesIsotropicDiag(lambda, layerNormZ[j]);
            result = W.dot(expDiag.dot(c));
            _Ex.push(result.x[0]);
            _Ey.push(result.x[1]);
            _Hx.push(result.x[2]);
            _Hy.push(result.x[3]);
        }
        
        _z = _z.concat(numeric.linspace(currentLeftZ, currentRightZ, 
                        numPoints*this.crystal.layers[i].length));
        if (i+1 < N) {
            currentLeftZ += this.crystal.layers[i].length;
            currentRightZ += this.crystal.layers[i+1].length;
            c = this.crystal.transferMatrices[i].dot(c);
        } 
    }
    return {z: _z, Ex: _Ex, Ey: _Ey, Hx: _Hx, Hy: _Hy};
};

/**
Returns the dispersion relationship for each layer and the ambient mediums on the left 
and right side of the strucutre. The relation is calculated over the free space wavelength 
range [-khi, khi] 'numPoints' is the number of in the interval for k0 that the relation is 
calculated for. 
The return value is an object with two fields:
dispersion.k0: the range of k0's that the dispersion relation was calcualted for 
dispersion.layersDispersions: an array of the calculated k_z values for the k_0 values 
                              for each layer. layersDispersions[0] is the left ambient 
                              and layersDispersions[n+1] is the right ambient. The array 
                              is length n+2 where n is the number of layers.
*/
// TODO: Return an object with methods to access data instead of one with only properties 
// where the user must use the hardcoded names of the properties to access the data
// TODO: Return as a function of kz instead of k0 (aka the frequency) -- the reverse of the current 
emScattering.PhotonicStructure1D.prototype.dispersionRelationship = function(kx, ky, khi, numPoints) {
    var _kz = new Array();
    var kz_pos = numeric.linspace(0, khi, Math.ceil(numPoints/2));
    for (var i = kz_pos.length-1; i > 0; i--)
        _kz.push(-1*kz_pos[i]);
        
    for (i = 0; i < kz_pos.length; i++)
        _kz.push(kz_pos[i]);
    
    var dispersionRelationships = new Array();
    for (i = 0; i < this.N; i++) {
        var k0 = new Array();
        var den = this.layers[i].mu*this.layers[i].epsilon;
        for (var j = 0; j < _kz.length; j++) {
            var num = Math.pow(_kz[j],2) + Math.pow(kx,2) + Math.pow(ky,2);
            k0.push(emScattering.C * Math.sqrt(num / den));
        }
        dispersionRelationships.push(k0);
    }
    // ACTUALLY RETURNS THE ANGULAR FREQUENCY -- NEEDS TO ALSO BE CHANGED IN GUI CODE
    // TODO: Return a DispersionRelationGraph object similar to how transmission() works
    return {kz: _kz, layersDispersions: dispersionRelationships};
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
emScattering.PhotonicStructure1D.prototype.materialInterfaces = function() {
    var interfaces = new Array();
    interfaces.push(0);
    for (var i = 0; i < this.N; i++)
        interfaces.push(interfaces[i] + this.layers[i].length);
    
    return interfaces;
};

/**
Returns the transmission vs. frequency graph object. For the specified k_z values (along with the k_x and k_y values), the
function returns the transmission for each of the k_z values over the specified omega range.
 */
emScattering.PhotonicStructure1D.prototype.transmission = function(kzList, kx, ky, omegaLo, omegaHi, numPoints)  {
    omegaRange = numeric.linspace(omegaLo, omegaHi, numPoints);

    var ambientLeftEigenvalue, ambientRightEigenvalue;
    var transmissionCoeffArrays = new Array();

    for (var i = 0; i < kzList.length; i++) {
        var num = Math.pow(kzList[i],2) + Math.pow(kx,2) + Math.pow(ky,2);
        var den = this.layers[0].mu*this.layers[0].epsilon;
        var omegaMin = emScattering.C * Math.sqrt(num/den);

        var transmissionCoeff = new Array();
        for (var j = 0, N = omegaRange.length; j < N; j++) {
            if (omegaRange[j] >= omegaMin) {
                var S = this.crystal.calcScatteringMatrix(omegaRange[j] / emScattering.C, kx, ky);
                var transmission = Math.pow(S.get([2, 0]).norm2(), 2);
                //emScattering.print("omega = " + omegaRange[j] + ", S = ", S)
                transmissionCoeff.push(transmission);
            }
            else {
                transmissionCoeff.push(0);
            }
        }
        transmissionCoeffArrays.push(transmissionCoeff);
    }

    return new emScattering.TransmissionCoeffGraph(kzList, kx, ky, omegaRange, transmissionCoeffArrays);

};

/**
Returns the wavelengths of the incoming modes supported by the structure. As this is an isotropic, 
1D structure, there will be 4 modes in total. Two forward travelling modes coming from the left 
of the structure (travelling to the right), and two backwards travelling modes coming from the 
right of the structure (travelling to the left). The array will be ordered so that the first two 
are in the ambient material on the left of the structure. The last two are in the ambient material 
on the right of the structure. 
*/
// TODO: Change for normal incidence. Deduce incoming from the transfer matrices.
emScattering.PhotonicStructure1D.prototype.incomingModes = function(k0,kx,ky) {
    this.scattering(k0, kx, ky, [0,0,0,0]);
    return emScattering.Mode1D.createModes(this.crystal.eigenvalues[0],this.crystal.eigenvectors[0]);
};


/**
 TransmissionCoeffGraph Object
 ------------------------------------------------------------------------------------
 Represents the transmission vs. frequency for a list of specified kz values at set kx and ky values over a range of
 frequencies. Able to query the transmission for a given wavenumber (if it is available).
 */
emScattering.TransmissionCoeffGraph = function(kzList, kx, ky, omegaRange, transmissionCoeffArrays) {
    this.kzList = kzList;
    this.kx = kx;
    this.ky = ky;
    this.omegaRange = omegaRange;
    this.transmissionCoeffArrays = transmissionCoeffArrays;
};

emScattering.TransmissionCoeffGraph.prototype.getOmegaRange = function() {
    return this.omegaRange;
};

emScattering.TransmissionCoeffGraph.prototype.getTransmission = function(kz) {
    if (this.kzInGraph(kz))
        return this.transmissionCoeffArrays(kzList.indexOf(kz));
    else
        return [];
};

emScattering.TransmissionCoeffGraph.prototype.k1 = function() {
    return this.kx;
};

emScattering.TransmissionCoeffGraph.prototype.k2 = function() {
    return this.ky;
};

emScattering.TransmissionCoeffGraph.prototype.kzList = function() {
    return this.kzList;
};

emScattering.TransmissionCoeffGraph.prototype.kzInGraph = function(kz) {
    return this.kzList.indexOf(kz) != -1;
};

/**
Mode Object
------------------------------------------------------------------------------------
Represents a 1D mode. Contains the eigenvalue, eigenvector, and if it is forward travelling or not.
*/
emScattering.Mode1D = function(eigenvalue, eigenvector) {
    this.eigenvalue = eigenvalue;
    this.eigenvector = eigenvector;
    this.forward = emScattering.isForwardPropogatingMode(eigenvector);
};

emScattering.Mode1D.prototype.wavenumber = function() {
    var wn = new numeric.T(0,-1);
    return wn.dot(eigenvalue);
};

emScattering.Mode1D.prototype.modeVector = function() {
    return eigenvector;
};

emScattering.Mode1D.createModes = function(eigenvalues, eigenvectors) {
        var modes = new Array(4);
        for (var i = 0; i < 4; i++)
            modes[i] = new emScattering.Mode1D(eigenvalues.get([i,i]), eigenvectors.getBlock([0,i], [3,i]));
        return modes;
};

// TODO: Move calculation for Poynting vector into its own function
emScattering.isForwardPropogatingMode = function(eigenvector) {
    var coeff, Ex, Ey, Hx, Hy, poynting;
    coeff = new numeric.T(0, 1 / emScattering.ETTA_0);
    Ex = eigenvector.get([0,0]);
    Ey = eigenvector.get([1,0]);
    Hx = eigenvector.get([2,0]);
    Hy = eigenvector.get([3,0]);
    poynting = Ex.mul(Hy).sub(Ey.mul(Hx));
    poynting = coeff.dot(poynting);
    //emScattering.print("",eigenvector);
    //emScattering.print("",poynting);
    return poynting.x > 0;
};


/*
For testing.
------------------------------------------------------------------------------------
*/
//epsilon = [1,5,10,5,1];
//mu = [1,10,1,10,1];
//lengths = [6,3,3,3,6];

//crystal = new emScattering.PhotonicStructure1D(epsilon, mu, lengths);

//transmissionGraph = crystal.transmission([1,2,3], 0, 0, 0, 5, 100);
//emScattering.printTransmissionGraph(transmissionGraph);

//fields = crystal.determineField(5, 0, 0, [1,0,0,0]);
//emScattering.printFields2(fields);

//S = crystal.crystal.calcScatteringMatrix(5 / emScattering.C, 0, 0);
//transmission = Math.pow(S.get([2, 0]).norm2(), 2);
//emScattering.print("", transmission);

//dispersion = crystal.dispersionRelationship(1, 1, 4, 100);
//emScattering.printDispersion(dispersion);

//interfaces = crystal.materialInterfaces();
//emScattering.print("", interfaces);

//incoming = crystal.incomingModes(1,.2,.4);
//emScattering.printModes(incoming);


