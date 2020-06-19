import math
import time
import numpy as np
from scipy.linalg import null_space, solve, lstsq
from numpy.linalg import inv

# Set precision for printing arrays
np.set_printoptions(precision=6, suppress=True)

def expDiagonal(eigVals, zNorm):
    tmp = [eigVals[0], eigVals[1], eigVals[2], eigVals[3]]
    return np.diag(np.exp(np.multiply(tmp,zNorm)))

def swapArrayIndices(a, i, j):
    #print('Swapping ' + str(a[i]) + ' and ' + str(a[j]))
    a[i], a[j] = a[j], a[i]
    return a

def swapMatrixRows(a, i, j):
    #a[[i, j]] = a[[j,i]]
    #print('row swap ' + str(i) + ' and ' + str(j))
    a[:,[i,j]] = a[:,[j,i]]
    return a

def isNumZero(num):
    precision = 1e-6
    return (num < precision and num > -precision) or num == 0

def organizeEigen(val, vec):
    #print(val)
    for i in range(4):
        v = val[i]
        if (not isNumZero(v.imag)) and (isNumZero(v.real)):
            #print('Eigenvalue ' + str(v) + ' is imaginary')
            if (v.imag > 0) and (not i == 0):
                #print('Eigenvalue ' + str(v) + ' is positive imaginary. Swapping to front!')
                val = swapArrayIndices(val, i, 0)
                vec = swapMatrixRows(vec, i, 0)
            #elif (v.imag < 0) and (not i == 2):
            #    print('Eigenvalue ' + str(v) + ' is negative imaginary. Swapping to third entry!')
            #    val = swapArrayIndices(val, i, 2)
            #    vec = swapMatrixRows(vec, i, 2)
        elif (isNumZero(v.imag)) and (not isNumZero(v.real)):
            #print('Eigenvalue ' + str(v) + ' is real')
            if (v.real < 0) and (not i == 1):
                #print('Eigenvalue ' + str(v) + ' is negative real. Swapping to second index!')
                val = swapArrayIndices(val, i, 1)
                vec = swapMatrixRows(vec, i, 1)
            #elif (v.real > 0) and (not i == 3):
            #    print('Eigenvalue ' + str(v) + ' is positive real. Swapping to last index!')
            #    val = swapArrayIndices(val, i, 3)
            #    vec = swapMatrixRows(vec, i, 3)
    t1 = val[2]
    t2 = val[3]
    if (isNumZero(t1.imag)) and (not isNumZero(t1.real)):
        #print('Third entry is real')
        if (t1.real > 0):
            #print('Third entry is positive real')
            val = swapArrayIndices(val, 2, 3)
            vec = swapMatrixRows(vec, 2, 3)
    return val, vec

# Classes for storing structure data
class Structure:

    # Sub class for defining each layer
    class Layer:
        # Instance variables for each Layer object
        def __init__(self, name, length, epsilon, mu):
            print('     Instanciating Layer')
            self.name = name
            self.length = length
            self.epsilon = epsilon
            self.mu = mu
            self.solution = np.zeros(4, dtype=complex)

        def __str__(self):
            try:
                return self.name + ': ' + str(self.length) + '\nEigen: ' + self.eigVal + self.eigVec
            except:
                return self.name + ': ' + str(self.length) 

    # Instance variables for each Structure object
    def __init__(self, num, omega, k1, k2):
        print('Instanciating Structure')
        self.num = num
        self.omega = omega
        self.c = 1
        self.kap = omega/self.c
        self.k1 = k1/self.kap
        self.k2 = k2/self.kap
        self.layers = []
        self.transferMatrices = []

    def printLayers(self):
        print('\nLAYERS: ')
        for layer in self.layers:
            print(layer)

    def removeLayer(self, n):
        print('Removing Layer')
        self.layers.pop(n)

    # Method for adding a layer to the structure
    def addLayer(self, name, length, epsilon, mu):
        print('Adding Layer')
        l = self.Layer(name, length, epsilon, mu)
        self.layers.append(l)

    # Method for inserting a layer into given index n
    def insertLayer(self, name, length, epsilon, mu, n):
        print('Inserting Layer')
        l = self.Layer(name, length, epsilon, mu)
        self.layers.insert(n, l)

    # Create the maxwell matrices
    def buildMatrices(self):
        print('\nBuilding Maxwell')
        maxwell_matrices = []
        for layer in self.layers:
            e = layer.epsilon
            u = layer.mu
            k1 = self.k1
            k2 = self.k2
            imaginary = complex(0, 1)
            omega = self.omega * imaginary
            m11 = omega * (-(e[2][0]*k1/e[2][2]) - (k2*u[1][2]/u[2][2]))
            m12 = omega * (-(e[2][1]*k1/e[2][2]) + (k1*u[1][2]/u[2][2]))
            m13 = omega * ((k1*k2/e[2][2]) + u[1][0] -
                        (u[1][2]*u[2][0]/u[2][2]))
            m14 = omega * ((-k1 ** 2/e[2][2]) + u[1]
                        [1] - (u[1][2]*u[2][1]/u[2][2]))
            m21 = omega * (-(e[2][0]*k2/e[2][2]) - (k2*(u[0][2]/u[2][2])))
            m22 = omega * (-(e[2][1]*k2/e[2][2]) - (k1*u[0][2]/u[2][2]))
            m23 = omega * ((k2 ** 2/e[2][2]) - u[0]
                        [0] + (u[0][2]*u[2][0]/u[2][2]))
            m24 = omega * (-(k1*k2/e[2][2]) - u[0]
                        [1] + (u[0][2]*u[2][1]/u[2][2]))
            m31 = omega * (-e[1][0] + (e[1][2]*e[2][0] /
                                    e[2][2]) - (k1*k2/u[2][2]))
            m32 = omega * (-e[1][1] + (e[1][2]*e[2][1] /
                                    e[2][2]) + (k1 ** 2/u[2][2]))
            m33 = omega * (-(e[1][2]*k2/e[2][2]) - (k1*u[2][0]/u[2][2]))
            m34 = omega * ((e[1][2]*k1/e[2][2]) - (k1*u[2][1]/u[2][2]))
            m41 = omega * (e[0][0] - (e[0][2]*e[2][0] /
                                    e[2][2]) - (k2 ** 2/u[2][2]))
            m42 = omega * (e[0][1] - (e[0][2]*e[2][1] /
                                    e[2][2]) + (k1*k2/u[2][2]))
            m43 = omega * ((e[0][2]*k2/e[2][2]) - (k2*u[2][0]/u[2][2]))
            m44 = omega * (-(e[0][2]*k1/e[2][2]) - (k2*u[2][1]/u[2][2]))
            maxwell_matrix = np.matrix([[m11,  m12, m13, m14],
                                        [m21,  m22, m23, m24],
                                        [m31,  m32, m33, m34],
                                        [m41,  m42, m43, m44]], dtype=complex)
            maxwell_matrices.append(maxwell_matrix)
        self.maxwell = maxwell_matrices
        return maxwell_matrices

    # Import maxwell matrices into struct
    def importMatrices(self, matrices):
        print('\nImporting Maxwell from Data')
        self.maxwell = matrices

    # Calculate the eigenproblem for all layers of the structure
    def calcEig(self):
        print('\nCalculating Eigen Problem')
        for n in range(len(self.layers)):
            print('For layer ' + str(n+1))
            start = time.perf_counter()
            eigVal, eigVec = np.linalg.eig(self.maxwell[n])
            eigVal, eigVec = organizeEigen(eigVal, eigVec)
            #eigVal, eigVec = organizeEigen(eigVal, eigVec)
            end = time.perf_counter()
            print(f'Time to calculate Eigensystem: {end-start:0.5f} seconds')
            self.layers[n].eigVal = eigVal
            self.layers[n].eigVec = np.transpose(eigVec)
            print(f'Values:\n{self.layers[n].eigVal}')
            print(f'Vectors:\n{self.layers[n].eigVec}')

    # Import preexisting eigendata
    def importEig(self, e_vals, e_vecs):
        print('\nImporting previously created eigendata')
        for n in range(len(self.layers)):
            self.layers[n].eigVal = e_vals[n]
            self.layers[n].eigVec = e_vecs[n]

    # Calculate the modes. Result will not contain constant multiplication
    def calcModes(self):
        print('\nCalculating Modes')
        for layer in self.layers:
            mode = np.zeros((1,4), dtype=complex)
            for n in range(4):
                mode += np.real(layer.eigVec[n]) * math.exp(np.real(layer.eigVal[n]) * layer.length)
            layer.modes = mode
            print(str(layer.name) + ' Modes:\nc*' + str(layer.modes))

    def printMaxwell(self):
        print('Maxwells:')
        for m in self.maxwell:
            print(m)

    def calcTransfer(self):
        interfaces = self.num-1
        transferMatrices = [] * interfaces
        if (self.num > 1):
            for i in range(interfaces):
                wNext = inv(np.transpose(self.layers[i+1].eigVec))
                w = np.transpose(self.layers[i].eigVec)
                if (i == 0):
                    zNorm = 0
                else:
                    zNorm = self.omega * self.layers[i].length
                expDia = expDiagonal(self.layers[i].eigVal, zNorm)
                tmp = w * expDia
                transferMatrices.append(wNext*tmp)
        self.transferMatrices = transferMatrices
        return transferMatrices

    def interfaces(self):
        interfaces = []
        interfaces.append(-self.layers[0].length)
        for i in range(self.num):
            interfaces.append(interfaces[i] + self.layers[i].length)
        return interfaces

    def calcScattering(self):
        layers = self.num
        I = layers - 1
        s = np.zeros((4*I,4*layers), dtype=complex)
        interfaces = np.zeros(layers)
        ifaces = self.interfaces()
        leftPsi = [np.zeros(4)] * I
        rightPsi = [np.zeros(4)] * I
        for i in range(layers):
            if i < I:
                interfaces[i] = 0
            else:
                interfaces[i] = ifaces[i] - ifaces[i-1]
        for i in range(I):
            expVecLeft = np.exp(np.multiply(self.layers[i].eigVal,(interfaces[i+1] - interfaces[i])))
            expVecRight = np.exp(self.layers[i+1].eigVal * (interfaces[i+1] - interfaces[i+1]))
            leftPsi[i] = np.matmul(np.transpose(self.layers[i].eigVec),np.diag(expVecLeft))
            rightPsi[i] = np.matmul(np.transpose(self.layers[i+1].eigVec),np.diag(expVecRight))
        for i in range(I):
            for j in range(4):
                for k in range(4):
                    s[4*i+j][4*i+k] = leftPsi[i].item(j,k)
                    s[4*i+j][4*(i+1)+k] = np.negative(rightPsi[i].item(j,k))
        self.scattering = s
        return s
        
    def calcConstants(self, c1, c2, c3, c4):
        scattering = self.calcScattering()
        layers = self.num
        interfaces = layers-1
        s = np.zeros((4*interfaces,4*interfaces), dtype=complex)
        f = np.zeros(4*interfaces, dtype=complex)
        for i in range(4*interfaces):
            for j in range(4*interfaces):
                s[i][j] = scattering[i][j+2]
        for i in range(4):
            f[i] = np.subtract(f.item(i), np.subtract(np.multiply(scattering[i][0],c1),np.multiply(scattering[i][1],c2)))
            aug = 4 * (interfaces-1) + i
            aug1 = 4 * (interfaces+1) - 1 
            aug2 = 4 * (interfaces+1) - 2 
            f[aug] = np.subtract(f.item(aug), np.subtract(np.multiply(scattering[aug][aug2],c3),np.multiply(scattering[aug][aug1],c4)))
        bPrime = solve(s,f)
        #bPrime = lstsq(s, f)[0] # May want to use solve instead
        b = np.zeros(4*layers, dtype=complex) # Constants vector
        b[0] = c1
        b[1] = c2
        b[4*layers-2] = c3
        b[4*layers-1] = c4
        for i in range(4*interfaces):
            b[i+2] = bPrime[i]
        self.constants = b
        for i in range(self.num):
            for j in range(4):
                # Set solution equal to the mode times the constant
                sol = self.layers[i].modes[0][j] * b[i*4+j]
                # Store solution in the structure layer
                self.layers[i].solution[j] = sol 
        return b
   
    # Get all the solutions for the structure
    def solution(self):
        solutions = []
        for n in range(self.num):
            solutions.append(self.layers[n].solution)
        return solutions

    # Check the sum of all the solutions is 0
    def checkSol(self):
        total = 0
        sol = self.solution()
        for n in range(self.num):
            for i in range(4):
                total =+ sol[n][i]
                
        return total == 0

    def printSol(self):
        print('Solutions:')
        for s in self.solution():
            print(s)

    # Determine the field of a structure
    def determineField(self, num_points=200):
        num_layers = self.num
        z_ends = self.interfaces()
        interfaces = [0] * num_layers
        z_arr = []
        Ex = []
        Ey = []
        Hx = []
        Hy = []

        for z in range(num_layers):
            if z < num_layers - 1:
                interfaces[z] = 0
            else:
                interfaces[z] = z_ends[z] - z_ends[z-1]

        for layer in range(num_layers):
            length = z_ends[layer+1] - z_ends[layer]

            current_c = self.constants[layer*4:4+(layer*4)]

            print("Constant vector at layer " + str(layer))
            print(current_c)

            for i in range(num_points):
                z = z_ends[layer] + (i * length) / num_points

                piScalar = np.pi * 0.4

                scalar = np.exp(np.multiply(np.complex(0.0, 1.0), piScalar))
                scalarMat = np.multiply(scalar, self.layers[layer].eigVec.transpose())
                expDiag = np.diag(np.exp(np.multiply(self.layers[layer].eigVal, (z - interfaces[layer]))))
                expMat = np.matmul(scalarMat, expDiag)
                fieldVec = np.matmul(expMat, current_c)

                if i == 0:
                    print("Field vec at " + str(z))
                    print(fieldVec)

                z_arr.append(z)
                Ex.append(fieldVec.item(0, 0).real)
                Ey.append(fieldVec.item(0, 1).real)
                Hx.append(fieldVec.item(0, 2).real)
                Hy.append(fieldVec.item(0, 3).real)

        field = {
            'z': z_arr,
            'Ex': Ex,
            'Ey': Ey,
            'Hx': Hx,
            'Hy': Hy
        }

        self.field = field
        return field

    # The structure string method
    def __str__(self):
        return 'Omega: ' + str(self.omega) + '\n(k1,k2): (' + str(self.k1*self.omega) + ',' + str(self.k2*self.omega) + ')\n'


# Test code
def test():
    print('Starting Test:')
    size = 3
    omega = 0.398
    k1 = 0.5
    k2 = 0.22
    s = Structure(size, omega, k1, k2)
    e = np.array([[1.5,0,0],
                [0,8,0],
                [0,0,1]])
    ee = np.array([[8,0,0],
                [0,1.5,0],
                [0,0,1]])
    u = np.array([[4,0,0],
                [0,1,0],
                [0,0,1]])
    uu = np.array([[1,0,0],
                [0,4,0],
                [0,0,1]])
    s.addLayer('Ambient Left', 10, e, u)
    s.addLayer('Layer 1', 7, e, u)
    s.addLayer('Ambient Right', 10, e, u)
    s.printLayers()
    s.removeLayer(1)
    s.printLayers()
    s.insertLayer('Layer 1', 7, ee, uu, 1)
    s.printLayers()
    m = s.buildMatrices()
    print('Number of Layers: ' + str(len(m)))
    print('Dimension of 1st layer Maxwell: ' + str(m[0].shape))
    print('Dimension of 2nd layer Maxwell: ' + str(m[1].shape))
    print('Dimension of 3rd layer Maxwell: ' + str(m[2].shape))
    print(s)
    s.printMaxwell()
    s.calcEig()
    s.calcModes()
    c1 = -1
    c2 = 0
    c3 = 0
    c4 = 0
    const = s.calcConstants(c1,c2,c3,c4)
    print('With incoming coefficients (' + str(c1)+ ', ' + str(c2)+ ') on the left and (' + str(c3)+ ', ' + str(c4)+ ') on the right')
    s.printSol()
    print('\nFinal Constants: \n' + str(const))
    print('\nFinal Scattering Matrix:\n' + str(s.scattering))
    print('Sum of solutions in all layers is 0: ' + str(s.checkSol()))
    print('\n\nEnd of test\n\n')
    
    
start = time.perf_counter()
test()
end = time.perf_counter()
print(f'Ran test in: {end-start:0.4f} seconds')
