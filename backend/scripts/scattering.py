import math
import time
import numpy as np
from scipy.linalg import null_space, solve

# Set precision for printing arrays
np.set_printoptions(precision=6, suppress=True)

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
            self.solution = np.zeros(4)

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

    def insertLayer(self, name, length, epsilon, mu, n):
        print('Inserting Layer')
        l = self.Layer(name, length, epsilon, mu)
        self.layers.insert(n, l)

    # Create the maxwell matrices
    def buildMatrices(self):
        print('\nBuilding Maxwell')
        maxwell_matrices = []
        for layer in self.layers:
            # layer = layers[num]
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
                                    [m41,  m42, m43, m44]])
            # maxwell_matrix = map(lambda x: x * omega, maxwell_matrix)
            # maxwell_matrix = map(lambda x: map(
            #     lambda y: complex(y, 1), x), maxwell_matrix)
            maxwell_matrices.append(maxwell_matrix)
            # m = list(maxwell_matrix)
            # print("layer" + str(num+1) + ": ")
            # print("\n")
        self.maxwell = maxwell_matrices
        return maxwell_matrices

    def calcEig(self):
        print('\nCalculating Eigen Problem')
        for n in range(len(self.layers)):
            print('For layer ' + str(n+1))
            start = time.perf_counter()
            eig = np.linalg.eig(self.maxwell[n])
            end = time.perf_counter()
            print(f'Time to calculate Eigensystem: {end-start:0.5f} seconds')
            self.layers[n].eigVal = eig[0]
            self.layers[n].eigVec = eig[1]
            print(f'Values:\n{self.layers[n].eigVal}')
            print(f'Vectors:\n{self.layers[n].eigVec}')

    def calcModes(self):
        print('\nCalculating Modes')
        for layer in self.layers:
            mode = np.zeros((1,4))
            for n in range(4):
                mode += np.real(layer.eigVec[n]) * math.exp(np.real(layer.eigVal[n]) * layer.length)
            layer.modes = mode
            print(str(layer.name) + ' Modes:\nc*' + str(layer.modes))

    def printMaxwell(self):
        print('Maxwells:')
        for m in self.maxwell:
            print(m)

    def __str__(self):
        return 'Omega: ' + str(self.omega) + '\n(k1,k2): (' + str(self.k1*self.omega) + ',' + str(self.k2*self.omega) + ')\n'
    def calcConstants(self, c1, c2, c3, c4):
        c = np.zeros((self.num*4,1))
        c[0][0] = c1
        c[1][0] = c2
        c[self.num*4-2][0] = c3
        c[self.num*4-1][0] = c4
        a = np.zeros((self.num*4,self.num*4))
        f = c
        for n in range(self.num):
            aug = 4 * n
            for i in range(4):
                for j in range(4):
                    # Set the diagonals of A equal to the maxwell for each layer
                    a[i+aug][j+aug] = self.maxwell[n].A[i][j].imag
                    #if n != self.num-1:
                        # Set the off diagonal entries equal to the complex conjugate
                    #    a[i+aug][j+4+aug] = -self.maxwell[n+1].real.getH().item((i,j))
        for n in range(4):
            f[n][0] = f[n][0] - ((a[n][0]*c1) - (a[n][1]*c2))
            m = 4 * (self.num-1) + n
            i = 4 * (self.num-1) - 2
            j = 4 * (self.num-1) - 1
            f[m][0] = f[m][0] - ((a[m][i]*c3) - (a[m][j]*c4))
        x = solve(a,f)
        for i in range(self.num):
            for j in range(4):
                sol = self.layers[i].modes[0][j] * f[i*4+j][0]
                self.layers[i].solution[j] = sol 
        return x
   
    def solution(self):
        solutions = []
        for n in range(self.num):
            solutions.append(self.layers[n].solution)
        return solutions

    def checkSol(self):
        total = 0
        sol = self.solution()
        for n in range(self.num):
            for i in range(4):
                total =+ sol[n][i]
                
        return total


    def printSol(self):
        print('Solutions:')
        for s in self.solution():
            print(s)

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
    print('Final Constants: \n' + str(s.calcConstants(c1,c2,c3,c4)))
    print('With incoming coefficients (' + str(c1)+ ', ' + str(c2)+ ') on the left and (' + str(c3)+ ', ' + str(c4)+ ') on the right')
    s.printSol()
    print('Sum of solutions in all layers: ' + str(s.checkSol()))
    print('\nEnd of test\n\n')
    
    
start = time.perf_counter()
test()
end = time.perf_counter()
print(f'Ran test in: {end-start:0.4f} seconds')
