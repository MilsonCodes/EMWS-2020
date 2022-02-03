#imports 
import numpy as np 
import time 
import scattering

arr = []
# Test code
def test():
    print('starting test')
    
    for i in range(1, 10, 1): 
        omega = i/10
        print("starting test for omega " , omega)
        size = 3
        k1 = 0.5
        k2 = 0.22
        global arr
        s = scattering.Structure(size, omega, k1, k2)
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
        c1 = 1
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
        arr.append(s)

def tester():
    start = time.perf_counter()
    test()
    end = time.perf_counter()
    print(f'Ran test in: {end-start:0.4f} seconds')

tester()
