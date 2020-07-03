from flask import Flask
from flask import request
from flask import json
import sys, os
sys.path.append('..')
from backend.scripts.scattering import Structure as s
import numpy

# Run server by calling python app.py
app = Flask(__name__)

base = '''
Welcome to the EMWS API!

Source code and documentation can be found here:
    https://github.com/MilsonCodes/EMWS-2020

The live site can be found here:
    https://www.math.lsu.edu/~shipman/EMWS/html/dashboard.4.html
'''

# Base route example
@app.route('/')
def hello_world():
    return base

# Function to encode complex numbers into tuples to allow for json serialization
def encode_complex(z):
    try:
        return { 're': z.real, 'im': z.imag }
    except:
        return z

def encode_vector(v):
    vec = []
    for n in range(len(v)):
        vec.append(encode_complex(v[n]))
    return vec

# Encode 3d array
def encode_matrix(m):
    size = m.size
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                m[i][j][k] = encode_complex(m[i][j][k])
    return m
                

def decode_complex(val):
    z
    try:
        z.real = val.re
        z.imag = val.im
        return z
    except:
        z.real = val
        return z

# Function to replace all elements of a 4x4 array with tuples
def encode_maxwell(m):
    n = [[0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0]]
    for i in range(4):
        for j in range(4):
            n[i][j] = encode_complex(m.item(i,j))
    return n

def decode_maxwell(m):
    n = [[0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0]]
    for i in range(4):
        for j in range(4):
            n[i][j] = decode_complex(m.item(i,j))
    return n


# Function to replace all elements of a 4 item vector
def encode_eigen(m):
    n = [0, 0, 0, 0]
    for i in range(len(m)):
        n[i] = encode_complex(m[i])
    return n

def decode_eigen(m):
    n = [0, 0, 0, 0]
    for i in range(len(m)):
        n[i] = decode_complex(m[i])
    return n

# Function to replace all elements of each 4 item vector
def encode_evecs(m):
    n = [0, 0, 0, 0]
    for i in range(len(m)):
        n[i] = encode_eigen(m[i])
    return n

def decode_evecs(m):
    n = [0, 0, 0, 0]
    for i in range(len(m)):
        n[i] = decode_eigen(m[i])
    return n

# Function set to encode/decode a scattering matrix
def encode_scattering(m):
    n = [[0 for i in range(len(m[j]))] for j in range(len(m))]

    for i in range(len(m)):
        for j in range(len(m[i])):
            n[i][j] = encode_complex(m[i][j])

    return n

def decode_scattering(m):
    n = [[0 for i in range(len(m[j]))] for j in range(len(m))]

    for i in range(len(m)):
        for j in range(len(m[i])):
            n[i][j] = decode_complex(m[i][j])

    return n

# Function set to encode/decode a constant vector
def encode_constants(m):
    n = [0] * len(m)

    for i in range(len(m)):
        n[i] = encode_complex(m[i])
    
    return n

def decode_constants(m):
    n = [0] * len(m)

    for i in range(len(m)):
        n[i] = decode_complex(m[i])
    
    return n

# Route for creating a crystal structure and calculating eigen problem
@app.route('/structure/modes', methods=['POST'])
def modes():
    assert request.method == 'POST'
    # Parse data
    req = json.loads(request.data)
    omega = req['omega']
    k1 = req['k1']
    k2 = req['k2']
    layers = req['layers']
    num = len(layers)
    struct = s(num, omega, k1, k2)
    for layer in layers:
        struct.addLayer(layer['name'], layer['length'], layer['epsilon'], layer['mu'])
    struct.buildMatrices()
    struct.calcEig()
    maxwells = []
    e_vals = []
    e_vecs = []
    i = 0
    for layer in struct.layers:
        print(layer.eigVec)
        m = encode_maxwell(struct.maxwell[i])
        n = encode_eigen(layer.eigVal.tolist())
        o = encode_evecs(layer.eigVec.tolist())


        maxwells.append(m)
        e_vals.append(n)
        e_vecs.append(o)
        i += 1

    data = {
        'maxwell_matrices': maxwells,
        'eigenvalues': e_vals,
        'eigenvectors': e_vecs
    }

    return json.jsonify(data)

# Route for getting data points
@app.route('/structure/field', methods=['POST'])
def field():
    assert request.method == 'POST'
    req = json.loads(request.data)
    omega = req['omega']
    k1 = req['k1']
    k2 = req['k2']
    layers = req['layers']
    num = len(layers)
    struct = s(num, omega, k1, k2)
    for layer in layers:
        struct.addLayer(layer['name'], layer['length'], layer['epsilon'], layer['mu'])
    
    # Setup return data
    data = {}

    # Get existing data
    maxwell_matrices = None
    eigenvalues = None
    eigenvectors = None

    try:
        maxwell_matrices = req['maxwell_matrices']
        eigenvalues = req['eigenvalues']
        eigenvectors = req['eigenvectors']
    except Exception:
        print('\nFailed to maxwell or eigendata. Will calculate data')
    
    # Handle maxwells
    if maxwell_matrices == None:
        struct.buildMatrices()
        maxwells = []
        for maxwell in struct.maxwell:
            m = encode_maxwell(maxwell)
            maxwells.append(m)
        data['maxwell_matrices'] = maxwells
    else:
        maxwells = []
        for maxwell in maxwell_matrices:
            maxwells.append(decode_maxwell(maxwell))
        struct.importMatrices(maxwells)

    #Handle eigendata
    if eigenvalues == None or eigenvectors == None:
        struct.calcEig()
        struct.calcModes()
        e_vals = []
        e_vecs = []
        i = 0
        for layer in struct.layers:
            n = encode_eigen(layer.eigVal.tolist())
            o = encode_evecs(layer.eigVec.tolist())

            e_vals.append(n)
            e_vecs.append(o)
            i += 1
        data['eigenvalues'] = e_vals
        data['eigenvectors'] = e_vecs
    else:
        e_vals = []
        e_vecs = []
        for vals in eigenvalues:
            e_vals.append(decode_eigen(vals))
        for vecs in eigenvectors:
            e_vecs.append(decode_evecs(vecs))
        struct.importEig(e_vals, e_vecs)
        struct.calcModes()

    # Calculate Scattering Matrix and Constants
    incoming = None
    try:
        incoming = decode_eigen(req['incoming'])
    except Exception:
        print('\nDid not find incoming constants! Using defaults...')
        incoming = [1, 0, 0, 0]

    struct.calcScattering()
    struct.calcConstants(incoming[0], incoming[1], incoming[2], incoming[3])
    data['scattering'] = encode_scattering(struct.scattering)
    data['constants'] = encode_constants(struct.constants)

    num_points = None
    try:
        num_points = req['num_points']
    except Exception:
        num_points = 200

    field = struct.determineField(num_points)
    data['field'] = field

    return json.jsonify(data)


@app.route('/structure/constants', methods=['POST'])
def constants():
    assert request.method == 'POST'
    # Parse data
    req = json.loads(request.data)
    omega = req['omega']
    k1 = req['k1']
    k2 = req['k2']
    layers = req['layers']
    try:
        c = req['incoming']
    except:
        print('No incoming coeffecients found, using defaults.')
        c = [1, 0, 0, 0]
    maxwell = False
    eigen = False
    e_vals = []
    e_vecs = []
    num = len(layers)
    struct = s(num, omega, k1, k2)
    try:    
        maxwells = decode_maxwell(req['maxwell'])
        maxwell = True
    except:
        maxwells = []
        print('No maxwell matrix included')
    for layer in layers:
        struct.addLayer(layer['name'], layer['length'], layer['epsilon'], layer['mu'])
        try:
            e_vals.append(decode_eigen(layer['eig_values']))
            e_vecs.append(decode_evecs(layer['eig_vectors']))
            eigen = True
        except:
            print('No eigen data included')
    if (maxwell == False):
        struct.buildMatrices()
    else:
        struct.maxwell = maxwells[n]
    if (eigen == False):
        struct.calcEig()
    else:
        struct.importEig(e_vals, e_vecs)
    maxwells = []
    e_vals = []
    e_vecs = []
    i = 0
    for layer in struct.layers:
        print(layer.eigVec)
        m = encode_maxwell(struct.maxwell[i])
        n = encode_eigen(layer.eigVal)
        o = encode_evecs(layer.eigVec.tolist())

        maxwells.append(m)
        e_vals.append(n)
        e_vecs.append(o)
        i += 1

    struct.calcModes()
    const = struct.calcConstants(c[0], c[1], c[2], c[3])
    constants = encode_vector(const)

    scattering_matrix = []
    for n in range(len(struct.scattering)):
        scat = encode_vector(struct.scattering[n])
        scattering_matrix.append(scat)

    data = {
        'maxwell_matrices': maxwells,
        'eigenvalues': e_vals,
        'eigenvectors': e_vecs,
        'scattering': scattering_matrix,
        'constants': constants
    }

    return json.jsonify(data)

if __name__ == '__main__':
    app.run()
