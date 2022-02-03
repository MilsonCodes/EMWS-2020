from flask import Flask, render_template
from flask import request
from flask import json
import numpy as np 
from flask_cors import CORS, cross_origin
import sys
sys.path.append('/Users/joelkeller/EMWS-2020/backend')
from scattering import Structure as s

# Run server by calling python app.py
app = Flask(__name__)
# List of accepted origins
origins = ["http://localhost:8000", "https://www.math.lsu.edu"]
# Change origins to '*' if this solution gives issues
CORS(app, resources={r"/structure": {"origins": origins}})

# Greeting message, currently unused
base =  '''
        Welcome to the EMWS API!
        Source code and documentation can be found here:
        \n\thttps://github.com/MilsonCodes/EMWS-2020
        \n\n
        The live site can be found here:
        \n\thttps://www.math.lsu.edu/~shipman/EMWS/html/dashboard.4.html
 
 
    '''

# Function to encode complex numbers into tuples to allow for json serialization
def encode_complex(z):
    try:
        return { 're': z.real, 'im': z.imag }
    except:
        return z

def decode_complex(val):
    z = None
    try:
        z = complex(val['re'], val['im'])
        return z
    except:
        z = val
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
    n = np.zeros((4,4), dtype=complex)
    for i in range(4):
        for j in range(4):
            n[i][j] = decode_complex(m[i][j])
    return n


# Function to replace all elements of a 4 item vector
def encode_eigen(m):
    n = [0, 0, 0, 0]
    for i in range(len(m)):
        n[i] = encode_complex(m[i])
    return n

def decode_eigen(m):
    n = np.zeros(4, dtype=complex)
    for i in range(4):
        n[i] = decode_complex(m[i])
    return n

# Function to replace all elements of each 4 item vector
def encode_evecs(m):
    n = [0, 0, 0, 0]
    for i in range(len(m)):
        n[i] = encode_eigen(m[i])
    return n

def decode_evecs(m):
    n = np.zeros((4,4), dtype=complex)
    for i in range(4):
        for j in range(4):
            n[i][j] = decode_complex(m[i][j])
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
    n = np.zeros(len(m), dtype=complex)

    for i in range(len(m)):
        n[i] = decode_complex(m[i])

    return n

# Route for creating a crystal structure and calculating eigen problem
@app.route('/structure/omega', methods=['POST'])
@cross_origin()
def omegas():
    assert request.method == 'POST'
    # print('Mode route beginning')
    # Parse data
    req = json.loads(request.data)
    initial_omega = req['initial_omega']
    final_omega = req['final_omega']
    k1 = req['k1']
    k2 = req['k2']
    layers = req['layers']
    step = req['step']
    num = len(layers)

    #for loop omega 
    for omega in range(initial_omega, final_omega, step):
    # Create structure
        struct = s(num, omega, k1, k2)
        for layer in layers:
            struct.addLayer(layer['name'], layer['length'], layer['epsilon'], layer['mu'])
        # Calculate and build structure
        struct.buildMatrices()
        struct.calcEig()
        struct.calcModes()

    # Create list of values for response
    maxwells = []
    e_vals = []
    e_vecs = []
    modes = []
    i = 0
    for layer in struct.layers:
        m = encode_maxwell(struct.maxwell[i])
        n = encode_eigen(layer.eigVal.tolist())
        o = encode_evecs(layer.eigVec.tolist())
        mm = encode_evecs(layer.modes.tolist())


        maxwells.append(m)
        e_vals.append(n)
        e_vecs.append(o)
        modes.append(mm)
        i += 1

    # Prepare response data
    data = {
        'maxwell_matrices': maxwells,
        'eigenvalues': e_vals,
        'eigenvectors': e_vecs,
        'modes': modes
    }
    
    # print('Mode route completed')

    return json.jsonify(data)


