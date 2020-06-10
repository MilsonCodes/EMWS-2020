from flask import Flask
from flask import request
from flask import json
import sys, os
sys.path.append('..')
from backend.scripts.scattering import Structure as s
import numpy

# Run server by calling python app.py
app = Flask(__name__)

# Base route example
@app.route('/')
def hello_world():
    return 'Hello, World!'

# Function to encode complex numbers into tuples to allow for json serialization
def encode_complex(z):
    try:
        return (z.real, z.imag)
    except:
        return (z, 0)

def decode_complex((x,y)):
    z
    try:
        z.real = x
        z.imag = y
        return z
    except:
        z.real = x
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

def encode_eigen(m):
    n = [0, 0, 0, 0]
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
        m = encode_maxwell(struct.maxwell[i])
        n = encode_eigen(layer.eigVal.tolist())
        o = encode_eigen(layer.eigVal.tolist())


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

if __name__ == '__main__':
    app.run()
