from flask import Flask
from flask import request
import sys, os
sys.path.append(os.path.abspath('../scripts'))
from scattering import Structure as s

# Run server by calling python app.py
app = Flask(__name__)

@app.route('/')
def hello_world():
    return 'Hello, World!'

@app.route('/structure/modes', methods=['POST'])
def modes():
    assert request.method == 'POST'
    print(request.form)
    omega = request.json['omega']
    k1 = request.json['k1']
    k2 = request.json['k2']
    layers = request.json['layers']
    num = len(layers)
    struct = s(num, omega, k1, k2)
    for layer in layers:
        struct.addLayer(layer['name'], layer['length'], layer['epsilon'], layer['mu'])
    maxwells = struct.buildMatrices()
    struct.calcEig()
    e_vals = []
    e_vecs = []
    for layer in struct.layers:
        e_vals.append(layer.eigVal)
        e_vecs.append(layer.eigVec)

    data = {
        'maxwell_matrices': maxwells,
        'eigenvalues': e_vals,
        'eigenvectors': e_vecs
    }

    return data

if __name__ == '__main__':
    app.run()
