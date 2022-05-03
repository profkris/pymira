# Converts an object loaded using AmiraMesh from https://github.com/CABI-SWS/reanimate
# There must be a folder on the python path called reanimate with the amiramesh.py file in it
from pymira import amiramesh as am
import json

def convert(filepath):
    a = am.AmiraMesh()
    a.read(filepath)
    o = dict()
    # AmiraMesh object data held in fields by name:
    # 'VertexCoordinates', 'EdgeConnectivity', 'NumEdgePoints', 'EdgePointCoordinates', 'thickness'
    # Dislike the naming of thickness, so capitalize.
    # Data stored as numpy array, so call tolist()
    for field in a.fields:
        name = field['name']
        name = name[0].upper() + name[1:]
        o[name] = field['data'].tolist()
    f = open(filepath + '.json', 'w')
    json.dump(o, f)
    f.close()
