# Converts an object loaded using AmiraMesh from https://github.com/CABI-SWS/reanimate
# There must be a folder on the python path called reanimate with the amiramesh.py file in it
from pymira import amiramesh as am
import json
from pathlib import Path
import os
join = os.path.join

def convert(filepath,opath=None):
    a = am.AmiraMesh()
    a.read(filepath,quiet=True)
    o = dict()
    # AmiraMesh object data held in fields by name:
    # 'VertexCoordinates', 'EdgeConnectivity', 'NumEdgePoints', 'EdgePointCoordinates', 'thickness'
    # Dislike the naming of thickness, so capitalize.
    # Data stored as numpy array, so call tolist()
    for field in a.fields:
        name = field['name']
        name = name[0].upper() + name[1:]
        if field['data'] is not None:
            o[name] = field['data'].tolist()
            
    if opath is not None:
        f = join(opath,Path(filepath).stem+'.json')
    else:
        f = filepath.replace('.am','.json')
    
    print(f)
    with open(f, 'w') as handle:
        json.dump(o, handle)
        
    return f
