# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 14:50:11 2018

@author: simon
"""

from pymira import spatialgraph
import os
import copy
import numpy as np

#path = r'C:\Users\simon\Dropbox'
#file = os.path.join(path,'ct_output_exponential.am')
#opath = os.path.join(path,'ct_output_stack')

path = r'C:\Users\simon\Dropbox\Tumour'
file = os.path.join(path,'ct_output_parker_half_tumour.am')
opath = os.path.join(path,'ct_output_stack_parker_half')

graph = spatialgraph.SpatialGraph()
print('Reading graph: {}'.format(file))
graph.read(file)
print('Graph read')

tinds = [i for i,n in enumerate(sg.fieldNames) if 'TP' in n]
inds = [i for i,n in enumerate(sg.fieldNames) if 'TP' not in n]

times = [int(i.replace('TP','')) for i in graph.fieldNames if 'TP' in i]

for t in times:
    print('Copying time point {}'.format(t))
    gc = copy.deepcopy(graph)
    fields = [gc.fields[j] for j,n in enumerate(gc.fieldNames) if n=='TP{}'.format(t) or 'TP' not in n]
    fieldNames = [f['name'] if 'TP' not in f['name'] else 'Concentration' for f in fields]
    for i,n in enumerate(fieldNames):
        fields[i]['name'] = n
        
    gc.fields = fields
    gc.fieldNames = fieldNames
    
    ofile = os.path.join(opath,'concentration_t{}.am'.format(t))
    print(('Writing {}'.format(ofile)))
    #gc.add_parameter('Time',np.asscalar(t))
    gc.parameters[0]['Time'] = t #np.asscalar(t)
    try:
        gc.write(ofile)
    except Exception as e:
        print(e)