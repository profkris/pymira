# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 09:37:47 2017

@author: simon

Merge two spatial grphs
Required for converting Paul Sweeney's files

"""

import pymira.spatialgraph as sp

def merge_graphs(graph1,graph2):

    dif1  = list(set(graph1.fieldNames) - set(graph2.fieldNames))
    dif2  = list(set(graph2.fieldNames) - set(graph1.fieldNames))
    
    for fName in dif2:
        f = graph2.get_field(fName)
        marker = graph1.generate_next_marker()
        f['marker'] = marker
        print('Adding {} {}...'.format(marker,fName))
        graph1.fields.append(f)

#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\'  
#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\'
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
fp = dir_+'Press2Amira.txt'
ff = dir_+'Flow2Amira.txt'
fs = dir_+'Stress2Amira.txt'
fv = dir_+'Velocity2Amira.txt'
mFiles = [fp,ff,fs]#,fv]

fo = dir_+'spatialGraph.am'

graph = sp.SpatialGraph()
print('Reading source graph: {}'.format(mFiles[0]))
graph.read(mFiles[0])

for f in mFiles[1:]:
    graph_to_add = sp.SpatialGraph()
    graph_to_add.read(f)
    merge_graphs(graph,graph_to_add)

graph.write(fo)