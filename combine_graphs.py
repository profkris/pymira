# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 09:37:47 2017

@author: simon

Merge two spatial grphs
Required for converting Paul Sweeney's files

"""

import pymira.spatialgraph as sp
import os
import numpy as np

def merge_graphs(graph1,graph2):

    dif1  = list(set(graph1.fieldNames) - set(graph2.fieldNames))
    dif2  = list(set(graph2.fieldNames) - set(graph1.fieldNames))
    
    for fName in dif2:
        f = graph2.get_field(fName)
        marker = graph1.generate_next_marker()
        f['marker'] = marker
        print('Adding {} {}...'.format(marker,fName))
        graph1.fields.append(f)

for i in np.arange(2,13,1):
    #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\'  
    dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\{}\\'.format(i)
    #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
    #dir_ = r'C:\Users\simon\Dropbox\170606_Ben Vessel Networks\C1M3\2%'
    #dir_ = r'C:\Users\simon\Dropbox\Glioma'
    fn = os.path.join(dir_,'Press2Amira.txt')
    #fp = os.path.join(dir_,'Press2Amira.txt')
    ff = os.path.join(dir_,'Flow2Amira.txt')
    fs = os.path.join(dir_,'Stress2Amira.txt')
    fv = os.path.join(dir_,'Velocity2Amira.txt')
    #mFiles = [fp,ff,fs]#,fv]
    mFiles = [fn,ff,fv,fs]
    
    fo = os.path.join(dir_,'spatialGraph.am')
    
    graph = sp.SpatialGraph()
    print('Reading source graph: {}'.format(mFiles[0]))
    graph.read(mFiles[0])
    
    for f in mFiles[1:]:
        graph_to_add = sp.SpatialGraph()
        graph_to_add.read(f)
        merge_graphs(graph,graph_to_add)
    
    graph.write(fo)