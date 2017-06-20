# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 17:18:37 2017

@author: simon
"""

from pymira import spatialgraph
import os
import dill as pickle
import numpy as np

path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T'
f = os.path.join(path,'fix_graph_30_LRG_NET.am')
graph = spatialgraph.SpatialGraph()
graph.read(f)

flag = graph.get_field(name='flag')['data']

# Get node list
nodeFile = os.path.join(path,'nodeList.dill')
if not os.path.isfile(nodeFile):
    print 'Generating node list...'
    nodeList = graph.node_list()
            
    print('Pickling node list...')
    with open(nodeFile,'wb') as fo:
        pickle.dump(nodeList,fo)
else:
    with open(nodeFile ,'rb') as fo:
        nodeList = pickle.load(fo)

# Find nodes with branches with different flags
red_nodes = []
flagInd = [i for i,x in enumerate(nodeList[0].edges[0].scalarNames) if x=='flag']
assert len(flagInd)==1
flagInd = flagInd[0]

nred = 0
red_nodes = []
for node in nodeList:
    edges = node.edges
    curFlag = []
    for edge in edges:
        curFlag.append(edge.scalars[flagInd][0])
        
    curFlag = np.asarray(curFlag)
    if len(np.unique(curFlag))>1:
        node.add_scalar('red bond',1)
        nred += 1
        red_nodes.append(node)
    else:
        node.add_scalar('red bond',0)

# Load flowdata
path2= r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1'
graph2 = spatialgraph.SpatialGraph()
graph2.read(os.path.join(path2,'spatialGraph_RIN.am'))

import pdb
pdb.set_trace()
#graph.write(os.path.join(path,'red_bonds.am'))