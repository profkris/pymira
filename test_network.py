# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:26:51 2017

@author: simon
"""

import numpy as np
from pymira import spatialgraph

nodes = [[0.,0.,0.],
         [1.,0.,0.],
         [2.,0.,0.],
         [2.,1.,0.],
         [2.,1.,1.]]
nodes = np.asarray(nodes,dtype='float')

edgeConn = [[0,1],
            [1,2],
            [2,3],
            [2,4]]
edgeConn = np.asarray(edgeConn,dtype='int')            
            
edgePoints = [nodes[0],
              [0.1,0.,0.],
              [0.3,0.,0.],
              [0.5,0.1,0.],
              [0.75,0.05,0.],
              nodes[1],

              nodes[1],
              nodes[2],

              nodes[2],
              nodes[3],

              nodes[2],
              nodes[4] 
              ]
edgePoints = np.asarray(edgePoints,dtype='float')
              
radii = np.zeros(edgePoints.shape[0]) + 0.05
              
nedgepoints = [6,
               2,
               2,
               2]
nedgepoints = np.asarray(nedgepoints,dtype='int')               
               
graph = spatialgraph.SpatialGraph(initialise=True,scalars=['Radii'])
graph.set_definition_size('VERTEX',nodes.shape[0])
graph.set_definition_size('EDGE',edgeConn.shape[0])
graph.set_definition_size('POINT',edgePoints.shape[0])
graph.set_data(nodes,name='VertexCoordinates')
graph.set_data(edgeConn,name='EdgeConnectivity')
graph.set_data(nedgepoints,name='NumEdgePoints')
graph.set_data(edgePoints,name='EdgePointCoordinates')
graph.set_data(radii,name='Radii')
ofile = 'C:\\Anaconda2\\Lib\\site-packages\\pymira\\test_graph.am'
graph.write(ofile)