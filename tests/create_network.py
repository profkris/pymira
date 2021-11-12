# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 07:18:10 2017

@author: simon
"""

import numpy as np
import spatialgraph

nodes = [[0,0,0],
         [1,0,0],
         [0.5,0,0],
         [0,1,0],
         [2,2,0],
         [5,5,0],
         [2,5,0],
        ]
        
connectivity = [[0,1],
                [2,3],
                [1,4],
                [4,5],
                [4,6],
                [0,0]
                ]
        
points = [[0,0,0],
         [0.5,0.,0.],
         [1,0,0],
         [0,0,0],
         [0,0.5,0],
         [0,1,0],
         [1,0,0],
         [1.5,0,0],
         [2,2,0],
         [2,2,0],
         [3,3,0],
         [4,4,0],
         [5,5,0],
         [2,2,0],
         [2,3,0],
         [2,4,0],
         [2,5,0],
         [0,0,0],
         [1.5,0,0],
         [1.5,1.5,0],
         [0,1.5,0],
         [0,0,0]
        ]
        
radii = [3,
          3,
          3,
          1,
          1,
          1,
          2,
          2,
          2,
          1,
          1,
          1,
          1,
          1.5,
          1.5,
          1.5,
          1.5,
          0.5,
          0.5,
          0.5,
          0.5,
          0.5,
          ]
        
nedgepoints = [3,
               3,
               3,
               4,
               4,
               5
               ]

nodes = np.asarray(nodes,dtype='float')
nedgepoints = np.asarray(nedgepoints,dtype='int')
connectivity = np.asarray(connectivity,dtype='int')
edgepoints = np.asarray(points,dtype='float')
radii = np.asarray(radii,dtype='float')

nnode = nodes.shape[0]
nedge = nedgepoints.shape[0]
npoint = edgepoints.shape[0]
    
graph = spatialgraph.SpatialGraph(initialise=True,scalars=['Radii'])
graph.set_definition_size('VERTEX',nnode)
graph.set_definition_size('EDGE',nedge)
graph.set_definition_size('POINT',npoint)
graph.set_data(nodes,name='VertexCoordinates')
graph.set_data(connectivity,name='EdgeConnectivity')
graph.set_data(nedgepoints,name='NumEdgePoints')
graph.set_data(edgepoints,name='EdgePointCoordinates')
graph.set_data(radii,name='Radii')

#graph.sanity_check(deep=True)

#import pdb
#pdb.set_trace()
#graph.constrain_nodes(xrange=[2,5])
graph.delete_nodes([2,5])

graph.sanity_check(deep=True)

#graph._print()
ofile = 'C:\\Anaconda2\\Lib\\site-packages\\pymira\\test_network.am'
graph.write(ofile)
