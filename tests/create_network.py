# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 07:18:10 2017

@author: simon
"""

import numpy as np
arr = np.asarray
from pymira import spatialgraph

def bifurcation():

    nodes = [[0,0,0],
             [100,0,0],
             [200.,100,0],
             [200.,-100.,0],
            ]
            
    connectivity = [[0,1],
                    [1,2],
                    [1,3],
                    ]
                    
                    
    vessel_type = [0,
                   0,
                   0,
                   ]                    
                    
    edgeradii = [ 30.,
                  20.,
                  20.,
                 ] 
                
    return nodes, connectivity, edgeradii, vessel_type
    
def double_bifurcation():

    nodes = [[0,0,0],
             [100,0,0],
             [200.,100,0],
             [200.,-100.,0],
             [300.,200,0],
             [300.,100.,0],
             [300.,-200,0],
             [300.,-100.,0],
            ]
            
    connectivity = [[0,1],
                    [1,2],
                    [1,3],
                    [2,4],
                    [2,5],
                    [3,6],
                    [3,7],
                    ]
                    
                    
    vessel_type = [0,
                   0,
                   0,
                   
                   2,
                   2,
                   2,
                   2,
                   ]                    
                    
    edgeradii = [ 30.,
                  20.,
                  20.,
                  10.,
                  10.,
                  10.,
                  10.,
                 ] 
                
    return nodes, connectivity, edgeradii, vessel_type
    
def double_bifurcation_reconnected():

    nodes = [[0,0,0],
             [100,0,0],
             [200.,100,0],
             [200.,-100.,0],
             
             [300.,200,0],
             [300.,100.,0],
             [300.,-200,0],
             [300.,-100.,0],
             
             [400.,100,0],
             [400.,-100.,0],
             [500,0,0],
             [600,0,0],
            ]
            
    connectivity = [[0,1],
                    [1,2],
                    [1,3],
                    
                    [2,4],
                    [2,5],
                    [3,6],
                    [3,7],
                    
                    [4,8],
                    [5,8],
                    [6,9],
                    [7,9],
                    
                    [8,10],
                    [9,10],
                    [10,11],
                    ]
                    
    vessel_type = [0,
                   0,
                   0,
                   
                   2,
                   2,
                   2,
                   2,
                   
                   2,
                   2,
                   2,
                   2,
                   
                   1,
                   1,
                   1,
                   ]
                   
                    
    edgeradii = [ 30.,
                  20.,
                  20.,
                  
                  10.,
                  10.,
                  10.,
                  10.,
                  
                  10.,
                  10.,
                  10.,
                  10.,
                  
                  20.,
                  20.,
                  30.,
                 ] 
                
    return nodes, connectivity, edgeradii, vessel_type

def custom_graph():

    nodes = [[0,0,0],
             [100,0,0],
             [50.,0,0],
             [0,100.,0],
             [200,200,0],
             [500,500,0],
             [200,500,0],
            ]
            
    connectivity = [[0,1],
                    [2,3],
                    [1,4],
                    [4,5],
                    [4,6],
                    [0,0]
                    ]
                    
    vessel_type = np.zeros(len(connectivity))
                    
    edgeradii = [ 30.,
                  10.,
                  20.,
                  10.,
                  15.,
                  5.,
                 ] 
                 
    return nodes, connectivity, edgeradii, vessel_type
    
nodes, connectivity, edgeradii, vessel_type = double_bifurcation_reconnected()
                
edgepoints,nedgepoints,radii,category = [],[],[],[]
for i,conn in enumerate(connectivity):
    edgepoints.append(nodes[conn[0]])
    edgepoints.append(nodes[conn[1]])
    nedgepoints.append(2)
    radii.append([edgeradii[i],edgeradii[i]])
    category.append([vessel_type[i],vessel_type[i]])
edgepoints = arr(edgepoints)
nedgepoints = arr(nedgepoints)
radii = arr(radii).flatten()
category = arr(category).flatten()

nodes = np.asarray(nodes,dtype='float')
nedgepoints = np.asarray(nedgepoints,dtype='int')
connectivity = np.asarray(connectivity,dtype='int')
edgepoints = np.asarray(edgepoints,dtype='float')
radii = np.asarray(radii,dtype='float')
category = np.asarray(category,dtype='int')

nnode = nodes.shape[0]
nedge = nedgepoints.shape[0]
npoint = edgepoints.shape[0]
    
graph = spatialgraph.SpatialGraph(initialise=True,scalars=['Radii','VesselType'])
graph.set_definition_size('VERTEX',nnode)
graph.set_definition_size('EDGE',nedge)
graph.set_definition_size('POINT',npoint)
graph.set_data(nodes,name='VertexCoordinates')
graph.set_data(connectivity,name='EdgeConnectivity')
graph.set_data(nedgepoints,name='NumEdgePoints')
graph.set_data(edgepoints,name='EdgePointCoordinates')
graph.set_data(radii,name='Radii')
graph.set_data(category,name='VesselType')

graph.sanity_check(deep=True)

#graph._print()
ofile = '/mnt/data2/retinasim/data/cco_circ_domain/graph/test_network.am'
graph.plot_graph()
graph.write(ofile)
