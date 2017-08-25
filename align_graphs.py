# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 07:45:36 2017

@author: simon
"""

path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T'
f_flag = os.path.join(path,'spatialGraph_flag_RIN_scaled.am')
graph_flag = spatialgraph.SpatialGraph()
graph_flag.read(f_flag)

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
        
# Load flowdata
path2= r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1'
graph_flow = spatialgraph.SpatialGraph()
graph_flow.read(os.path.join(path2,'spatialGraph_RIN.am'))

# Get node list
nodeFile = os.path.join(path2,'nodeList.dill')
if not os.path.isfile(nodeFile):
    print 'Generating node list...'
    nodeList_flow = graph.node_list()
            
    print('Pickling node list...')
    with open(nodeFile,'wb') as fo:
        pickle.dump(nodeList_flow,fo)
else:
    with open(nodeFile ,'rb') as fo:
        nodeList_flow = pickle.load(fo)
        
#import pdb
#pdb.set_trace()