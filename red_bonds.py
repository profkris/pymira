# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 17:18:37 2017

@author: simon
"""

from pymira import spatialgraph
import os
import dill as pickle
import numpy as np

#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T'
path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T'
f = os.path.join(path,'spatialGraph_flag_RIN_scaled.am')
#f = os.path.join(path,'spatialGraph_RIN.am')
graph = spatialgraph.SpatialGraph()
graph.read(f)

#path_post = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T - Post-VDA\1'
#fpost = os.path.join(path,'spatialGraph_RIN.am')
#graphPost = spatialgraph.SpatialGraph()
#graphPost.read(f)

flag = graph.get_field(name='flag')['data']

# Get node list
nodeFile = os.path.join(path,'nodeList.dill')
#if not os.path.isfile(nodeFile):
if True:
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
        
#trMatParam = [x for x in graph.parameters if x['parameter']=='TransformationMatrix']
#trMat = trMatParam[0]['value']
trMat = [[0.914375, 0.0, 0.0, 0.0], [0.0, 0.914375, 0.0, 0.0], [0.0, 0.0, 0.914375, 0.0], [-1636.25, -2544.86, -80.6958, 1.0]]
trMat = np.asmatrix(trMat)
#id = np.matlib.identity(4)
#for i in enumerate(id):
#    id[i] = trMat    

# Load flowdata
path2 = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1'
graph2 = spatialgraph.SpatialGraph()
graph2.read(os.path.join(path2,'spatialGraph_RIN.am'))

# Get node list
nodeFile2 = os.path.join(path2,'nodeList.dill')
if True:
#if not os.path.isfile(nodeFile2):
    print 'Generating node list...'
    nodeList2 = graph2.node_list()
            
    print('Pickling node list...')
    with open(nodeFile2,'wb') as fo:
        pickle.dump(nodeList2,fo)
else:
    with open(nodeFile2 ,'rb') as fo:
        nodeList2 = pickle.load(fo)

coords = graph.get_data('VertexCoordinates')
coords2 = graph2.get_data('VertexCoordinates')
for i,co in enumerate(coords):
    co = [co[0]*trMat[0,0],co[1]*trMat[1,1],co[2]*trMat[2,2]]
    co = [co[0]+trMat[3,0],co[1]+trMat[3,1],co[0]+trMat[3,2]]
    #co = [co[0]*trMat[0,0],co[1]*trMat[1,1],co[2]*trMat[2,2]]
    coords[i] = co

m1 = np.mean(coords,axis=0)
m2 = np.mean(coords2,axis=0)
dif = m2 - m1

coords = coords + dif

#import pdb
#pdb.set_trace()

red_node2 = []
mapping = []
for i,node in enumerate(nodeList):
    #cShift = rn.coords + dif
    print '{} of {}'.format(i,len(nodeList))
    #co = rn.coords
    #co = [co[0]*trMat[0,0],co[1]*trMat[1,1],co[2]*trMat[2,2]]
    #co = [co[0]+trMat[3,0],co[1]+trMat[3,1],co[0]+trMat[3,2]]
    co = coords[node.index]
    for i,co2 in enumerate(coords2):
        curDif = np.asarray(co2) - np.asarray(co)
        mg = np.linalg.norm(curDif)
        if i==0:
            mn = mg
            mnInd = i
        if mg<mn:
            mn = mg
            mnInd = i
    mapping.append([node.index,mnInd])
    if node in red_nodes:
        red_node2.append(nodeList2[mnInd])

redIndex = np.asarray([x.index for x in red_node2])
for node in nodeList2:
    if node in red_node2:
        node.add_scalar('red bond',1)
        nred += 1
        red_nodes.append(node)
    else:
        node.add_scalar('red bond',0)

#import pdb
#pdb.set_trace()
graph3 = graph2.node_list_to_graph(nodeList2)
graph3.write(os.path.join(path,'red_bonds.am'))