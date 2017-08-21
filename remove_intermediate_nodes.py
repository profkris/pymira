# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 18:54:55 2017

@author: simon

Script for removing non-branching nodes and replacing them with edgepoints
Required when converting from Paul Sweeney's data format 

"""
from pymira import spatialgraph
import pickle
import os

amdata = spatialgraph.SpatialGraph()

#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\' 
#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\' 
#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'

dir_ = r'C:\Users\simon\Dropbox\170606_Ben Vessel Networks\C1M3\2%'
fo = os.path.join(dir_,'spatialGraph.am')

#dir_ = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T'
#fo = os.path.join(dir_,'fix_graph_30_LRG_NET.am')

ofile = os.path.join(dir_ , 'spatialGraph_RIN.am')

#dir_ = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\LS#1\\'
#fo = dir_+'LS1_spatialGraph_scaled.am'
#ofile = dir_+'LS1_spatialGraph_RIN.am'

print 'Loading graph...'
amdata.read(fo)
print('Graph loaded')

editor = spatialgraph.Editor()
new_graph = editor.remove_intermediate_nodes(amdata,path=dir_)
import pdb
pdb.set_trace()
#new_graph.sanity_check(deep=True)
new_graph.write(ofile)
#new_graph.write_node_list()