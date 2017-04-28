# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 18:54:55 2017

@author: simon

Script for removing non-branching nodes and replacing them with edgepoints
Required when converting from Paul Sweeney's data format 

"""
from pymira import spatialgraph
import pickle

amdata = spatialgraph.SpatialGraph()

dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\' 
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\' 
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
fo = dir_+'spatialGraph.am'
ofile = dir_ + 'spatialGraph_RIN.am'

dir_ = 'C:\\Users\\simon\\Dropbox\\VDA_1_lectin\\Control\LS#1\\'
fo = dir_+'LS1_spatialGraph_scaled.am'
ofile = dir_+'LS1_spatialGraph_RIN.am'

print 'Loading graph...'
amdata.read(fo)
print('Graph loaded')

editor = spatialgraph.Editor()
new_graph = editor.remove_intermediate_nodes(amdata,path=dir_)
new_graph.sanity_check(deep=True)
new_graph.write(ofile)
#new_graph.write_node_list()