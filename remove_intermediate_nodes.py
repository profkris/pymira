# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 18:54:55 2017

@author: simon
"""
from pymira import spatialgraph
import pickle

amdata = spatialgraph.SpatialGraph()

dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\' 
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\' 
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
fo = dir_+'spatialGraph.am'
ofile = dir_ + 'spatialGraph_RIN.am'

print 'Loading graph...'
amdata.read(fo)
print('Graph loaded')

editor = spatialgraph.Editor()
new_graph = editor.remove_intermediate_nodes(amdata)
new_graph.sanity_check(deep=True)
new_graph.write(ofile)