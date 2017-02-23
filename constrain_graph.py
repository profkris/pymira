# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 18:54:55 2017

@author: simon
"""
from pymira import spatialgraph
import pickle

amdata = spatialgraph.SpatialGraph()

dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\' 
fo = dir_+'spatialGraph.am'
ofile = dir_ + 'spatialGraph_constrained.am'

#dir_ = 'C:\\Anaconda2\\Lib\\site-packages\\pymira\\'
#fo = dir_ + 'circle.am'
#ofile = dir_ + 'circle_RIN.am'

print 'Loading graph...'
amdata.read(fo)
print('Graph loaded')
#print('Performing sanity check...')
#amdata.sanity_check(deep=True)

amdata.constrain_nodes(xrange=[1000.,2000.],yrange=[0.,2000.])
#amdata.constrain_nodes(xrange=[0.,10.],yrange=[0.,10.])
#amdata.remove_field('Pressure')
#amdata.remove_field('Stress')
#amdata.remove_field('Flow')
#amdata.remove_field('Radii')
amdata.sanity_check(deep=True)
#ofile = dir_+'circle_constrained.am'
amdata.write(ofile)