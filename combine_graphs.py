# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 09:37:47 2017

@author: simon
"""

def merge_graphs(f1,f2,fo):
    import pymira.amiramesh as am
    a1 = am.AmiraMesh()
    a1.read(f1)
    if a1.get_parameter_value('ContentType')!='HxSpatialGraph':
        print('{} is not a SpatialGraph file! {}'.format(f1,a1.get_parameter_value('ContentType')))
        return
        
    a2 = am.AmiraMesh()
    a2.read(f2)
    if a2.get_parameter_value('ContentType')!='HxSpatialGraph':
        print('{} is not a SpatialGraph file! {}'.format(f2,a1.get_parameter_value('ContentType')))
        return
    
    dif1  = list(set(a1.fieldNames) - set(a2.fieldNames))
    dif2  = list(set(a2.fieldNames) - set(a1.fieldNames))
    
    for fName in dif2:
        f = a2.get_field(fName)
        a1.add_field(f)
        
    a1.write(fo)

dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\' 
fp = dir_+'Press2Amira.txt'
ff = dir_+'Flow2Amira.txt'
fs = dir_+'Stress2Amira.txt'
fv = dir_+'Velocity2Amira.txt'
fo = dir_+'spatialGraph.am'
merge_graphs(fp,ff,fo)
merge_graphs(fo,fv,fo)
merge_graphs(fo,fs,fo)