# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:18:40 2017

@author: simon
"""

import os
import numpy as np

path = r'C:\Users\simon\Dropbox\45-20mmHg_perf154\IFV'
path = r'C:\Users\simon\Downloads\LSpre_IFV_Vector'
path = r'C:\Users\simon\Dropbox\Pre-VDA SW1222\IFV_velocity'

v_order = ['u','v','w']
#v_order = ['v','w','u']
path_uvw = [os.path.join(path,vo) for vo in v_order]

# Get data size
for i,p in enumerate(path_uvw):
    files = [os.path.join(p,f) for f in os.listdir(p) if f.endswith('.txt')]
    files.sort(key=lambda f: int(filter(str.isdigit, f)))

    nslice = len(files)
    
    for j,f in enumerate(files):
        with open(f,'r') as fo:
            print 'Reading {}'.format(f)
            cur = fo.read()
        cur = [xt for xt in cur.split('\n') if xt!='']
        tmp = [xt for xt in cur[0].split('\t') if xt!='']
        
        if i==0 and j==0:
            ncol = len(tmp)
            nrow = len(cur)
            pixSize = [50.0, 50.0, 50.0]
            #dims = [nslice,nrow,ncol,3] # x,y,z,3
            dims = [nrow,ncol,nslice,3] # y,z,x,3
            bbox = [0.,pixSize[0]*dims[0], 0.,pixSize[1]*dims[1], 0.,pixSize[2]*dims[2]]
            bboxStr = '{} {} {} {} {} {}'.format(bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5])
        
        #print nrow,len(cur)
        if nrow!=len(cur):
            import pdb
            pdb.set_trace()
        assert len(cur)==nrow

# Read data
grid = np.zeros(dims,dtype='float')
for i,p in enumerate(path_uvw):
    files = [os.path.join(p,f) for f in os.listdir(p) if f.endswith('.txt')]
    files.sort(key=lambda f: int(filter(str.isdigit, f)))

    nslice = len(files)
    
    for j,f in enumerate(files):
        with open(f,'r') as fo:
            print 'Reading {}'.format(f)
            cur = fo.read()
        cur = [x for x in cur.split('\n') if x!='']
        tmp = [x for x in cur[0].split('\t') if x!='']
        
        for k in xrange(0,nrow,1):
    
            try:
                tmp = np.asarray([x for x in cur[k].split('\t') if x!=''],dtype='float')
                if i in [0,1]:
                    #grid[i,j,:,i2] = tmp #np.abs(tmp)
                    grid[k,:,j,i] = -tmp # y,z,x,3
                else:
                    #grid[i,j,:,i2] = tmp #np.abs(tmp)
                    grid[k,:,j,i] = tmp # y,z,x,3
            except Exception,e:
                print 'Error ({},{},{}): {}'.format(i,j,k,e)
           
#from matplotlib import pyplot
#pyplot.figure()
#pyplot.imshow(np.squeeze(grid[:,30,:]))
           
    #grid[grid<1e-12] = 0.
           
#print 'Min/max: {} {}'.format(np.min(grid),np.max(grid))
           
# Write Amira mesh
from pymira import vector_field
m = vector_field.VectorField()
m.set_lattice_data(grid)
m.set_bounding_box(bbox)
ofile = os.path.join(path,'interstitial_velocity.am')
m.write(ofile)
print 'Written to {}'.format(ofile)
