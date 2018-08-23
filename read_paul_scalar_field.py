# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:18:40 2017

@author: simon
"""

import os
import numpy as np

#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\Pressure mmHg'
#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T - Post-VDA\1\Pressure mmHg'
#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\Perfusion ml_min_100g'
#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T - Post-VDA\1\Perfusion ml_min_100g'

#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\SW1222\1\Perfusion'
#path = r'C:\Users\simon\Dropbox\Perfusion'
path = r'C:\Users\simon\Dropbox\SW1222\Pre-VDA\Pressure - mmHg'
path = r'C:\Users\simon\Dropbox\SW1222\Pre-VDA\Perfusion - ml_min_100g'
path = r'C:\Users\simon\Dropbox\Glioma\Interstitial_Perfusion'
path = r'C:\Users\simon\Dropbox\LS147T (1)\Pre-VDA\Pressure (mmHg)'
path = r'C:\Users\simon\Downloads\LSpre_IFV\Data'

path = r'C:\Users\simon\Dropbox\Pre-VDA SW1222\IFP\Data'
#r'C:\Users\simon\Dropbox\Glioma\Perfusion'
#r'C:\Users\simon\Dropbox\Glioma\Velocity\Data'

files = [os.path.join(path,f) for f in os.listdir(path) if f.endswith('.txt')]
files.sort(key=lambda f: int(list(filter(str.isdigit, f))))

nslice = len(files)

#pixSize = [148.8, 150.2, 150.0]
pixSize = [140.0, 140.0, 140.0]

for i,f in enumerate(files):
    with open(f,'r') as fo:
        print('Reading {}'.format(f))
        cur = fo.read()
    cur = [x for x in cur.split('\n') if x!='']
    tmp = [x for x in cur[0].split('\t') if x!='']
    
    if i==0:
        ncol = len(tmp)
        nrow = len(cur)
        grid = np.zeros([nslice,nrow,ncol],dtype='float')
        bbox = [0.,pixSize[0]*nslice, 0.,pixSize[1]*nrow, 0.,pixSize[2]*ncol]
        bboxStr = '{} {} {} {} {} {}'.format(0.,pixSize[0]*nslice,0.,pixSize[1]*nrow,0.,pixSize[2]*ncol)
    
<<<<<<< HEAD
    #print nrow,len(cur)
=======
    print(nrow,len(cur))
>>>>>>> a680e5459e9a481499a8c12d0867a6763350cfc5
    if nrow!=len(cur):
        import pdb
        pdb.set_trace()
    assert len(cur)==nrow
    
    for j in range(0,nrow,1):

        try:
            tmp = np.asarray([x for x in cur[j].split('\t') if x!=''],dtype='float')
            assert len(tmp)==ncol
            #import pdb
            #pdb.set_trace()
<<<<<<< HEAD
            #print np.max(tmp)
            #tmp = log(tmp);
            #tmp[isnan(tmp)] = 0.;
            grid[i,j,:] = tmp #np.abs(tmp)
        except Exception,e:
            print 'Error ({},{}): {}'.format(i,j,e)
=======
            print(np.max(tmp))
            grid[i,j,:] = np.abs(tmp)
        except Exception as e:
            print('Error ({},{}): {}'.format(i,j,e))
>>>>>>> a680e5459e9a481499a8c12d0867a6763350cfc5
           
#from matplotlib import pyplot
#pyplot.figure()
#pyplot.imshow(np.squeeze(grid[:,30,:]))
           
    #grid[grid<1e-12] = 0.
           
print('Min/max: {} {}'.format(np.min(grid),np.max(grid)))

# Write Amira mesh
from pymira import mesh
m = mesh.Mesh()
m.set_lattice_data(grid)
m.set_bounding_box(bbox)
ofile = os.path.join(path,'pressure.am')
m.write(ofile)
print('Written to {}'.format(ofile))
