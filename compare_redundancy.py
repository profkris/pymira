# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:43:32 2017

@author: simon
"""

import dill as pickle
import numpy as np
import os
from matplotlib import pyplot as plt

dir_ = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\SW1222\1'
#dir_ = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1'
odir = os.path.join(dir_,'test')

eDir = os.path.join(odir,'crawl_calcs')
files = os.listdir(eDir)
nfiles = len(files)
print('Loading crawl calc results ({} files)...'.format(len(files)))

mdif = np.zeros(0)
mrdif = np.zeros(0)
mn = np.zeros(0)
conn = np.zeros(0)

for fi,f in enumerate(files):
    print('Reading file {} of {}: {}'.format(fi+1,nfiles,f))
    if True: #try:
        with open(os.path.join(eDir,f),'rb') as fo:
            rec = pickle.load(fo)
        if len(rec)==2:
            curEdges,ind = rec[0],rec[1]
            redundancy = None
        elif len(rec)==4:
            curEdges,ind,redundancy,connectivity = rec[0],rec[1],rec[2],rec[3]
            
        mdif = np.append(mdif,redundancy[:,1])
        mrdif = np.append(mrdif,redundancy[:,2])
        mn = np.append(mn,redundancy[:,3])
        conn = np.append(conn,connectivity)

import pdb
pdb.set_trace()
print 'Mean difference: {} +/- {}'.format(np.mean(mdif),np.std(mdif))
print 'Mean rel difference: {} +/- {}'.format(np.mean(mrdif),np.std(mrdif))
print 'Mean n: {} +/- {}'.format(np.mean(mn),np.std(mn))
print 'Mean conn: {} +/- {}'.format(np.mean(conn[conn>0.]),np.std(conn[conn>0.]))

plt.figure()
hist = plt.hist(mrdif,bins=np.arange(0,2.5,0.01))
ax = plt.gca()
ax.set_yscale('log')

plt.figure()
hist = plt.hist(conn[conn>0.],bins=np.arange(0,0.25,0.001))
ax = plt.gca()
ax.set_yscale('log')

hist = plt.hist([mrdifls,mrdifsw],bins=np.arange(0,2.5,0.05),normed=True)
ax = plt.gca()
ax.set_yscale('log')

hist = plt.hist([mdifls,mdifsw],bins=np.arange(0,10000,100),normed=True)
ax = plt.gca()
ax.set_yscale('log')

x1 = np.log(connls[connls>0.])
x2 = np.log(connsw[connsw>0.])
hist = plt.hist([x1,x2],bins=np.arange(-8.,0,0.2),normed=True)
ax = plt.gca()
#ax.set_xscale('log')
ax.set_yscale('log')