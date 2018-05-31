# -*- coding: utf-8 -*-
"""
Created on Tue Aug 01 13:50:02 2017

@author: simon
"""

import numpy as np
import matplotlib.pyplot as plt

conc = np.exp(np.asarray([m.get_field('Lattice')['data'] for m in meshes]))
concl = conc[-1,:,:,:]
inds = np.where(concl>0.)
#print len(inds[0])
coords = np.asarray([[inds[0][i],inds[1][i],inds[2][i]] for i,tmp in enumerate(inds[0])])
from scipy.spatial import ConvexHull
hull = ConvexHull(coords)
edges = list(zip(*coords))
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

#nt = conc.shape[0]

shTr = np.identity(3)
kv = [1.,0.75,0.5,0.25,0.1]
#colV = ['blue','red','green','black','grey']
#colV = [1.,0.9,0.8,0.7,0.6]
colV = [str(x) for x in np.linspace(0.6,0.,num=len(kv))]
fig = plt.figure()
import scipy.ndimage as ndimage
#pltOffset = 0 #5

mxTarg = 0.5e-9
sdTarg = 0.1e-9

for j,k in enumerate(kv):
    #k = 0.5
    shTr[0,0],shTr[1,1],shTr[2,2] = k,k,k
    com = np.mean(coords,axis=0)
    dif = coords - com
    difP = np.dot(dif,shTr)
    coordsP = difP + com
    coordsP = coordsP.astype('int')
    vec = np.zeros(nt)
    count = 0.
    #mask_er = ndimage.binary_erosion(mask_er)#.astype(a.dtype)
    for i in hull.vertices:
        ax.scatter(coordsP[i,0],coordsP[i,1],coordsP[i,2],color='red')
        vec += conc[:,coordsP[i,0],coordsP[i,1],coordsP[i,2]]
        count += 1

    vec /= float(count)
    vec /= np.max(vec)
    vec *= np.random.normal(loc=mxTarg,scale=sdTarg)
    
    #vox_size = np.asarray([150.,150.,150.]) # um
    #vox_vol = np.product(vox_size*1e-6) * 1000. # L
    #vec *= 1e-3 / vox_vol #mmol/L
    #plt.plot(time_min,perim)
    #yoff = np.mean(vec[pltOffset:pltOffset+6])
    #vec -= yoff
    #print yoff

    from scipy import optimize
    p0 = [0.02, 0.1, 0.001, 250.] # Initial guess for the parameters
    fitfunc = lambda p, x: np.clip(p[0]*(1.-np.exp(-p[1]*(x-p[3]))) * np.exp(-p[2]*(x-p[3])),0.,1e6)
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    p1, success = optimize.leastsq(errfunc, p0[:], args=(time,vec))
    timeStr = 3. #0.6
    timep = np.linspace(0.,time[-1],num=len(time))
    time_min = np.asarray(timep)*timeStr / 60.
    #plt.plot(time_min-pltOffset,fitfunc(p1, time), color=colV[j], linewidth=2) 
    plt.plot(time_min+2,vec,color=colV[j], linewidth=2)
plt.xlim([0,15]) #np.max(time_min)])
#plt.ylim([-1e-8,1.3e-7])