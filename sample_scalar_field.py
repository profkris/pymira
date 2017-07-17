# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 15:02:21 2017

@author: simon
"""

import os
from pymira import spatialgraph,mesh
from matplotlib import pyplot
import numpy as np
import scipy.signal

import matplotlib
matplotlib.rcParams.update({'font.size': 22})

#path = r'D:\160113_paul_simulation_results\LS147T\1\ca1_kt0p00001\interstitial_concentration_recon'
#path = r'D:\160113_paul_simulation_results\LS147T\1'
#path = r'D:\gd\interstitial_concentration_recon'
path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\parker\interstitial_concentration_recon'
#path = r'D:\gd_kt001\interstitial_concentration_recon'
files = []
for f in os.listdir(path):
    if f.endswith('.am'):
        files.append(f)

nfile = len(files)
#Sort files
files.sort(key=lambda f: int(filter(str.isdigit, f)))

meshes = []
time = []
for f in files:
    m = mesh.Mesh()
    m.read(os.path.join(path,f))
    meshes.append(m)
    time.append(m.get_parameter_value('Time'))
    
time = [float(t) for t in time]
nt = len(time)
    
d0 = meshes[0].get_field('Lattice')['data']
dims = d0.shape

data = meshes[-1].get_field('Lattice')['data']

mxInds = list(np.unravel_index(data.argmax(), data.shape))
targ = 10 # nth largest value
targVal = np.partition(data.flatten(), -targ)[-targ]
targInds = np.where(data==targVal)
targInds = [targInds[0][0],targInds[1][0],targInds[2][0]]

sample = [int(dims[0]/2.),int(dims[1]/2.),int(dims[2]/2.)]
sample = [4,24,39] # Hot spot
sample = [24,39,4]
sample = mxInds
sample = targInds

offset = [-189.5,-200.,-200.]
pixsize = [150.519,152.659,151.243]
coords = [4448,2292,9818]
#coords = [3050.,2293.,7944.] # central low uptake
sample = np.asarray([int((x-offset[i])/pixsize[i]) for i,x in enumerate(coords)])

def process_sample(data,meshes,time,sample,logs=True,smooth=False,width=9):
    
    nt = len(time)
    vals = np.zeros(nt)
    vals = [meshes[i].get_field('Lattice')['data'][sample[0],sample[1],sample[2]] for i,t in enumerate(time)] # umol
    
    vox_size = np.asarray([150.,150.,150.]) # um
    vox_vol = np.product(vox_size*1e-6) * 1000. # L
        
    #pyplot.figure()
    srt = np.argsort(time)
    time = np.asarray(time)#[srt]
    timeSrt = time[srt] # [t for i,srt in srt]
    y = np.asarray(vals)[srt]
    
    if logs:
        y = np.exp(y)
    
    y *= 1e-3 / vox_vol #mmol/L
    
    if smooth:
        y = scipy.signal.medfilt(y,width)
        
    return timeSrt,y
    

def plot_lattice_conc(data,meshes,time,sample,logs=True,smooth=False):

    vals = np.zeros(nt)
    vals = [meshes[i].get_field('Lattice')['data'][sample[0],sample[1],sample[2]] for i,t in enumerate(time)] # umol
    
    vox_size = np.asarray([150.,150.,150.]) # um
    vox_vol = np.product(vox_size*1e-6) * 1000. # L
        
    #pyplot.figure()
    srt = np.argsort(time)
    time = np.asarray(time)#[srt]
    timeSrt = time[srt] # [t for i,srt in srt]
    valsSrt = np.asarray(vals)[srt]
    
    #subs = [i for i,t in enumerate(timeSrt) if t<45 or t>60]
    #valsSrt = valsSrt[subs]
    #timeSrt = timeSrt[subs]
    
    fig, ax = pyplot.subplots(figsize=(7,5))
    x = np.append(np.asarray([-10.,0.]),timeSrt)

    if logs:
        y = np.append(np.asarray([0.,0.]),np.exp(valsSrt))
    else:
        y = np.append(np.asarray([0.,0.]),valsSrt)
    y *= 1e-3 / vox_vol #mmol/L
    
    if smooth:
        y = scipy.signal.medfilt(y,9)
    
    pyplot.plot(x/60.,y)
    ax.set_ylabel('Concentration (mM)')
    ax.set_xlabel('Time (minutes)')
    pyplot.xlim(0,20)
    pyplot.show()
    return x,y
    
def fit_sample(time,cur,sample):
    from scipy import optimize
    fitfunc = lambda p, x: np.clip(p[0]*(1.-np.exp(-p[1]*(x-p[3]))) * np.exp(-p[2]*(x-p[3])),0.,1e6)
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    p0 = [0.02, 0.1, 0.001, 250.] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(time,cur))
    pyplot.figure()
    pyplot.plot(time/60., cur, "ro", time/60., fitfunc(p1, time), "r-") 
    #print('Sample: {}'.format(s))
    print('Fit params: {}'.format(p1))
        
#x,y = plot_lattice_conc(data,meshes,time,sample,smooth=True)
x,y = process_sample(data,meshes,time,sample,logs=True,smooth=False,width=9)
fit_sample(x,y/10.,sample)