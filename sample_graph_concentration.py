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

#path = r'D:\160113_paul_simulation_results\LS147T\1\ca1\vascular_recon_log'
#path = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\gd\vascular_recon'
#nlpath = r'D:\160113_paul_simulation_results\LS147T\1\ca1'
#nlpath = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\gd'
path = r'D:'
nlpath = r'D:\gd'

files = []
for f in os.listdir(path):
    if f.endswith('.am'):
        files.append(f)

nfile = len(files)
#Sort files
files.sort(key=lambda f: int(list(filter(str.isdigit, f))))

meshes = []
time = []
for f in files:
    m = spatialgraph.SpatialGraph()
    m.read(os.path.join(path,f))
    meshes.append(m)
    time.append(m.get_parameter_value('Time'))
    
time = [float(t) for t in time]
nt = len(time)

sample = 10

def plot_conc(sample,meshes,time):
    conc = np.zeros(len(meshes))
    for i,m in enumerate(meshes):
        curnode = m.get_node(sample)
        for edge in curnode.edges:
            try:
                conc[i] += edge.get_scalar('Concentration')[0]
            except Exception as e:
                import pdb
                pdb.set_trace()
        conc[i] /= float(len(curnode.edges))
        
    vals = conc
    
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
    y = np.append(np.asarray([0.,0.]),np.exp(valsSrt))
    y *= 1e-3 / vox_vol #mmol/L
    
    #ysm = scipy.signal.medfilt(y,9)
    
    pyplot.plot(x,y)
    ax.set_ylabel('Concentration (mM)')
    ax.set_xlabel('Time (minutes)')
    #pyplot.xlim(0,100)
    pyplot.show()
    
plot_conc(sample,meshes,time)