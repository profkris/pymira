# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:59:46 2017

@author: simon

Performs blind optical deconvolution (BOD) by
maximising the high-frequency regions of the power spectrum for each z-projection
Data must be ordered with the cutting dimension ('z') as the first dimension of the numpy array
and x,y (in plane images) as the second and third dimensions. Z should be first cut at 0 and last
cut at the end (-1).

"""

import numpy as np
import matplotlib.pyplot as pyplot
from scipy import signal
from scipy.optimize import minimize

def obj_func(beta,cur=None,verbose=False):
    
    if verbose:
        print beta
    cp = decon(cur,beta)
    
    ft = (np.abs(np.fft.fft(cp)))**2
    dt = 1./dims[0]
    freqs = np.fft.fftfreq(cp.size,dt)
    idx = np.argsort(freqs)
    ft = ft[idx]
    freqs = freqs[idx]
    ft0 = np.sum(ft[-20:-1])
    return -ft0
    
def decon(cur,beta):
    
    z = np.arange(0,cur.size)
    ir = np.exp(-np.abs(z)*beta)
    
    # Pad and flip data then deconvolve impulse response
    curPad = np.append(np.zeros(cur.size-1),cur)
    curPad = np.flip(curPad,0)
    cp, remainder = signal.deconvolve(curPad,ir)
    cpflip = np.flip(cp,0)
    cpflip = cpflip[0:cur.size]
    cpflip = np.clip(cpflip,0,np.max(cpflip))
    
    return cpflip

def simulate_bod_data():
    
    dims = [128,128,128]
    grid = np.zeros(dims)
    
    # Embed sphere into grid
    radius = 10
    offset = np.asarray([100,64,64])
    
    for i in range(radius*2):
        for j in range(radius*2):
            for k in range(radius*2):
                #r = radius - i/2.
                x = np.asarray([offset[0]-radius+i,offset[1]-radius+j,offset[2]-radius+k])
                if np.abs(np.linalg.norm(x-offset))<=radius:
                    grid[np.int(x[0]),np.int(x[1]),np.int(x[2])] = 1
    
    # Bleed-through simulation
    grid_p = np.zeros(dims)
    beta = np.zeros([dims[1],dims[2]]) + 0.1
    
    # Step through z-stack backwards and project signal onto slices
    # further up stack, with exponential decay
    for i in range(dims[0]-1,0,-1):
        for ip in range(i,0,-1):
            grid_p[ip,:,:] += grid[i,:,:]*np.exp(-np.abs(ip-i)*beta)
            
    return grid_p,beta
        
def bod(data,show_result=False,verbose=False):
    
    dims = data.shape    
    grid_c = np.zeros(dims)
    beta_est = np.zeros([dims[1],dims[2]])

    for i in range(dims[1]):
        for j in range(dims[2]):
            
            cur = data[:,i,j]
            
            res = minimize(obj_func,1.0,args=(cur,verbose),method='L-BFGS-B',bounds=[(0.01,2.)])
            if not res.success:
                print('Error at {},{}: {} '.format(i,j,res.message))
            else:
                grid_c[:,i,j] = decon(cur,res.x)
                beta_est[i,j] = res.x

    if show_result:
        ind = 64
        fig = pyplot.figure()
        a = fig.add_subplot(1,2,1)
        pyplot.imshow(data[:,ind,:])
        pyplot.axis('off')
        a = fig.add_subplot(1,2,2)
        pyplot.imshow(grid_c[:,ind,:])
        pyplot.axis('off')
        
    return grid_c,beta_est

def main():
    data,beta_true = simulate_bod_data()
    bod_data,beta = bod(data)
    
    beta_dif = np.abs(beta-beta_true)
    
if __name__ == "__main__":
    main()