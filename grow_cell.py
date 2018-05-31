# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 08:26:24 2017

@author: simon

First attempt at balloon cell growth simulation (INCOMPLETE)

"""

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt

def cell_volume(coords):
    #x = [c[0] for c in coords]
    #y = [c[1] for c in coords]
    #z = np.zeros(len(x))
    #z = [c[2] for c in coords]
    return ConvexHull(coords).volume
    
def plot_cell(coords):
    plt.figure()
    x = [c[0] for c in coords]
    y = [c[1] for c in coords]
    plt.plot(x,y)

init_pressure = 1.

nCoords = [0.,0.,0.]
center = [0.,0.]
nperim = 200
alpha = np.arange(0,2*np.pi,2*np.pi/np.float(nperim))
perimVelocity = []
perimPosition = []
perimConn = []

for i in range(nperim+1):
    if i==nperim:
        perimConn.append([i,0])
    else:
        perimConn.append([i,i+1])

v_init = 1. # um/s
for i in range(nperim):
    vx = v_init * np.cos(alpha[i])
    vy = v_init * np.sin(alpha[i])
    #vz = 0.
    perimVelocity.append([vx,vy])#,vz])
    perimPosition.append([0.+center[0],0.+center[1]])#,0.])
    
perimVelocity = np.asarray(perimVelocity)
perimPosition = np.asarray(perimPosition)

nstep = 100
dt = 1 #s
time = np.arange(0,nstep,dt)

box_size = np.asarray([50,50,50])
box_center = np.asarray([0,0,0])

for tInd,t in enumerate(time):
    newPos = []
    for pos,vel in zip(perimPosition,perimVelocity):
        newPos.append(pos + vel*dt)
    perimPosition = np.asarray(newPos)
    
    eta = []
    strain_energy = []
    for i in range(nperim):
        if i==0:
            neigh = [nperim-1,i+1]
        elif i==nperim-1:
            neigh = [nperim-2,0]
        else:
            neigh = [i-1,i+1]
            
        pn0 = np.asmatrix(perimPosition[neigh[0]])
        pn1 = np.asmatrix(perimPosition[neigh[1]])
        pn = np.asmatrix(perimPosition[i])
        dx = np.abs(pn1-pn0)
        dx[dx==0] = 1e-6
        
        #displ = np.asarray([pn0-pn,pn1-pn])
        displ = (pn-center).T
    
        B = np.asmatrix([ 
              [1./dx[0,0], 0     ],
              [0., 1./dx[0,1]    ],
              [1./dx[0,1],1./dx[0,0] ]
        ])
        
        eta.append(np.dot(B,displ)) # strain
        
        strain_energy.append(eta[i]*dx[0,0]*dx[0,1])
        
        # Pressure
        volume = cell_volume(perimPosition)
        pressure = init_pressure / volume
        import pdb
        pdb.set_trace()
        #tang = np.mean([pn1-pn,pn-pn0])
        normal = [ [-(pn1[0,1]-pn0[0,1]),   pn1[0,0]-pn0[0,0]  ],
                   [  pn1[0,1]-pn0[0,1] , -(pn1[0,0]-pn0[0,0]) ]
                   ]
        #normal = normal / np.linalg(normal)
        pForce = pressure*dx
    
    eta = np.asmatrix(np.squeeze(eta))
    
    #eta = np.asarray(eta)
    #strain_energy = 0.5 * np.sum(strain_energy)
        
    #stress
    lam = 1.
    mu = 1.
    D = np.asmatrix([
          [lam+2.*mu, lam, 0.],
          [lam, lam+2.*mu, 0.],
          [0., 0., mu],
    ])
    sigma = D*eta.T
        
    volume = cell_volume(perimPosition)
    pressure = init_pressure / volume
    pEnergy = pressure * volume
    
    energy = 0.
    for i in range(nperim):
        pn = np.asmatrix(perimPosition[i])
        dx = np.abs(pn1-pn0)
        dx[dx==0] = 1e-6
        
        displ = pn.T
        dx = np.abs(pn1-pn0)
        dx[dx==0] = 1e-6
        energy += (displ.T * B.T * D * B * displ) * dx[0,0] * dx[0,1] * 0.5  
        
    energy -= pEnergy
    
    print(('Volume {} Pressure {} Energy {} pEnergy {}'.format(volume,pressure,energy,pEnergy)))

plot_cell(perimPosition)
