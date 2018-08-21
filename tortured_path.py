# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 07:37:28 2018

@author: simon
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
#import mathutils #.geometry as geometry

import math as m
def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az

def asSpherical(xyz):
    #takes list xyz (single coord)
    x       = xyz[0]
    y       = xyz[1]
    z       = xyz[2]
    r       =  np.sqrt(x*x + y*y + z*z)
    theta   =  np.arccos(z/r) #*180/ np.pi #to degrees
    phi     =  np.arctan2(y,x) #*180/ np.pi
    return [r,theta,phi]
    
def asCartesian(rthetaphi):
    #takes list rthetaphi (single coord)
    r       = rthetaphi[0]
    theta   = rthetaphi[1] #* pi/180 # to radian
    phi     = rthetaphi[2] #* pi/180
    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    return [x,y,z]

def tortured_path(a,b,tortuosity=None):
    
    def point_to_line(r,p,q):
        x = p-q
        n = q-p
        if np.dot(n,p-r)>0.: # Closest point is p
            return -1,np.dot(p-r,p-r)
        if np.dot(n,r-q)>0.: # Closest point is q
            return 1,np.dot(r-q,r-q)
        t = np.dot(r-q, x)/np.dot(x, x)
        return 0,np.linalg.norm(t*(p-q)+q-r)
        
    norm = np.linalg.norm
        
    a = np.asarray(a)
    b = np.asarray(b)
        
    length = np.linalg.norm(a-b)
    if tortuosity is None:
        tortuosity = length/50.
    #if sigma is None:
    #    sigma = length / 10. #float(npoints)
    
    path = [] #np.zeros([npoints,3])
    path.append(a)
    #path[1,:] = b
    
    dbPrev = length
    lim1 = 0. #sigma/10. #length/100. # 0. #length / 10. # Gradient of line - needs to head downhill (>=0)
    lim2 = length/2.#3. # Distance from central line
    #for i in range(1,npoints-1):
    i = 1
    count = 0
    count_lim = 1e5
    while True:
        prevPoint = path[i-1]
        theta_r = np.deg2rad(np.random.normal(loc=0.,scale=30))
        phi_r = np.deg2rad(np.random.normal(loc=0.,scale=30))
        l_r = np.random.uniform(low=length/50.,high=length/10.)
        #l_r = length / 10.
        
        #newPoint = prevPoint + np.random.normal(scale=sigma,size=3)
        #newPoint = np.asarray([l * np.sin(theta) * np.cos(phi), 
        #        l * np.sin(theta) * np.sin(phi), 
        #        l * np.cos(phi)])
        #newPoint += prevPoint
        #import pdb
        #pdb.set_trace()
        if i>1:
            v1 = path[i-1] - path[i-2]
        else:
            v1 = b - a
            
        #r,az,el = cart2sph(v1[0],v1[1],v1[2])
        r,az,el = asSpherical(v1)
        theta = az + theta_r
        phi = el + phi_r
        #import pdb
        #pdb.set_trace()
        #displacement = np.asarray([l_r * np.sin(theta) * np.cos(phi), 
        #                       l_r * np.sin(theta) * np.sin(phi), 
        #                       l_r * np.cos(phi)])
        displacement = asCartesian([l_r,theta,phi])
        newPoint = prevPoint + displacement
        #import pdb
        #pdb.set_trace()
        #print newPoint
        #v2 = newPoint - path[i-1]
            #theta = np.rad2deg(np.arccos(np.dot(v1,v2) / np.dot(norm(v1),norm(v2))))
            
        #Calculate angle between proposed vessel segment and end point
        #vb = a - b #- a #newPoint
        #theta_b = np.rad2deg(np.arccos(np.dot(v2,vb) / np.dot(norm(v2),norm(vb))))
            
        loc,dpl = point_to_line(newPoint,a,b)

        db = np.linalg.norm(newPoint-b)
        acc = dpl < np.random.normal(loc=0.,scale=tortuosity)#/10.)
        acc2 = dbPrev-db > 0. #np.random.normal(loc=sigma,scale=2.*sigma)
        
        if acc and acc2: # and acc2: #(dbPrev-db)>lim1 and acc: #dpl<lim2 and acc: # and \
            #print('Closer: {}, length: {}'.format(db<dbPrev,np.linalg.norm(displacement)))
            path.append(newPoint)
            i += 1
            dbPrev = db
        if db<l_r: # and theta_b>180-t2 and theta_b<180+t2:
            path.append(b)
            break
        count += 1
        if count>count_lim:
            print('Path generation failed!')
            return None
    return np.asarray(path)
    
#a = [0.,0.,0.]
#b = [1.,1.,1.]
#sigma = 0.01
#path = tortured_path(a,b)
#print(len(path))
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#x = np.asarray([t[0] for t in path])
#y = np.asarray([t[1] for t in path])
#z = np.asarray([t[2] for t in path])
#ax.plot(x,y,z,color='r')
#ax.scatter(x,y,z,color='r')
#ax.plot([a[0],b[0]],[a[1],b[1]],[a[2],b[2]],color='b')
#ax.scatter([0],[0],[0],facecolors='b')
#ax.scatter([1],[1],[1],facecolors='r')
#ax.set_xlim([0,1])
#ax.set_ylim([0,1])
#ax.set_zlim([0,1])
#plt.show()