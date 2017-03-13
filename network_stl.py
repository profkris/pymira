# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:58:10 2017

@author: simon
"""

import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import math

def make_epilpsoid(rad1,rad2)

def make_cylinder(radius, length, nlength, alpha, nalpha, center, orientation):

    #Create the length array
    I = np.linspace(0, length, nlength)

    #Create alpha array avoid duplication of endpoints
    #Conditional should be changed to meet your requirements
    if int(alpha) == 360:
        A = np.linspace(0, alpha, num=nalpha, endpoint=False)/180*np.pi
    else:
        A = np.linspace(0, alpha, num=nalpha)/180*np.pi

    #Calculate X and Y
    X = radius * np.cos(A)
    Y = radius * np.sin(A)

    #Tile/repeat indices so all unique pairs are present
    #pz = np.tile(I, nalpha)
    #px = np.repeat(X, nlength)
    #py = np.repeat(Y, nlength)
    pz = np.repeat(I, nalpha)
    px = np.tile(X, nlength)
    py = np.tile(Y, nlength)

    points = np.vstack(( pz, px, py )).T
    
    # Define triangular faces
    faces = []
    for i in range(nalpha):
        if i<nalpha-1:
            faces.append([i,i+1,nalpha+i])
            faces.append([nalpha+i,nalpha+i+1,i+1])
        else:
            faces.append([i,0,nalpha+i])
            faces.append([nalpha+i,nalpha,0])
    faces = np.asarray(faces)
        
#    oldfaces = np.array([\
#           [0,1,nalpha],
#           [nalpha,nalpha+1,1],
#           [1,2,nalpha+1],
#           [nalpha+1,nalpha+2,2],
#           [2,3,nalpha+2],
#           [nalpha+2,nalpha+3,3],
#           [3,0,nalpha+3],
#           [nalpha+3,nalpha,0],
#                    ])

    #Shift to center
    shift = np.array(center) - np.mean(points, axis=0)
    points += shift

    #Orient tube to new vector

    #Grabbed from an old unutbu answer
    def rotation_matrix(axis,theta):
        a = np.cos(theta/2)
        b,c,d = -axis*np.sin(theta/2)
        return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    ovec = orientation / np.linalg.norm(orientation)
    cylvec = np.array([1,0,0])

    if np.allclose(cylvec, ovec):
        return points,faces

    #Get orthogonal axis and rotation
    oaxis = np.cross(ovec, cylvec)
    rot = np.arccos(np.dot(ovec, cylvec))

    R = rotation_matrix(oaxis, rot)
    return points.dot(R),faces
    
def make_tube(radius, thickness, length):
    
    center = [0,0,0]
    orient = [1,0,0]
    nlength = 2
    alpha = 360.
    nalpha = 16
    
    outer_cylinder = make_cylinder(radius,length,nlength,alpha,nalpha,center,orient)
    inner_cylinder = make_cylinder(radius+thickness,length,nlength,alpha,nalpha,center,orient)
    return [outer_cylinder,inner_cylinder]

# Define the 8 vertices of the cube
vertices = np.array([\
    [-1, -1, -1],
    [+1, -1, -1],
    [+1, +1, -1],
    [-1, +1, -1],
    [-1, -1, +1],
    [+1, -1, +1],
    [+1, +1, +1],
    [-1, +1, +1]])
# Define the 12 triangles composing the cube
faces = np.array([\
    [0,3,1],
    [1,3,2],
    [0,4,7],
    [0,7,3],
    [4,5,6],
    [4,6,7],
    [5,1,2],
    [5,2,6],
    [2,3,6],
    [3,7,6],
    [0,1,5],
    [0,5,4]])
    
#CYLINDER
vertices,faces = make_cylinder(1, 5, 2, 360, 16, [0,0,0], [1,0,0])

# Create the mesh
cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        cube.vectors[i][j] = vertices[f[j],:]                

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# Render the cube faces
#for m in meshes:
#    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(m.vectors))
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(cube.vectors))

# Auto scale to the mesh size
meshes = [cube]
scale = np.concatenate([m.points for m in meshes]).flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show() 

dir_ = 'C:\\Users\\simon\\Dropbox\\'
cube.save(dir_+'cylinder.stl')