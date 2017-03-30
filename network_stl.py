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

#def make_epilpsoid(rad1,rad2)

def graph_to_stl(graph):

    #coords = graph.get_data('VertexCoordinates')
    #edgeConn = graph.get_data('EdgeConnectivity')
    #edgeCoords = graph.get_data('EdgePointCoordinates')
    nodeList = graph.node_list()
    edges = graph.edges_from_node_list(nodeList)
    
    thickness = 1.
    
    #nedge = edgeConn.shape[0]
    
    tubes = []
    for edge in edges:
        pts = edge.coordinates
        radius = edge.get_scalar('Radii')
        strtNode = nodeList[edge.start_node_index]
        endNode = nodeList[edge.end_node_index]
        
        strtBranchInd = [i for i,e in enumerate(strtNode.edges) if e.index==edge.index][0]
        strtRev = strtNode.edge_indices_rev[strtBranchInd]
        endBranchInd = [i for i,e in enumerate(endNode.edges) if e.index==edge.index][0]
        endRev = endNode.edge_indices_rev[endBranchInd]
        
        if strtNode.nconn==1:
            strt_ba = 0.
        elif strtNode.nconn==2:
            strt_ba = np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
            deg = np.rad2deg(rad)
        else:
            import pdb
            pdb.set_trace()
        
        for i in range(pts.shape[0]-1):
            length = np.linalg.norm(pts[i]-pts[i+1])
            vec = (pts[i+1]-pts[i])/length
            center = np.mean(np.vstack((pts[i],pts[i+1])), axis=0)
            import pdb
            pdb.set_trace()
            #verts,faces = make_tube(radius[i:i+1]*2.,[thickness,thickness],lengthorientation=vec,center=center,outer_only=True)
            
            #print('Points:{}{}, center:{}, orient:{}'.format(pts[i],pts[i+1],center,vec))
            tubes.append(tube_mesh(np.mean(radius[i:i+1]*2),thickness,length,orientation=vec,center=center,outer_only=True))

    netMesh = mesh.Mesh(np.concatenate([t.data for t in tubes]))
    plot_mesh(netMesh) 
    return netMesh
    
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)    
    
def angle_between(v1, v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    
def GetSkew(x):
    return np.asmatrix(
           [[0, -x[0,2], x[0,1]],
            [x[0,2], 0, -x[0,0]],
            [-x[0,1], x[0,0], 0]])

def vector_rotation(v1, v2):
    
    mat = np.asmatrix
    norm = np.linalg.norm 
    
    v1 = np.asmatrix(v1)
    v2 = np.asmatrix(v2)
    
    # rotation vector
    w = np.cross(v1,v2)
    if norm(w)!=0.:
        w = w / norm(w)
    else:
        w = np.asarray([0.,0.,0.])
    
    w_hat = GetSkew(w)
    
    #rotation angle
    cos_tht = v1.T*v2/(norm(v1)*norm(v2))
    tht = np.squeeze(mat(np.arccos(cos_tht)))
    sin_tht = np.sin(tht)
    w_hat2 = np.square(w_hat)
    tht1 = 1. - np.cos(tht)
    R = np.identity(3) + w_hat*sin_tht.T + (w_hat2*tht1.T)
    
    return R
    
def rotate_points(points,R):
    return np.dot(points,np.transpose(U))
    
def align_vector_rotation(A,B):
    
    """ To align vector A to B:
    B = np.transpose(U*A.T))
    """
    
    mat = np.asmatrix
    norm = np.linalg.norm
    
    A = mat(A)/norm(A)
    B = mat(B)/norm(B)
    #print('A={}'.format(A))
    #print('B={}'.format(B))
    
    if np.all(A==B):
        #print('All equal')
        return np.identity(3)
    elif np.all(A==-B):
        #print('All negative')
        return -np.identity(3)
    
    G = np.asmatrix(
         [[A*B.T, mat(-norm(np.cross(A,B))), 0.],
         [norm(np.cross(A,B)), A*B.T,  0.],
         [0., 0., 1.],
         ])

    Fi = np.asarray([ A , (B - (A*B.T)*A) / norm((B - (A*B.T)*A)) , np.cross(B,A) ])
    Fi = np.matrix(np.squeeze(Fi))
    Fi = np.nan_to_num(Fi)
    try:
        U = Fi * G * np.linalg.inv(Fi)
    except Exception,e:
        print(e)
        import pdb
        pdb.set_trace()
    
    return U
    
def rot_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    >>> R = rotation_matrix(math.pi/2, [0, 0, 1], [1, 0, 0])
    >>> numpy.allclose(numpy.dot(R, [0, 0, 0, 1]), [1, -1, 0, 1])
    True
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = numpy.random.random(3) - 0.5
    >>> point = numpy.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(-angle, -direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> I = numpy.identity(4, numpy.float64)
    >>> numpy.allclose(I, rotation_matrix(math.pi*2, direc))
    True
    >>> numpy.allclose(2, numpy.trace(rotation_matrix(math.pi/2,
    ...                                               direc, point)))
    True

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = unit_vector(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[ 0.0,         -direction[2],  direction[1]],
                      [ direction[2], 0.0,          -direction[0]],
                      [-direction[1], direction[0],  0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = np.array(point[:3], dtype=np.float64, copy=False)
        M[:3, 3] = point - np.dot(R, point)
    return M

def make_cylinder(radius, length, nlength, alpha, nalpha, center, orientation):

    #Create the length array
    if type(length) is not list:
        length = [length]
    if len(length)==1:
        I = np.linspace(0, length, nlength)
    elif len(length)==nlength-1:
        I = length
    else:
        print('Incorrect number of elements in length!')
        return None

    #Create alpha array avoid duplication of endpoints
    #Conditional should be changed to meet your requirements
    if int(alpha) == 360:
        A = np.linspace(0, alpha, num=nalpha, endpoint=False)/180*np.pi
    else:
        A = np.linspace(0, alpha, num=nalpha)/180*np.pi

    #Calculate X and Y
    if type(radius) is not list:
        radius = [radius]
    if len(radius)==1:
        X = np.asarray(radius) * np.cos(A)
        Y = np.asarray(radius) * np.sin(A)

        #Tile/repeat indices so all unique pairs are present
        px = np.repeat(I, nalpha)
        py = np.tile(X, nlength)
        pz = np.tile(Y, nlength)
    elif len(radius)==nlength:
        X = np.zeros(0)
        Y = np.zeros(0)
        for r in radius:
            X = np.append(X,np.asarray(r) * np.cos(A))
            Y = np.append(Y,np.asarray(r) * np.sin(A))
    
            #Tile/repeat indices so all unique pairs are present
        px = np.repeat(I, nalpha)
        py = X
        pz = Y
    else:
        print('Incorrect number of elements in radius!')
        return None

    points = np.vstack(( px, py, pz )).T
    
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
    
    #import pdb
    #pdb.set_trace()

    #Orient tube to new vector
    origin, xaxis, yaxis, zaxis = [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
    U = align_vector_rotation(xaxis,orientation)
    points_p = np.dot(U,np.transpose(points)).T
    #points_p = points
    #print U
        
    #Shift to center
    shift = np.array(center) - np.mean(points_p, axis=0)
    points_p += shift
    
    return points_p,faces
      
def tube_mesh(radius, thickness, length,orientation=[1,0,0],center=[0,0,0],outer_only=False):
    
    nlength = 2
    alpha = 360.
    nalpha = 16
    
    oVerts,oFaces = make_cylinder(radius,length,nlength,alpha,nalpha,center,orientation)
    if not outer_only:
        iVerts,iFaces = make_cylinder(radius+thickness,length,nlength,alpha,nalpha,center,orientation)
    
    # Create the mesh
    outer_mesh = mesh.Mesh(np.zeros(oFaces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(oFaces):
        for j in range(3):
            outer_mesh.vectors[i][j] = oVerts[f[j],:]                    
        
    if not outer_only:                
        inner_mesh = mesh.Mesh(np.zeros(iFaces.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate(iFaces):
            for j in range(3):
                inner_mesh.vectors[i][j] = iVerts[f[j],:]
        
        combinedMesh = mesh.Mesh(np.concatenate([inner_mesh.data, outer_mesh.data]))
    else:
        combinedMesh = outer_mesh

    return combinedMesh

def plot_mesh(meshes):
    
    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    
    # Render the cube faces
    #for m in meshes:
    #    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(m.vectors))
    #axes.add_collection3d(mplot3d.art3d.Poly3DCollection(inner_mesh.vectors))
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(meshes.vectors))
    
    # Auto scale to the mesh size
    scale = np.concatenate(meshes.points).flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)
    
    # Show the plot to the screen
    pyplot.show() 

gfile = 'C:\\Users\\simon\\Dropbox\\Mesentery\\Flow2AmiraPressure.am'
#gfile = 'C:\\Anaconda2\\Lib\\site-packages\\pymira\\test_graph.am'
from pymira import spatialgraph
graph = spatialgraph.SpatialGraph()
graph.read(gfile)
combined = graph_to_stl(graph)

#CYLINDER
#radius = 20.
#thickness = 5.
#length = 50.
#it,ot = make_tube(radius,thickness,length)
#inner_vertices = it[0]
#inner_faces = it[1]
#outer_vertices = ot[0]
#outer_faces = ot[1]
#
#
#plot_mesh(combined)
#
dir_ = 'C:\\Users\\simon\\Dropbox\\'
##combined.save(dir_+'cylinder.stl')
combined.save(dir_+'mesentry.stl')