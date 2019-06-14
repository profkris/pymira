# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:58:10 2017

@author: simon

Utility for converting Amira SpatiaGraphs into STL files (DEVELOPMENTAL)

"""

import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import math

#def make_epilpsoid(rad1,rad2)

def graph_to_stl(graph):

    coords = graph.get_data('VertexCoordinates')
    edgeConn = graph.get_data('EdgeConnectivity')
    nconn = edgeConn.shape[0]
    #edgeCoords = graph.get_data('EdgePointCoordinates')
    nodeList = graph.node_list()
    edges = graph.edges_from_node_list(nodeList)
    
    thickness = 1.
    
    #nedge = edgeConn.shape[0]
    #verts = np.empty(0)
    verts = None
    faces = []
    
    tubes = []
    connectionDefined = np.zeros(nconn,dtype='int')
    
    # Loop through edges
    for ei,edge in enumerate(edges):
        print('Edge {} of {}'.format(ei,len(edges)))
        pts = edge.coordinates
        radius = edge.get_scalar('Radii')
        strtNode = nodeList[edge.start_node_index]
        endNode = nodeList[edge.end_node_index]
        
        strtBranchInd = [i for i,e in enumerate(strtNode.edges) if e.index==edge.index][0]
        strtRev = strtNode.edge_indices_rev[strtBranchInd]
        endBranchInd = [i for i,e in enumerate(endNode.edges) if e.index==edge.index][0]
        endRev = endNode.edge_indices_rev[endBranchInd]
        #print(strtRev,endRev,pts)
        
        #print strtNode.nconn
        vec1,vec2 = None,None        
        if True: #strtNode.nconn==1: # Terminating edge
            strt_face_angle = 0.
        elif strtNode.nconn==2:
            # Find the next node that the edge connects to (and which isn't the end node)
            #other_conn = strtNode.connecting_node[strtNode.connecting_node!=endNode.index]
            if strtNode.connecting_node[0]==strtNode.connecting_node[1]:
                other_conn = strtNode.connecting_node[0]
            else:
                other_conn = [x for x in strtNode.connecting_node if x!=endNode.index][0]
            # Calculate the distance of next and previous edges
            length1 = np.linalg.norm(coords[other_conn]-coords[strtNode.index])
            vec1 = (coords[other_conn]-coords[strtNode.index])/length1
            length2 = np.linalg.norm(coords[strtNode.index]-coords[endNode.index])
            vec2 = (coords[strtNode.index]-coords[endNode.index])/length2
            # Calculate angle between the two connecting edges
            costh = np.dot(vec2,vec1)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
            strt_face_angle = -np.arccos(np.clip(costh,-1.,1.))
            strt_face_angle = np.rad2deg(strt_face_angle) / 2.
            #print(other_conn,strtNode.index,endNode.index,coords[other_conn],coords[strtNode.index],coords[endNode.index])
            #import pdb
            #pdb.set_trace()
        else:
            strt_face_angle = 0.
        #if strt_face_angle!=90.:
        #    pass
        print('Start angle: {} , v1 {} v2 {} '.format(strt_face_angle,vec1,vec2))
        if not np.isfinite(strt_face_angle):
            import pdb
            pdb.set_trace()
            
        #print endNode.nconn
        vec1,vec2 = None,None
        if True: #endNode.nconn==1:
            end_face_angle = 0.
        elif endNode.nconn==2:
            if endNode.connecting_node[0]==endNode.connecting_node[1]:
                other_conn = endNode.connecting_node[0]
            else:
                other_conn = [x for x in endNode.connecting_node if x!=strtNode.index][0]
            length1 = np.linalg.norm(coords[endNode.index]-coords[other_conn])
            vec1 = (coords[endNode.index]-coords[other_conn])/length1
            length2 = np.linalg.norm(coords[endNode.index]-coords[strtNode.index])
            vec2 = (coords[endNode.index]-coords[strtNode.index])/length2
            costh = np.dot(vec2,vec1)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
            end_face_angle = np.arccos(np.clip(costh,-1.,1.))
            end_face_angle = np.rad2deg(end_face_angle)
            #import pdb
            #pdb.set_trace()
        else:
            end_face_angle = 0.
        print('End angle: {} , v1 {} v2 {} '.format(end_face_angle,vec1,vec2))
        if not np.isfinite(end_face_angle):
            import pdb
            pdb.set_trace()

        fa = np.zeros(pts.shape[0])
        fa[0] = strt_face_angle
        fa[pts.shape[0]-1] = end_face_angle
        grad = np.zeros([pts.shape[0],3])
        
        for i in range(pts.shape[0]-1):
            length = np.linalg.norm(pts[i]-pts[i+1])
            vec = (pts[i+1]-pts[i])/length
            center = np.mean(np.vstack((pts[i],pts[i+1])), axis=0)
            #print strt_face_angle,end_face_angle
            
            if False: #i<pts.shape[0]-2:
                length2 = np.linalg.norm(pts[i+1]-pts[i+2])
                vec2 = (pts[i+2]-pts[i+1])/length2
                fa[i+1] = np.arccos(np.dot(vec,vec2)/(np.linalg.norm(vec)*np.linalg.norm(vec2)))
                fa[i+1] = -np.rad2deg(fa[i+1])
                
            if i==0:
                grad[i] = vec
            elif i<pts.shape[0]-2:
                grad[i] = pts[i+1]-pts[i]
            #verts,faces = make_tube(radius[i:i+1]*2.,[thickness,thickness],lengthorientation=vec,center=center,outer_only=True)
            #print(fa[i],fa[i+1])
            if np.isnan(fa[i]) or np.isnan(fa[i+1]):
                import pdb
                pdb.set_trace()
            
            #print('Points:{}{}, center:{}, orient:{}'.format(pts[i],pts[i+1],center,vec))
            curRad = np.mean(radius[i:i+1]*2)

            #newMesh = tube_mesh(curRad,thickness,length,orientation=vec,center=center,outer_only=True,start_face_angle=-fa[i],end_face_angle=fa[i+1])
            nlength = 2
            alpha = 360.
            nalpha = 16
            oVerts,oFaces = make_cylinder(curRad,length*1.2,nlength,alpha,nalpha,center,vec,start_face_angle=-fa[i],end_face_angle=fa[i+1])

            if verts is None:
                verts = oVerts.astype('float')
                faces = oFaces
            else:
                #import pdb
                #pdb.set_trace()
                verts = np.append(verts,oVerts.astype('float'),axis=0)
                nf = np.max(faces)
                faces = np.append(faces,oFaces+nf+1,axis=0)
            #verts.append(oVerts)
            #faces.append(oFaces)

    faces = np.asarray(faces)
    verts = np.asarray(verts)
            
    newMesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            newMesh.vectors[i][j] = verts[f[j],:]                    
    
    #tubes.append(newMesh)

    #netMesh = mesh.Mesh(np.concatenate([t.data for t in newMesh]))
    #plot_mesh(newMesh) 
    return newMesh
    
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
    return np.dot(points,np.transpose(R))
    
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
    except Exception as e:
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

def make_cylinder(radius, length, nlength, alpha, nalpha, center, orientation,start_face_angle=0.,end_face_angle=0.,doubleSided=True):
    
    origin, xaxis, yaxis, zaxis = [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
    
    if end_face_angle>180:
        return None
    if end_face_angle>90.:
        end_face_angle -= 180.    

    #Create the length array
    if type(length) is not list:
        length = [length]
    if len(length)==1:
        I = np.linspace(0, length[0], nlength)
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
        posInd = px
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
        posInd = px
        py = X
        pz = Y
    else:
        print('Incorrect number of elements in radius!')
        return None

    points = np.vstack(( px, py, pz )).T
    if doubleSided:
        points = np.concatenate((points,points),axis=0)
    
    # Define triangular faces
    faces = []
    for i in range(nalpha):
        if i<nalpha-1:
            faces.append([i,i+1,nalpha+i])
            faces.append([i+1,nalpha+i+1,nalpha+i])
        else:
            faces.append([i,0,nalpha+i])
            faces.append([0,nalpha,nalpha+i])
    if doubleSided:
        for i in range(nalpha,2*nalpha):
            if i<nalpha-1:
                faces.append([2*nalpha+i,i+1,i])
                faces.append([2*nalpha+i,2*nalpha+i+1,i+1])
            else:
                faces.append([2*nalpha+i,0,i])
                faces.append([2*nalpha+i,2*nalpha,0])
    
    faces = np.asarray(faces)

    #Orient tube to new vector
    U = align_vector_rotation(xaxis,orientation)
    points_p = np.dot(U,np.transpose(points)).T
    #points_p = points
        
    #Shift to center
    shift = np.array(center) - np.mean(points_p, axis=0)
    points_p += shift
    
    U = U.astype('float')
    Up = np.linalg.inv(U)
    
    if start_face_angle!=0.:
        R = rot_matrix(np.deg2rad(start_face_angle),zaxis)
        pStartInds = [i for i,p in enumerate(px) if p==I[0]]
        pStart = points_p[pStartInds]
        rotPoint = np.mean(pStart,axis=0)
        pStart -= rotPoint
        #pStartR = rotate_points(pStart,Up[:3,:3])
        pStartR = rotate_points(pStart,R[:3,:3]) 
        #pStartR = rotate_points(pStartR,U[:3,:3]) 
        pStartR += rotPoint
        
        #pStart -= rotPoint
        #pStartR = rotate_points(pStart,R[:3,:3]) + rotPoint
        points_p[pStartInds] = pStartR
        
    if end_face_angle!=0.:
        R = rot_matrix(np.deg2rad(end_face_angle),zaxis)  
        pEndInds = [i for i,p in enumerate(px) if p==I[-1]]
        pEnd = points_p[pEndInds]
        rotPoint = np.mean(pEnd,axis=0)
        pEnd -= rotPoint
        #pEndR = rotate_points(pEnd,Up[:3,:3])
        pEndR = rotate_points(pEnd,R[:3,:3]) 
        #pEndR = rotate_points(pEndR,U[:3,:3]) 
        pEndR += rotPoint
        points_p[pEndInds] = pEndR
    
    return points_p,faces
      
def tube_mesh(radius, thickness, length,orientation=[1,0,0],center=[0,0,0],outer_only=False,start_face_angle=0.,end_face_angle=0.):
    
    nlength = 2
    alpha = 360.
    nalpha = 16
    
    oVerts,oFaces = make_cylinder(radius,length,nlength,alpha,nalpha,center,orientation,start_face_angle=start_face_angle,end_face_angle=end_face_angle)
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

#gfile = 'C:\\Users\\simon\\Dropbox\\Mesentery\\Flow2AmiraPressure.am'
#gfile = 'C:\\Anaconda2\\Lib\\site-packages\\pymira\\test_graph.am'
pyplot.close('all')
#gfile = r'C:\Anaconda2\Lib\site-packages\pymira\test_join.am'
#gfile = r'C:\Users\simon\Dropbox\segmentation_database\vessels\VDA Colorectal cancer\Control\LS\LS#2\LS2_spatial_graph.am'
gfile = r'C:\Users\simon\Dropbox\RoyalSocSims\subvol_4339_3691_4147.am'
#gfile = r'C:\Users\simon\Dropbox\RoyalSocSims\spatialGraph_RIN_flow_flag.am'
#gfile = r'C:\Anaconda2\Lib\site-packages\pymira\tests\test_network_simple.am'
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
#dir_ = 'C:\\Users\\simon\\Dropbox\\'
dir_ = r'C:\Users\simon\Dropbox\RoyalSocSims'
##combined.save(dir_+'cylinder.stl')
#combined.save(dir_+'mesentry.stl')
#ofile = r'C:\Users\simon\Dropbox\RoyalSocSims\spatialGraph_RIN_flow_flag.stl'
#ofile = r'C:\Users\simon\Dropbox\RoyalSocSims\spatialGraph_RIN_flow_flag.stl'
#ofile = r'C:\Anaconda2\Lib\site-packages\pymira\tests\test_network_simple.stl'
ofile = r'C:\Users\simon\Dropbox\RoyalSocSims\subvol.stl'
combined.save(ofile)
