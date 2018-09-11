# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 08:37:38 2018

@author: simon
"""

import numpy as np
import nibabel as nib
import scipy.ndimage.interpolation as interp
import scipy.ndimage.measurements as measurements 
from pymira import tortured_path
import sys
import os

class VesselNetwork(object):
    
    def __init__(self):
        self.start = [0.,0.,0.]
        self.l = 0.
        self.r = 0.
        self.n = 0.
        self.theta = 0.
        self.phi = 0.
        self.term = False
        self.parent = -1
        self.level = -1
        self.index = -1
        self.occ = False
        self.proc = None
        self.netIndex = 1
        
class GrowthProcess(object):
    
    def __init__(self,name,**kwargs):
        self.name = name
        self.gamma = kwargs.get('gamma')
        self.beta = kwargs.get('beta')
        self.f = kwargs.get('f')
        self.nf = kwargs.get('nf')
        self.theta = kwargs.get('theta')
        self.phi = kwargs.get('phi')
        self.sigmaL = kwargs.get('sigmaL')
        self.sigmaD = kwargs.get('sigmaD')
        self.sigmaT = kwargs.get('sigmaT')
        self.sigmaP = kwargs.get('sigmaP')

def asCartesian(rthetaphi):
    #takes list rthetaphi (single coord)
    r       = rthetaphi[0]
    theta   = rthetaphi[1] #* pi/180 # to radian
    phi     = rthetaphi[2] #* pi/180
    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    return [x,y,z]

def vessel_network_3d():

  # Number of levels in vessel hierarchy
  nlevel = 9
  # Number of branches per vessel
  nb = 2
  # Number of vessels to model
  nvessel = np.sum([nb**i for i in range(nlevel)])
  print('No. vessels: {}'.format(nvessel))
  
  # Define a structure to hold network information
  v = [VesselNetwork() for i in range(nvessel)]
  
  # Power law exponent for vessel diameters: D0^gamma = D1^gamma + D2^gamma
  # D0 = parent diameter, D1 & D2 = daughter vessel diameters.
  gamma = 3.
  # Symmetry variable: beta = D1 / D2 (0<beta<1)
  beta = 1. #0.95
  # Ratio of parent length to daughter length: l_daughter = f x l_parent
  f = 0.9
  # Fraction of length at which node occurs
  nf = 0.9
  # Branching angles
  theta = 27. * (np.pi/180.) # 27. rad (25.5<ba<28.5 for kidney; 25<ba<140 for tumour)
  phi = 10. * (np.pi/180.) # 110. rad
  
  # Randomness variables (uniform distribution)
  # Kidney
  sigmaL = 0.02
  sigmaD = 0.02
  sigmaT = 0.02
  sigmaP = 10.
  
  tortFactor = 0.5
  
  # Define objects for holding growth process information
  gp = GrowthProcess('normal',gamma=gamma,beta=beta,f=f,nf=nf,theta=theta,phi=phi, \
        sigmaL=sigmaL,sigmaD=sigmaD,sigmaT=sigmaT,sigmaP=sigmaP)
  
  # Populate network structure.
  # Starting vessel  
  v[0].occ = True
  v[0].l = 4500. #5000 #um
  v[0].r = 500. #um
  v[0].theta = 0. #rad
  v[0].phi = 0. #rad
  v[0].start = [0.,0.,0.] #x,y,z (um)
  v[0].startIndex = 0
  v[0].end = np.asarray(v[0].start) + asCartesian([v[0].l,v[0].theta,v[0].phi])
  v[0].endIndex = 1
  v[0].path = None
  v[0].term = True
  v[0].parent = -1
  v[0].index = 0
  v[0].level = 0
  if tortFactor<=0:
      path = np.asarray([v[0].start,v[0].end])
  else:
      tort = np.linalg.norm(np.asarray(v[0].start)-np.asarray(v[0].end))*tortFactor
      path = tortured_path.tortured_path(v[0].start,v[0].end,tortuosity=tort)
  v[0].path = path
  v[0].r = np.zeros(path.shape[0]) + v[0].r
  v[0].occ = True
  
  endloop = False
  count = 0
  nparent = 1
  pSubs = [0]
  level = 0
  level = level + 1

  count = 1
  node_count = 2

  while(endloop is False):
    for i in range(nparent):
      curV = v[pSubs[i]]
      curV.term = False
      node = curV.end
      tmp = (curV.r[-1]**gp.gamma) / (1.+(gp.beta**gp.gamma))
      d1 = np.exp(np.log(tmp)/gp.gamma)
      d2 = gp.beta * d1
      dnew = [d1,d2]
        
      nxt = np.where([x.occ is False for x in v])
      nFree = len(nxt[0])
      if nFree<nb:
          plot_vessel_network(v)
          endloop =True
          break
      nxt = nxt[0]
      
      for j in range(len(dnew)):
        # Generate random numbers
        rndmL = np.random.uniform(low=-gp.sigmaL,high=gp.sigmaL)
        rndmD = np.random.uniform(low=-gp.sigmaD,high=gp.sigmaD)
        rndmT = np.random.uniform(low=-gp.sigmaT,high=gp.sigmaT)
        rndmP = np.random.uniform(low=-gp.sigmaP,high=gp.sigmaP)      
        
        newV = v[nxt[j]]
        newV.proc = curV.proc        
        newl = curV.l * f
        newV.l = newl + (newl*rndmL) #um
        newV.r = dnew[j] + (dnew[j]*rndmD)
        if j % 2==0:
            newV.theta = curV.theta + gp.theta + (gp.theta*rndmT)
        else: 
            newV.theta = curV.theta - gp.theta + (gp.theta*rndmT)
        if j % 2==0:
            newV.phi = curV.phi + gp.phi + (gp.phi*rndmP)
        else:
            newV.phi = curV.phi - gp.phi + (gp.phi*rndmP)
        newV.start = node #x,y (um)
        newV.end = node + asCartesian([newV.l,newV.theta,newV.phi])
        newV.term = True
        newV.parent = curV.index
        newV.level = level
        newV.index = count
        newV.startIndex = curV.endIndex
        newV.endIndex = node_count
        
        node_count += 1
        ntry = 10
        count1 = 0
        while True: # Try and connect to node using constrained random walk
            if tortFactor<=0:
                path = np.asarray([newV.start,newV.end])
            else:
                tort = np.linalg.norm(np.asarray(newV.start)-np.asarray(newV.end)) * tortFactor
                path = tortured_path.tortured_path(newV.start,newV.end,tortuosity=tort) # np.asarray([curV.start,node])
            if path is not None:
                if np.any(path[0]!=newV.start) or np.any(path[-1]!=newV.end):
                    import pdb
                    pdb.set_trace()
                break
            else:
                count1 += 1
                if count1 >= ntry:
                    import pdb
                    pdb.set_trace()
        newV.path = path
        if not np.isscalar(newV.r):
            import pdb
            pdb.set_trace()
        newV.r = np.zeros(path.shape[0]) + newV.r
        newV.occ = True
        v[nxt[j]] = newV
        count += 1

    pSubs = np.where([x.term is True for x in v])
    nparent = len(pSubs[0])
    pSubs = pSubs[0]
    free = np.where([x.occ==0 for x in v])
    nfree = len(free[0])
    free = free[0]
    if nfree==0:
        endloop = True
    #IF nParent*nb GE nfree THEN endloop = 1 
    level += 1
    
  v = [x for x in v if x.occ is True]
    
  return v
  
def coord_to_pix(coord,extent_um,dims,clip=True):
    """ Convert 3D spatial coordinates into 3D pixel indices """
    if clip:
        return np.asarray([np.clip(np.round((coord[i]-extent_um[i,0])*dims[i]/(extent_um[i,1]-extent_um[i,0])),0,dims[i]) for i in [0,1,2]],dtype='int')
    else:
        return np.asarray([(coord[i]-extent_um[i,0])*dims[i]/(extent_um[i,1]-extent_um[i,0]) for i in [0,1,2]],dtype='int')

def list_to_amira(v):
    
    from pymira import spatialgraph
    sg = spatialgraph.SpatialGraph(initialise=True)
    sg.add_field(name='Radius',marker='@5',
                              definition='POINT',type='float',
                              nelements=1,nentries=[0],data=None)
                              
    nnode = 0
    nedge = 0
    npoint = 0
    # Get starting nodes
    nodes = np.asarray([x.start for x in v])
    inds = np.asarray([x.startIndex for x in v])
    # Get terminating nodes
    nodes = np.vstack((np.asarray([x.end for x in v if x.term]),nodes))
    inds = np.append([x.endIndex for x in v if x.term],inds)
    # Sort to get unique nodes
    unq_inds = np.unique(inds)
    nnodes = unq_inds.shape[0]
    nodes_unq = np.zeros((nnodes,3))
    #import pdb
    #pdb.set_trace()
    for i in unq_inds:
        j = np.where(inds==i)
        nodes_unq[i,:] = nodes[j[0][0]]
        
    sg.fields[0]['data'] = np.asarray(nodes_unq)
        
    # Get connections
    conns = [[x.startIndex,x.endIndex] for x in v]
    sg.fields[1]['data'] = np.asarray(conns)
    
    edges = [np.asarray(x.path) for x in v]
    edges = np.concatenate(edges)
    npoints = np.asarray([x.path.shape[0] for x in v])
    radii = np.asarray([x.r for x in v])
    radii = np.concatenate(radii)

    sg.fields[2]['data'] = npoints
    sg.fields[3]['data'] = edges
    sg.fields[4]['data'] = radii
    
    #import pdb
    #pdb.set_trace()
    
    sg.definitions[0]['size'] = [nnodes]
    sg.definitions[1]['size'] = [len(conns)]
    sg.definitions[2]['size'] = [np.sum(npoints)]
    return sg    
    
def vector_align(v1,v2):
    a,b = (v1/ np.linalg.norm(v1)).reshape(3), (v2/ np.linalg.norm(v2)).reshape(3)
    v = np.cross(a,b)
    c = np.dot(a,b)
    
    if np.all(v==0.):
        if c<0:
            return np.eye(3)
        else:
            return np.eye(3)
    
    s = np.linalg.norm(v)
    I = np.identity(3)
    vXStr = '{} {} {}; {} {} {}; {} {} {}'.format(0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0)
    k = np.matrix(vXStr)
    r = I + k + np.matmul(k,k) * ((1 -c)/(s**2))
    return r

    #for i in xrange(item.shape[0]):
    #    item[i] = (np.dot(r, item[i]).reshape(3,1)).reshape(3)    
    
def vector_norm(data, axis=None, out=None):
    import numpy
    import math
    """Return length, i.e. Euclidean norm, of ndarray along axis.

    >>> v = numpy.random.random(3)
    >>> n = vector_norm(v)
    >>> numpy.allclose(n, numpy.linalg.norm(v))
    True
    >>> v = numpy.random.rand(6, 5, 3)
    >>> n = vector_norm(v, axis=-1)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=2)))
    True
    >>> n = vector_norm(v, axis=1)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
    True
    >>> v = numpy.random.rand(5, 4, 3)
    >>> n = numpy.empty((5, 3))
    >>> vector_norm(v, axis=1, out=n)
    >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
    True
    >>> vector_norm([])
    0.0
    >>> vector_norm([1])
    1.0

    """
    data = numpy.array(data, dtype=numpy.float64, copy=True)
    if out is None:
        if data.ndim == 1:
            return math.sqrt(numpy.dot(data, data))
        data *= data
        out = numpy.atleast_1d(numpy.sum(data, axis=axis))
        numpy.sqrt(out, out)
        return out
    else:
        data *= data
        numpy.sum(data, axis=axis, out=out)
        numpy.sqrt(out, out)
        
def decompose_matrix(matrix):
    import numpy
    import math
    """Return sequence of transformations from transformation matrix.

    matrix : array_like
        Non-degenerative homogeneous transformation matrix

    Return tuple of:
        scale : vector of 3 scaling factors
        shear : list of shear factors for x-y, x-z, y-z axes
        angles : list of Euler angles about static x, y, z axes
        translate : translation vector along x, y, z axes
        perspective : perspective partition of matrix

    Raise ValueError if matrix is of wrong type or degenerative.

    >>> T0 = translation_matrix([1, 2, 3])
    >>> scale, shear, angles, trans, persp = decompose_matrix(T0)
    >>> T1 = translation_matrix(trans)
    >>> numpy.allclose(T0, T1)
    True
    >>> S = scale_matrix(0.123)
    >>> scale, shear, angles, trans, persp = decompose_matrix(S)
    >>> scale[0]
    0.123
    >>> R0 = euler_matrix(1, 2, 3)
    >>> scale, shear, angles, trans, persp = decompose_matrix(R0)
    >>> R1 = euler_matrix(*angles)
    >>> numpy.allclose(R0, R1)
    True

    """
    # epsilon for testing whether a number is close to zero
    _EPS = numpy.finfo(float).eps * 4.0
    
    M = numpy.array(matrix, dtype=numpy.float64, copy=True).T
    if abs(M[3, 3]) < _EPS:
        raise ValueError('M[3, 3] is zero')
    M /= M[3, 3]
    P = M.copy()
    P[:, 3] = 0.0, 0.0, 0.0, 1.0
    if not numpy.linalg.det(P):
        raise ValueError('matrix is singular')

    scale = numpy.zeros((3, ))
    shear = [0.0, 0.0, 0.0]
    angles = [0.0, 0.0, 0.0]

    if any(abs(M[:3, 3]) > _EPS):
        perspective = numpy.dot(M[:, 3], numpy.linalg.inv(P.T))
        M[:, 3] = 0.0, 0.0, 0.0, 1.0
    else:
        perspective = numpy.array([0.0, 0.0, 0.0, 1.0])

    translate = M[3, :3].copy()
    M[3, :3] = 0.0

    row = M[:3, :3].copy()
    scale[0] = vector_norm(row[0])
    row[0] /= scale[0]
    shear[0] = numpy.dot(row[0], row[1])
    row[1] -= row[0] * shear[0]
    scale[1] = vector_norm(row[1])
    row[1] /= scale[1]
    shear[0] /= scale[1]
    shear[1] = numpy.dot(row[0], row[2])
    row[2] -= row[0] * shear[1]
    shear[2] = numpy.dot(row[1], row[2])
    row[2] -= row[1] * shear[2]
    scale[2] = vector_norm(row[2])
    row[2] /= scale[2]
    shear[1:] /= scale[2]

    if numpy.dot(row[0], numpy.cross(row[1], row[2])) < 0:
        numpy.negative(scale, scale)
        numpy.negative(row, row)

    angles[1] = math.asin(-row[0, 2])
    if math.cos(angles[1]):
        angles[0] = math.atan2(row[1, 2], row[2, 2])
        angles[2] = math.atan2(row[0, 1], row[0, 0])
    else:
        # angles[0] = math.atan2(row[1, 0], row[1, 1])
        angles[0] = math.atan2(-row[2, 1], row[1, 1])
        angles[2] = 0.0

    return scale, shear, angles, translate, perspective
    
def rotate_3d(vol,R,centre=[0,0,0]):
    
    tr = np.eye(4)
    tr[0:3,3] = [-x for x in centre]
    tri = np.linalg.inv(tr)
    
    Rh = np.eye(4)
    Rh[0:3,0:3] = R
    Rhi = np.linalg.inv(Rh)
    
    res = vol.copy() * 0.
    
    subs = np.where(vol>0.)
    npix = len(subs[0])
    if npix==0:
        return res
    subs_n = np.asarray([[subs[0][i],subs[1][i],subs[2][i],1] for i in range(npix)])
    
    arr = np.asarray    
    subs_np1 = arr([np.dot(tr,s) for s in subs_n])
    subs_np2 = arr([np.dot(Rhi,s) for s in subs_np1])
    subs_np3 = arr([np.dot(tri,s) for s in subs_np2])
    subs_np = subs_np3.astype('int') # arr([subs_np3[:,0],subs_np3[:,1],subs_np3[:,2]])
    keepInd = arr([i for i,x in enumerate(subs_np) if np.all(x>=0) and np.all(x[0:3]<vol.shape)]) 
    
    try:
        subs_np = subs_np[keepInd]
        subs_n = subs_n[keepInd]
        res[subs_np[:,0],subs_np[:,1],subs_np[:,2]] = vol[subs_n[:,0],subs_n[:,1],subs_n[:,2]]
    except Exception as e:
        print(e)
        import pdb
        pdb.set_trace()
                
    return res
                
def vessel_grid(extent=None,dims=None,r=10.,c=np.asarray([0.,0.,0.]),R=None,rot=[0.,0.,0.],l=10.,inside=1,outside=0):
    
    """ 
    extent: bounding box for matrix (um)
    dims: grid dimensions
    """

    grid = np.zeros(dims,dtype='int')
    r_pix = [np.int(np.round(r*dims[i]/(extent[i,1]-extent[i,0]))) for i in range(3)]
    l_pix = [np.clip(np.int(np.round(l*dims[i]/(extent[i,1]-extent[i,0]))),0,dims[i]) for i in range(3)]
    print('Lpix {}, Rpix {}'.format(l_pix,r_pix))
    
    grid[:] = outside
    # Create a circle of size r (um), in centre of grid
    circ = np.zeros(dims[0:2])+outside
    ct = np.asarray([dims[i]/2 for i in range(3)])
    for x in range(dims[0]/2-r_pix[0],dims[0]/2+r_pix[0]):
        for y in range(dims[1]/2-r_pix[1],dims[0]/2+r_pix[1]):
            #for z in range(0,l_pix[2]):
            if True:
                xy = np.asarray([x,y])
                if np.linalg.norm(xy-ct[0:2])<=r_pix[0]: # and xyz[2]<l_pix[2]:
                    circ[x,y] = inside

    # Copy circle to create cylinder of length l, in centre of grid
    for z in range((dims[2]-l_pix[2])/2,(dims[2]+l_pix[2])/2):
        grid[:,:,z] = circ
        
    #subgrid = 

    #import scipy.ndimage
    #Ri = np.linalg.inv(R)
    #Ri = np.linalg.inv(R)
    offset = [(x/2)-1 for x in dims] #- grid.dot(Ri)
    #grid = scipy.ndimage.affine_transform(grid,Ri,offset=offset)
    grid = rotate_3d(grid,R,centre=offset)

    # Rotate cylinder to correct orientation
    #grid = interp.rotate(grid,rot[0],axes=[1,2],reshape=False)
    #grid = interp.rotate(grid,rot[1],axes=[0,2],reshape=False)
    #grid = interp.rotate(grid,rot[2],axes=[0,1],reshape=False)
    
    # Shift to correct position
    com_pix = coord_to_pix(c,extent,dims)
    com_grid = measurements.center_of_mass(grid)
    com_dif = com_pix - com_grid
    grid = interp.shift(grid,com_dif)
    
    grid[grid>inside/2] = inside
                    
    return grid

def vessel_grid_fast(extent=None,dims=None,r=10.,c=np.asarray([0.,0.,0.]),R=None,rot=[0.,0.,0.],l=10.,inside=1,outside=0):
    
    """ 
    extent: bounding box for matrix (um)
    dims: grid dimensions
    """

    r_pix = [np.int(np.round(r*dims[i]/(extent[i,1]-extent[i,0]))) for i in range(3)]
    l_pix = [np.clip(np.int(np.round(l*dims[i]/(extent[i,1]-extent[i,0]))),0,dims[i]) for i in range(3)]
    print('Lpix {}, Rpix {}'.format(l_pix,r_pix))
    
    #import pdb
    #pdb.set_trace()
    gsize = np.repeat(np.max([np.max(l_pix),np.max(2*r_pix)]) * 2,3)
    #gsize = np.asarray([x+1 if x%2==0 else x for x in gsize]) # Make dimensions odd
    grid = np.zeros(gsize,dtype='int')
    
    grid[:] = outside
    # Create a circle of size r (um), in centre of grid
    circ = np.zeros(gsize[0:2])+outside
    ct = np.asarray([gsize[i]/2 for i in range(3)])
    xl,xh = int(gsize[0]/2-r_pix[0]), int(gsize[0]/2+r_pix[0])
    yl,yh = int(gsize[1]/2-r_pix[1]), int(gsize[1]/2+r_pix[1])
    for x in range(xl,xh):
        for y in range(yl,yh):
            xy = np.asarray([x,y])
            if np.linalg.norm(xy-ct[0:2])<=r_pix[0]: # and xyz[2]<l_pix[2]:
                circ[x,y] = inside

    # Copy circle to create cylinder of length l, in centre of grid
    zl,zh = int((gsize[2]/2-l_pix[2])/2), int((gsize[2]/2+l_pix[2]/2))
    for z in range(zl,zh):
        grid[:,:,z] = circ

    # Rotation offset
    offset = [(x/2)-1 for x in gsize]
    grid = rotate_3d(grid,R,centre=offset)
    
    # Shift to correct position
    com_pix = coord_to_pix(c,extent,dims)
    #com_grid = measurements.center_of_mass(grid)
    res = np.zeros(dims,dtype='int')
    xl,xh = int(com_pix[0]-gsize[0]/2), int(com_pix[0]+gsize[0]/2)
    yl,yh = int(com_pix[1]-gsize[1]/2), int(com_pix[1]+gsize[1]/2)
    zl,zh = int(com_pix[2]-gsize[2]/2), int(com_pix[2]+gsize[2]/2)
    xlp,xhp = 0,gsize[0]
    ylp,yhp = 0,gsize[1]
    zlp,zhp = 0,gsize[2]
    if xh>dims[0]:
        xhp = xhp-(xh-dims[0])
        xh = dims[0]
    if yh>dims[1]:
        yhp = yhp-(yh-dims[1])
        yh = dims[1]
    if zh>dims[2]:
        zhp = zhp-(zh-dims[0])
        zh = dims[0]
    if xl<0:
        xlp = xlp - xl
        xl = 0
    if yl<0:
        ylp = ylp - yl
        yl = 0
    if zl<0:
        zlp = zlp - zl
        zl = 0
    try:
        res[xl:xh,yl:yh,zl:zh] = grid[xlp:xhp,ylp:yhp,zlp:zhp]
    except Exception as e:
        print(e)
        import pdb
        pdb.set_trace()
                    
    return res

if sys.platform=='win32':
    dropbox_dir = r'C:\Users\simon\Dropbox'
else:
    dropbox_dir = r'/media/simon/Dropbox/Dropbox'

v = vessel_network_3d()

gfile = r'C:\Users\simon\Dropbox\simulated_network.am'
sg = list_to_amira(v)
sg.write(gfile)
print('Written graph to: {}'.format(gfile))

parent = [x.parent for x in v]
coords = np.asarray([[v[p].start,x.start] for x,p in zip(v,parent) if p>=0])
coords = np.asarray([coords[x] for x in np.arange(0,len(coords),2)])
xmn,xmx = np.min(coords[:,:,0]),np.max(coords[:,:,0])
ymn,ymx = np.min(coords[:,:,1]),np.max(coords[:,:,1])
zmn,zmx = np.min(coords[:,:,2]),np.max(coords[:,:,2])

ms = 512 # matrix size
dims = [ms,ms,ms] # voxels
extent = [np.min([xmn,ymn,zmn]),np.max([xmx,ymx,zmx])] # um
extent_um = np.asarray([[extent[0],extent[1]],[extent[0],extent[1]],[extent[0],extent[1]]])
vsize = (extent_um[1]-extent_um[0])/float(ms) # um
grid = np.zeros(dims)

scale = np.eye(4)
tr = np.eye(4)
tr[0,3] = extent_um[0,0]
tr[1,3] = extent_um[1,0]
tr[2,3] = extent_um[2,0]
scale = np.eye(4)
for i in range(3):
    scale[i,i] = (extent_um[i,1]-extent_um[i,0])/dims[i]
trans = np.matmul(tr,scale)

#import pdb
#pdb.set_trace()

# Add vessels
vessel_val = 2
if True:
    print('Adding vessels to grid...')
    ofile = os.path.join(dropbox_dir,'simulated_network_nobg.nii')
    for i,e in enumerate(coords):
        path = v[i].path
        for j in range(1,len(path)):
            rad_c = v[i].r[j]            
            dif = path[j]-path[j-1]
            l = np.linalg.norm(dif)
        
            print('Coordinate {} of {}, radius {}'.format(i,len(coords),rad_c))
            R = vector_align(dif,[0.,0.,-1.])
            Rh = np.eye(4)
            Rh[0:3,0:3] = R
            scale, shear, angles, _, persp = decompose_matrix(Rh)
            angles = [np.rad2deg(a) for a in angles]
            
            com = np.mean(path[j-1:j+1,:],axis=0)
            print('L={},rot={},COM={},scale={}'.format(l,angles,com,scale))   
            ng = vessel_grid_fast(extent=extent_um,dims=dims,R=R,r=rad_c,c=com,l=l,rot=None,inside=vessel_val,outside=0)
            if np.max(ng)<=0:
                #import pdb
                #pdb.set_trace()
                pass
            grid += ng
    
        if (i % 10)==9:
            img = nib.Nifti1Image(grid,trans)
            nib.save(img,ofile)
     
    ofile = r'C:\Users\simon\Dropbox\simulated_network_nobg.nii'       
    img = nib.Nifti1Image(grid,trans)
    nib.save(img,ofile)
    print('Saved to: {}'.format(ofile)) 

if False:
    print('adding tissue to grid...')
    # Create 'tissue' extending as far as convex hull
    from scipy.spatial import ConvexHull
    
    def in_hull(pnt, hull):
        """
        Test if points in `p` are in `hull`
    
        `p` should be a `NxK` coordinates of `N` points in `K` dimensions
        `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
        coordinates of `M` points in `K`dimensions for which Delaunay triangulation
        will be computed
        """
        #from scipy.spatial import Delaunay
        #if not isinstance(hull,Delaunay):
        #    hull = Delaunay(hull)
        new_hull = ConvexHull(np.vstack((hull.points, pnt)))
        if np.array_equal(new_hull.simplices, hull.simplices): 
            return True
        return False
        #return hull.find_simplex(p)>=0
        
    img = nib.Nifti1Image(grid,trans)
    nib.save(img,ofile)
    
    tissue = np.zeros(dims)
    
    #import pdb
    #pdb.set_trace()
    if False:
        c_pix = np.asarray(np.where(grid>0))
        c_pix = np.transpose(c_pix, (1,0))
        conv_hull = ConvexHull(c_pix)
        
        conv_hull_coords = np.asarray([[conv_hull.points[v,0],conv_hull.points[v,1],conv_hull.points[v,2]] for v in conv_hull.vertices],dtype='int')
        mn_ch = np.min(conv_hull_coords,axis=0)
        mx_ch = np.max(conv_hull_coords,axis=0)
        count = 0
        for i in range(mn_ch[0],mx_ch[0]):
            for j in range(mn_ch[1],mx_ch[1]):
                for k in range(mn_ch[2],mx_ch[2]):
                    if in_hull([i,j,k],conv_hull):
                        print('Inside! {} {} {}'.format(i,j,k))
                        tissue[i,j,k] = 1
        #for i,cc in enumerate(inside):
        #    tissue[cc[0],cc[1],cc[2]] = 1
                        if (count % 100)==99:
                            img = nib.Nifti1Image(tissue+grid,trans)
                            nib.save(img,ofile)
                            print('Progress image written to {}'.format(ofile))
                        count += 1
    if False:
        c_pix = np.asarray(np.where(grid>0))
        c_pix = np.transpose(c_pix, (1,0))
        mn = np.min(c_pix,axis=0)
        mx = np.max(c_pix,axis=0)
        tissue[mn[0]:mx[0],mn[1]:mx[1],mn[2]:mx[2]] = 1
        
    grid += tissue
    
    ofile = os.path.join(dropbox_dir,'simulated_network.nii')
    img = nib.Nifti1Image(grid,trans)
    nib.save(img,ofile)
    print('Saved to: {}'.format(ofile))   
    
#import winsound
#duration = 1000  # millisecond
#freq = 440  # Hz
#winsound.Beep(freq, duration)