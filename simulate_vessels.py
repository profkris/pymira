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
from scipy.spatial import ConvexHull

class VesselSegment(object):
    
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
 
class VesselNetwork(object):

    def __init__(self,grow=True,**kwargs):
        if grow:
            self.network = self.create_network(**kwargs)
        else:
            self.network = None
        
    def asCartesian(self,rthetaphi):
        # Convert 3D polar to Cartesian coordinates
        r = rthetaphi[0]
        theta = rthetaphi[1] # degrees
        phi = rthetaphi[2] # degrees
        x = r * np.sin( theta ) * np.cos( phi )
        y = r * np.sin( theta ) * np.sin( phi )
        z = r * np.cos( theta )
        return [x,y,z]
        
    def create_network(self,nlevel=9,nb=2,gamma=3.,beta=1.,f=0.9,nf=0.9,theta=27.,phi=10., 
                       sigmaL=0.02,sigmaD=0.02,sigmaT=0.02,sigmaP=10.,tortFactor=0.5,verbose=False):
        """
        Creates a 3D blood vessel-like network in graph format
        
        Input parameters:
          nlevel = number of levels in vessel hierarchy
          nb = number of branches per vessel (only nb=2 for now!)
          gamma = power law exponent for vessel diameters: D0^gamma = D1^gamma + D2^gamma
          beta = symmetry variable: beta = D1 / D2 (0<beta<1)
          f = ratio of parent length to daughter length: l_daughter = f x l_parent
          nf = Fraction of length at which node occurs
          theta, phi = vessel branching angles. 
              25.5<theta<28.5 for kidney; 25<theta<140 for tumour
          sigma(L,D,T,P) = standard deviation of random variables (length,diameter,theta,phi)
          tortFactor = tortuosity factor used for constrained random walk to connect nodes.
                       If zero or negative, a straight line is used.
        """
        
        # Number of vessels to model
        nvessel = np.sum([nb**i for i in range(nlevel)])
        if verbose:
            print('No. vessels: {}'.format(nvessel))
      
        # Define a structure to hold network information
        v = [VesselSegment() for i in range(nvessel)]

        # Branching angles
        thetaR = np.deg2rad(theta)
        phiR = np.deg2rad(phi)

        # Define objects for holding growth process information
        gp = GrowthProcess('normal',gamma=gamma,beta=beta,f=f,nf=nf,theta=thetaR,phi=phiR, \
            sigmaL=sigmaL,sigmaD=sigmaD,sigmaT=sigmaT,sigmaP=sigmaP)
            
        seg_count = 0 # Vessel segment counter
        node_count = 0 # Node counter
      
        # Populate first segment structure, then begin growth 
        v[0].occ = True
        v[0].l = 4500. #5000 #um
        v[0].r = 500. #um
        v[0].theta = 0. #rad
        v[0].phi = 0. #rad
        v[0].start = [0.,0.,0.] #x,y,z (um)
        v[0].startIndex = 0
        v[0].end = np.asarray(v[0].start) + self.asCartesian([v[0].l,v[0].theta,v[0].phi])
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
        nparent = 1
        pSubs = [0]
        level = 0
        level = level + 1

        # Step counters to reflect 1 new segment and two new nodes
        seg_count += 1
        node_count += 2

        # Iteratively add segments to initial segment
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
                #    plot_vessel_network(v)
                #    endloop =True
                    break
                nxt = nxt[0]
           
                for j in range(len(dnew)):
                    # Generate random numbers
                    rndmL = np.random.uniform(low=-gp.sigmaL,high=gp.sigmaL)
                    rndmD = np.random.uniform(low=-gp.sigmaD,high=gp.sigmaD)
                    rndmT = np.random.uniform(low=-gp.sigmaT,high=gp.sigmaT)
                    rndmP = np.random.uniform(low=-gp.sigmaP,high=gp.sigmaP)      
            
                    newV = v[nxt[j]] # Grab and populate new segment structure
                    newV.proc = curV.proc        
                    newl = curV.l * f
                    newV.l = newl + (newl*rndmL) #um
                    newV.r = dnew[j] + (dnew[j]*rndmD)
                    if j % 2==0: # Every second branch (assumes nb=2!)
                        newV.theta = curV.theta + gp.theta + (gp.theta*rndmT)
                    else: 
                        newV.theta = curV.theta - gp.theta + (gp.theta*rndmT)
                    if j % 2==0:
                        newV.phi = curV.phi + gp.phi + (gp.phi*rndmP)
                    else:
                        newV.phi = curV.phi - gp.phi + (gp.phi*rndmP)
                    newV.start = node #x,y (um)
                    newV.end = node + self.asCartesian([newV.l,newV.theta,newV.phi])
                    newV.term = True
                    newV.parent = curV.index
                    newV.level = level
                    newV.index = seg_count
                    newV.startIndex = curV.endIndex
                    newV.endIndex = node_count
            
                    node_count += 1
                    ntry = 10 # If using tortured path, allow 10 attempts to connect nodes
                    count1 = 0
                    while True: # Try and connect to node using constrained random walk (if tortFactor>=0)
                        if tortFactor<=0:
                            path = np.asarray([newV.start,newV.end])
                        else:
                            tort = np.linalg.norm(np.asarray(newV.start)-np.asarray(newV.end)) * tortFactor
                            path = tortured_path.tortured_path(newV.start,newV.end,tortuosity=tort) # np.asarray([curV.start,node])
                        if path is not None:
                            if np.any(path[0]!=newV.start) or np.any(path[-1]!=newV.end): # Check start and end points are correct
                                import pdb
                                pdb.set_trace()
                            break
                        else:
                            count1 += 1
                            if count1 >= ntry: # If number of allowed attempts exceeded, debug...
                                import pdb
                                pdb.set_trace()
                    newV.path = path
                    if not np.isscalar(newV.r):
                        import pdb
                        pdb.set_trace()
                    newV.r = np.zeros(path.shape[0]) + newV.r
                    newV.occ = True
                    v[nxt[j]] = newV
                    seg_count += 1

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
        
    def list_to_amira(self,v):
    
        """
        Convert a vessel graph to Amira spatial graph format 
        """
        from pymira import spatialgraph
        sg = spatialgraph.SpatialGraph(initialise=True)
        sg.add_field(name='Radius',marker='@5',
                              definition='POINT',type='float',
                              nelements=1,nentries=[0],data=None)
                              
        #nnode = 0
        #nedge = 0
        #npoint = 0
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
        
        sg.definitions[0]['size'] = [nnodes]
        sg.definitions[1]['size'] = [len(conns)]
        sg.definitions[2]['size'] = [np.sum(npoints)]
        return sg            

class EmbedInGrid(object):

    def coord_to_pix(self,coord,extent_um,dims,clip=True):
        """ Convert 3D spatial coordinates into 3D pixel indices """
        if clip:
            return np.asarray([np.clip(np.round((coord[i]-extent_um[i,0])*dims[i]/(extent_um[i,1]-extent_um[i,0])),0,dims[i]) for i in [0,1,2]],dtype='int')
        else:
            return np.asarray([(coord[i]-extent_um[i,0])*dims[i]/(extent_um[i,1]-extent_um[i,0]) for i in [0,1,2]],dtype='int')
    
    def vector_align(self,v1,v2):
        """
        Calculate transform to align two 3D vectors
        """
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

    def rotate_3d(self,vol,R,centre=[0,0,0]):
        """
        Perform rotation of a 3D array about an arbitrary axis
        """
        
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
                
    def embed_cylinder(self,extent=None,dims=None,r=10.,c=np.asarray([0.,0.,0.]),
                       R=None,rot=[0.,0.,0.],l=10.,inside=1,outside=0,verbose=False):
        
        """ 
        Embed a cylinder inside a 3D array
        extent: bounding box for matrix (um)
        dims: grid dimensions
        r: radius (um)
        c: centre (um)
        l: length (um)
        inside: pixel value for inside vessel
        outside: pixel value for outside vessel
        """

        r_pix = [np.int(np.round(r*dims[i]/(extent[i,1]-extent[i,0]))) for i in range(3)]
        l_pix = [np.clip(np.int(np.round(l*dims[i]/(extent[i,1]-extent[i,0]))),0,dims[i]) for i in range(3)]
        if verbose:
            print('Lpix {}, Rpix {}'.format(l_pix,r_pix))
        
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
        grid = self.rotate_3d(grid,R,centre=offset)
        
        # Shift to correct position
        com_pix = self.coord_to_pix(c,extent,dims)
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
        
    def create_grid_transform(self,extent_um,dims):
        # Create transformation matrix for embedding grid
        scale = np.eye(4)
        tr = np.eye(4)
        tr[0,3] = extent_um[0,0]
        tr[1,3] = extent_um[1,0]
        tr[2,3] = extent_um[2,0]
        scale = np.eye(4)
        for i in range(3):
            scale[i,i] = (extent_um[i,1]-extent_um[i,0])/dims[i]
        trans = np.matmul(tr,scale)
        return trans
        
    def write(self,grid,ofile,trans=np.eye(4)):
        img = nib.Nifti1Image(grid,trans)
        nib.save(img,ofile)
        
    def embed_vessels(self,v,extent_um=None,ms=512,inside=2,outside=0,
                      ofile=None,save_every=10,verbose=False):
    
        """
        Embed a vessel network (graph format) in a 3D array
        Input args:
            ms: matrix size
            inside: value inside vessel
        """
    
        # Find parent nodes
        parent = [x.parent for x in v]
        coords = np.asarray([[v[p].start,x.start] for x,p in zip(v,parent) if p>=0])
        coords = np.asarray([coords[x] for x in np.arange(0,len(coords),2)])
        
        # Find network extent
        xmn,xmx = np.min(coords[:,:,0]),np.max(coords[:,:,0])
        ymn,ymx = np.min(coords[:,:,1]),np.max(coords[:,:,1])
        zmn,zmx = np.min(coords[:,:,2]),np.max(coords[:,:,2])

        dims = [ms,ms,ms] # voxels
        extent = [np.min([xmn,ymn,zmn]),np.max([xmx,ymx,zmx])] # um
        extent_um = np.asarray([[extent[0],extent[1]],[extent[0],extent[1]],[extent[0],extent[1]]])
        #vsize = (extent_um[1]-extent_um[0])/float(ms) # um
        
        # Create grid
        grid = np.zeros(dims)
        trans = self.create_grid_transform(extent_um,dims)

        # Add vessels
        if verbose:
            print('Adding vessels to grid...')
        
        for i,e in enumerate(coords):
            path = v[i].path
            for j in range(1,len(path)):
                rad_c = v[i].r[j]            
                dif = path[j]-path[j-1]
                l = np.linalg.norm(dif)
                
                if verbose:
                    print('Coordinate {} of {}, radius {}'.format(i,len(coords),rad_c))
                R = self.vector_align(dif,[0.,0.,-1.])
                #Rh = np.eye(4)
                #Rh[0:3,0:3] = R
                #scale, shear, angles, _, persp = decompose_matrix(Rh)
                #angles = [np.rad2deg(a) for a in angles]
                    
                com = np.mean(path[j-1:j+1,:],axis=0)
                #print('L={},rot={},COM={},scale={}'.format(l,angles,com,scale))   
                # Embed current segment
                ng = self.embed_cylinder(extent=extent_um,dims=dims,R=R,r=rad_c,c=com,l=l,rot=None,inside=inside,outside=outside)
                # Add segment to main grid
                grid += ng
            
            if save_every>0 and (i % save_every)==save_every-1 and ofile is not None:
                self.write(grid,ofile,trans=trans)
         
        if ofile is not None:
            self.write(grid,ofile,trans=trans)
            
        return grid,extent_um,trans
    
    def in_hull(self,pnt, hull):
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
        
    def add_convext_hull_tissue(self,grid):
        """
         Create 'tissue' extending as far as convex hull
         NOT WORKING!
        """
        print('adding tissue to grid...')
    
        #img = nib.Nifti1Image(grid,trans)
        #nib.save(img,ofile)
    
        dims = grid.shape
        tissue = np.zeros(dims)

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
                    if self.in_hull([i,j,k],conv_hull):
                        print('Inside! {} {} {}'.format(i,j,k))
                        tissue[i,j,k] = 1
        #for i,cc in enumerate(inside):
        #    tissue[cc[0],cc[1],cc[2]] = 1
                        if (count % 100)==99:
                            img = nib.Nifti1Image(tissue+grid,trans)
                            nib.save(img,ofile)
                            print('Progress image written to {}'.format(ofile))
                        count += 1
        grid += tissue
        return grid
                        
    def add_tissue(self,grid,inside=1):
        
        dims = grid.shape
        tissue = np.zeros(dims)

        c_pix = np.asarray(np.where(grid>0))
        c_pix = np.transpose(c_pix, (1,0))
        mn = np.min(c_pix,axis=0)
        mx = np.max(c_pix,axis=0)
        tissue[mn[0]:mx[0],mn[1]:mx[1],mn[2]:mx[2]] = inside
        
        grid += tissue
        return grid  

def main():
    if sys.platform=='win32':
        dropbox_dir = r'C:\Users\simon\Dropbox'
    else:
        dropbox_dir = r'/media/simon/Dropbox/Dropbox'
        
    opath = os.path.join(dropbox_dir,'test_sim')
    if not os.path.exists(opath):
        os.makedirs(opath)

    import pdb
    pdb.set_trace()
    vn = VesselNetwork(grow=True)
    v = vn.network

    gfile = os.path.join(opath,'simulated_network_{}.am'.format(1))
    sg = vn.list_to_amira(v)
    sg.write(gfile)
    print('Written graph to: {}'.format(gfile))
    
    emb = EmbedInGrid()
    ofile = os.path.join(dropbox_dir,'simulated_network_nobg_{}.nii'.format(1))     
    grid,extent_um,trans = emb.embed_vessels(v,ms=512,inside=2,outside=0)
    grid = emb.add_tissue(grid)
    ofile = os.path.join(dropbox_dir,'simulated_network.nii')
    emb.write(grid,ofile,trans=trans)
    print('Saved to: {}'.format(ofile)) 

if __name__ == "__main__":
    main()