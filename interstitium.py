# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 15:15:58 2017

@author: simon
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 03 09:38:47 2017

@author: simon
"""

from tqdm import tqdm # progress bar

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy.signal
from scipy.ndimage import filters
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def impulse(t,delay):

    length = 20. #us
    conc = np.zeros(t.shape[0])
    conc[t>=delay] = 1.
    conc[t>(delay+length)] = 0.
    return conc
    
def parker(t,delay):
    
    # PARKER (average arterial input function from human population)
    import numpy as np
    try:
        a1 = 0.833
        a2 = 0.336
        t1 = 0.171
        t2 = 0.364
        sigma1 = 0.055 
        sigma2 = 0.134
        alpha = 1.064
        beta = 0.166
        s = 37.772
        tau = 0.482
        tMin = (t-delay) / 60.
        r1 = (a1 / (sigma1*np.sqrt(2.*np.pi))) * np.exp(-np.square(tMin-t1)/(2.*np.square(sigma1)))
        r2 = (a2 / (sigma2*np.sqrt(2.*np.pi))) * np.exp(-np.square(tMin-t2)/(2.*np.square(sigma2)))
        wo = alpha * np.exp(-beta*tMin) / (1. + np.exp(-s*(tMin-tau)))
        conc = r1 + r2 + wo
        
        conc[tMin<0] = 0.
        conc[conc<0] = 0.
        conc[~np.isfinite(conc)] = 0.
        
    except Exception,e:
        print(e)
        import pdb
        pdb.set_trace()
        
    return conc

class Interstitium(object):
    
    def __init__(self):
        self.grid = None
        self.count = None
        self.pixSize = [50.,50.,20.] #um
        self.dx = None
        self.dy = None
        self.dz = None
        self.dt = None
        
    def align_vector_rotation(self,A,B):
        
        """ To align vector A to B:
        B = np.transpose(U*A.T))
        """
        
        #mat = np.asmatrix
        norm = np.linalg.norm
        
        A = A/norm(A)
        B = B/norm(B)
        #print('A={}'.format(A))
        #print('B={}'.format(B))
        
        if np.all(A==B):
            #print('All equal')
            return np.identity(3)
        elif np.all(A==-B):
            #print('All negative')
            return -np.identity(3)
        
        G = np.asarray(
             [ [np.dot(A,B.T), (-norm(np.cross(A,B))), 0.],
               [norm(np.cross(A,B)), np.dot(A,B.T),  0.],
               [0., 0., 1.],
             ])
        G = np.asmatrix(G)
    
        Fi = [ A , (B - (np.dot(A,B.T))*A) / norm((B - (np.dot(A,B.T))*A)) , np.cross(B,A) ]
        Fi = np.squeeze(Fi)
        Fi = np.nan_to_num(Fi)
        Fi = np.asmatrix(Fi)
        try:
            U = Fi * G * np.linalg.inv(Fi)
        except Exception,e:
            print(e)
            import pdb
            pdb.set_trace()
        
        return U
    
    def sphere_coordinates(self,radius,centre,n):
        phi = np.linspace(0,2*np.pi,n)
        theta = np.linspace(0,2*np.pi,n)
        coords = []
        for th in theta:
            for ph in phi:
                coords.append( [radius*np.sin(ph)*np.cos(th) + centre[0],
                                radius*np.sin(ph)*np.sin(th) + centre[1],
                                radius*np.cos(ph) + centre[2] ])
        return coords
        
    def cylinder_coordinates(self,radius,strt,fin,n):
        
        length = np.linalg.norm(strt-fin)
        orientation = (strt-fin)/length
        I = np.linspace(0, length, n)
        #x = np.linspace(0,length,num=n)
        #y = np.linspace(strt[1],fin[1],num=n)
        #z = np.linspace(strt[2],fin[2],num=n)
        
        A = np.linspace(0, 360., num=n, endpoint=False)/180*np.pi
        #import pdb
        #pdb.set_trace()
        coords = []
        #for r in radius:
        if True:
            r = radius
            if r>0:
                X = np.asarray(r) * np.cos(A)
                Y = np.asarray(r) * np.sin(A)
                #Tile/repeat indices so all unique pairs are present
                px = np.repeat(I, n)
                py = np.tile(X, n)
                pz = np.tile(Y, n)
            else:
                px = I
                py = np.repeat(0.,n)
                pz = np.repeat(0.,n)
            coords = np.vstack(( px, py, pz )).T
        
        origin, xaxis, yaxis, zaxis = [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
        U = self.align_vector_rotation([1,0,0],orientation)
        coords_p = np.dot(U,np.transpose(coords)).T
        coords_p = np.asarray(coords_p)
        
        #Shift to center
        center = np.mean([strt,fin],axis=0)
        shift = np.array(center) - np.mean(coords_p, axis=0)
        coords_p += shift
        
        return np.asarray(coords_p,dtype='float')
        
    def find_nearest(self,array,value):
        idx = (np.abs(array-value)).argmin()
        return idx #array[idx]
        
    def radial_diffusion(self,coords,c_v,D,ktrans,k,h,nr,nt,time):
        
        lam = D*k/np.square(h) # + np.power(dr,-2) + np.power(dr,-2))
        #print('Lambda: {}'.format(lam))
        A = np.zeros((nr,nr))
        for i in xrange(nr):
            A[i,i] = 1 - 2.*lam
            if i>0:
                A[i-1,i] = lam
            if i<nr-1:
                A[i+1,i] = lam
        
        c_i = np.zeros((nr,nt))
        R = np.zeros((nr,nt))
        
        #inds = np.where(r<=radius)
        inds = 0
        filt = ktrans*np.exp(-ktrans*time)
        conv = np.convolve(c_v,filt) # Convolve impulse response with AIF
        conv = conv[0:len(time)]*(time[1]-time[0])
        
        conv = c_v
        
        R[inds,0:len(c_v)] = conv
        
        for j in xrange(nt-1):
            c_i[:,j+1] = np.dot(A,c_i[:,j]) + k*R[:,j] #(np.dot(B,R[:,j]))
        
        return c_i
        
    def interstitial_diffusion(self,nodes,index,conc,time,plot_conc=False):
        
        nnode = len(nodes)
        nodes = np.asarray(nodes)
        mnx,mxx = np.min(nodes[:,0]),np.max(nodes[:,0])
        mny,mxy = np.min(nodes[:,1]),np.max(nodes[:,1])
        mnz,mxz = np.min(nodes[:,2]),np.max(nodes[:,2])
    
        #radius = [10.,10.,10.] #um
        #delay = [0.]*nnode #s
        ktrans = [0.2]*nnode #/min
        D = [500.]*nnode # um2/s
        
        ktrans = [k/60. for k in ktrans] #/s
        nr = 100
        max_r = 1000. #um #dr*nr
        dr = max_r / float(nr)
        r = np.linspace(0,max_r,num=nr)
        #r = np.logspace(-10,np.log10(max_r),num=nr,endpoint=True)
        
        #embedDims = [[0,40],[0,40],[5,15]] #um
        pad = max_r / 5.
        embedDims = [ [mnx-pad,mxx+pad],
                      [mny-pad,mxy+pad],
                      [mnz-pad,mxz+pad] ]
        self.embedDims = embedDims

        dt = time[1] - time[0]
        max_time = time[-1]
        nt = len(time)
        print('Max TIME: {} s, {} min'.format(max_time,max_time/60.))
        print('dt: {} s'.format(dt))
        print('nt: {} s'.format(nt))
        print('Max RADIUS: {} um'.format(max_r))
        print('dr: {} um'.format(dr))
        print('nr: {} um'.format(nr))
        print('D: {} s2/um'.format(D[0]))
        print('Ktrans: {} /s'.format(ktrans[0]))
        print('Embedding dims: {} um'.format(embedDims))
        
        k = dt
        h = dr
        
        self.dt = dt
        
        print('k<=h2/2D   k (dt) = {}, h (dr) = {}, h2/2D = {}'.format(k,h,np.square(h)/(2.*D[0])))
        assert all(k<=np.square(h)/(2.*Dc) for Dc in D)
        
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #plt.gca().set_aspect('equal', adjustable='box')
        
        # Regrid
        ntg = nt
        ng = [np.int(np.ceil((embedDims[i][1]-embedDims[i][0])/self.pixSize[i])) for i in range(3)]
        ng[2] = 1
        self.ng = ng

        self.grid = np.zeros([ntg,ng[0],ng[1]])#,ng[2]])
        self.count = np.zeros([ntg,ng[0],ng[1]])#,ng[2]])
        dx = (embedDims[0][1]-embedDims[0][0]) / float(ng[0])
        dy = (embedDims[1][1]-embedDims[1][0]) / float(ng[1])
        dz = (embedDims[2][1]-embedDims[2][0]) / float(ng[2])
        self.dx,self.dy,self.dz = dx,dy,dz
        xg = np.linspace(embedDims[0][0],embedDims[0][1],num=ng[0],dtype='float')
        yg = np.linspace(embedDims[1][0],embedDims[1][1],num=ng[1],dtype='float')
        zg = np.linspace(embedDims[2][0],embedDims[2][1],num=ng[2],dtype='float')
        tg = np.linspace(0,max_time,num=ntg,dtype='float')
        tgInd = np.linspace(0,nt-1,num=ntg,dtype='int')
        
        print('Grid dx,dy,dz: {} {} {}'.format(dx,dy,dz))
        
        #import pdb
        #pdb.set_trace()
        
        cm = plt.get_cmap('jet')
        mn = 0.
        mx = 1e-10
        cNorm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        
        pbar = tqdm(total=nnode)
        
        for ni,curNode in enumerate(nodes):
            
            pbar.update(1)
            
            if ni<nnode-1 and index[ni]==index[ni+1]: # for cylinder only...
            
                #ve = 0.1            
                c_v = conc[:,ni] # impulse(time,delay[ni])
                c_i = self.radial_diffusion(curNode,c_v,D[ni],ktrans[ni],k,h,nr,nt,time)
                
                if plot_conc:
                    fig = plt.figure()
                    plt.plot(time,c_v)
                    for i in xrange(np.max([nr,100])):
                        plt.plot(time,c_i[i,:])
                # Regrid
                c_i = c_i[:,tgInd] 
                
                xn,yn,zn = (self.find_nearest(xg,curNode[0]),
                           self.find_nearest(yg,curNode[1]),
                           self.find_nearest(zg,curNode[2]))
            
                for ri,radius in enumerate(r):
                    #sph = self.sphere_coordinates(radius,curNode,10)
                    sph = self.cylinder_coordinates(radius,curNode,nodes[ni+1],10)
                    xs = [x[0] for x in sph]
                    ys = [x[1] for x in sph]
                    zs = [x[2] for x in sph]
                    for i,v in enumerate(xs):
                        #print(xs[i])
                        inGrid = all([xs[i]>=embedDims[0][0], xs[i]<=embedDims[0][1],
                                      ys[i]>=embedDims[1][0], ys[i]<=embedDims[1][1],
                                      zs[i]>=embedDims[2][0], zs[i]<=embedDims[2][1]] )
                        if inGrid:
                            xsc,ysc,zsc = (self.find_nearest(xg,xs[i]),
                                           self.find_nearest(yg,ys[i]),
                                           self.find_nearest(zg,zs[i]) )
                            #if xsc!=xn or ysc!=yn:
                            #    print('Outside node! {} {}'.format(ri,c_i[ri,:]))
                            self.grid[:,xsc,ysc] += c_i[ri,:]
                            self.count[:,xsc,ysc] += 1

        for i in xrange(nt):
            self.grid[i] = scipy.ndimage.filters.median_filter(self.grid[i],size=5)
        
        pbar.close()
        #self.normalise()
        
    def normalise(self):
        self.grid[self.count>0] = self.grid[self.count>0] / self.count[self.count>0]
        
    def save_grid(self,path):
                
        print('Max C_i: {} mM'.format(np.max(self.grid)))
            
        import nibabel as nib
        tmp = np.swapaxes(self.grid,0,1)
        tmp = np.swapaxes(tmp,1,2)
        img = nib.Nifti1Image(tmp[:,:,None,:],affine=np.eye(4))
        hdr = img.header
        hdr['pixdim'][1] = self.dx
        hdr['pixdim'][2] = self.dy
        hdr['pixdim'][3] = self.dz
        hdr['pixdim'][4] = self.dt
        ofile = os.path.join(path,'interstitial.nii')
        print('Saving to {}'.format(ofile))
        nib.save(img,ofile)
    
    def display_grid(self,last=True):
        #scalarMap.set_array([mn,mx])
        #fig.colorbar(scalarMap)
        if last:
            plt.figure()
            im = self.grid[-1,:,:]
            plt.imshow(im,vmin=0,vmax=0.002)
        else:
            for i in xrange(0,self.grid.shape[0],10):
                plt.figure()
                im = self.grid[i,:,:]
                plt.imshow(im,vmin=0,vmax=0.01)
        
plt.close('all')

from pymira import spatialgraph
dir_ = 'C:\\Users\\simon\\Dropbox\\Mesentery\\'
f = dir_ + 'ct_output.am'
graph = spatialgraph.SpatialGraph()
print('Reading graph...')
graph.read(f)
print('Graph read')

points = graph.get_data('EdgePointCoordinates')
npoints = points.shape[0]
nEdgePoint = graph.get_data('NumEdgePoints')
edgePointIndex = np.zeros(npoints,dtype='int')

offset = 0
edgeCount = 0
for npCur in nEdgePoint:
    edgePointIndex[offset:offset+npCur] = edgeCount
    edgeCount += 1
    offset += npCur
    
assert offset==npoints

concFieldInds = [i for i,x in enumerate(graph.fieldNames) if 'Concentration' in x]
concFields = [graph.fieldNames[i] for i in concFieldInds]
time = np.asarray([float(x.replace('Concentration_','')) for x in concFields])

nt = len(time)
npoint = points.shape[0]

conc = np.zeros((nt,npoint),dtype='float')
for ci,concField in enumerate(concFieldInds):
    field = graph.fields[concField]
    if 'data' in field:
        conc[ci,:] = field['data']
    else:
        conc[ci,:] = conc[ci-1,:]
        print('Data missing: {}'.format(concFields[ci]))

#nodes = [[10.,10.,10.],
#         [10.,20.,10.],
#         [30.,30.,10.]]

#import pdb
#pdb.set_trace()

jump = 1
lim = -1#jump*50
time = time[0:lim:jump]
conc = conc[0:lim:jump,:]

#conc = conc[:,0:100]
#points = points[0:100,:]
#edgePointIndex = edgePointIndex[0:100]

inter = Interstitium()

try:
    inter.interstitial_diffusion(points,edgePointIndex,conc,time)
except KeyboardInterrupt:
    print('Keyboard interrupt')
    #inter.normalise()
    #pass
#import pdb
#pdb.set_trace()
#inter.display_grid(last=False)
inter.save_grid(dir_)