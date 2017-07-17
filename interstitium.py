# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 15:15:58 2017

@author: simon

Finite element analysis of interstitial diffusion and vascular leakage
Links to injectagent.py

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

class Interstitium(object):
    
    def __init__(self,paramSet=None):
        if paramSet is None:
            from pymira import inject_agent
            paramSet = inject_agent.ParameterSet()
        self.grid = None
        self.count = None
        self.pixSize = paramSet.pixSize #um
        self.dx = None
        self.dy = None
        self.dz = None
        self.dt = None
                        
        self.ktrans = paramSet.ktrans #0.00001 #/min
        self.ef = paramSet.ef #  not used
        #self.D = 1e-7 * 1e8 # 500. #um2/s
        self.D = paramSet.D #m2/s Gd-DTPA (https://books.google.es/books?id=6fZGf8ave3wC&pg=PA343&lpg=PA343&dq=diffusion+coefficient+of+gd-dtpa&source=bl&ots=Ceg432CWar&sig=4PuxViFn9lL7pwOAkFVGwtHRe4M&hl=en&sa=X&ved=0ahUKEwjs1O-Z-NPTAhVJShQKHa6PBKQQ6AEIODAD#v=onepage&q=diffusion%20coefficient%20of%20gd-dtpa&f=false)
        
        self.feNSample = paramSet.feNSample
        
    def align_vector_rotation(self,A,B):
        
        """ Calculates a rotation matrix to align A to B
        B = np.transpose(U*A.T))
        """

        norm = np.linalg.norm        
        A = A/norm(A)
        B = B/norm(B)
        
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
        
    def circle_coordinates(self,radius,strt,fin,n):
        
        length = np.linalg.norm(strt-fin)
        orientation = (strt-fin)/length
        I = 0
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
                px = np.repeat(0.,n)
                py = np.tile(X, n)
                pz = np.tile(Y, n)
                coords = np.vstack(( px, py, pz )).T
            else:
                px = np.repeat(0.,n)
                py = np.repeat(0.,n)
                pz = np.repeat(0.,n)
                #coords = [0.,0.,0.]
                coords = np.vstack(( px, py, pz )).T
        
        #origin, xaxis, yaxis, zaxis = [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
        U = self.align_vector_rotation([1,0,0],orientation)
        coords_p = np.dot(U,np.transpose(coords)).T
        coords_p = np.asarray(coords_p)
        
        #Shift to center
        center = np.mean([strt],axis=0)
        shift = np.array(center) - np.mean(coords_p, axis=0)
        coords_p += shift
        
        return np.asarray(coords_p,dtype='float')
        
    def find_nearest(self,array,value):
        idx = (np.abs(array-value)).argmin()
        return idx #array[idx]
        
    def radial_diffusion(self,coords,c_v,D,ktrans,ef,k,h,nr,nt,time):
        
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
        c_v_out = c_v
        R = np.zeros((nr,nt))
        
        #inds = np.where(r<=radius)
        inds = 0
        #filt = ef*ktrans*np.exp(-ktrans*time)
        #conv = np.convolve(c_v,filt) # Convolve impulse response with AIF
        #conv = conv[0:len(time)]*(time[1]-time[0])
        
        #conv = c_v
        
        #R[inds,0:len(c_v)] = c_v - c_i[0,j]
        
        #http://iopscience.iop.org/article/10.1088/0031-9155/59/17/5175
        #mu = ktrans
        #sb_vi_av = 333.e-4 #/um
        #ubi = 10.e-4 #um/s
        #ktrans = sb_vi_av * ubi
        M = (c_i.shape[0]/float(self.feNSample)) * ktrans #* (vessel_surfArea / interstitial_volume)
        
        for j in xrange(nt-1):
            c_i[:,j+1] = np.dot(A,c_i[:,j]) + M*(c_v[j]-c_i[:,j]) #(np.dot(B,R[:,j]))
            c_v_out[j+1] += np.sum(M*(c_i[0,j]-c_v[j]))
        c_v_out = np.clip(c_v_out,0.,1e100)
            
        #c_v_out = (1.-ef) * c_v
        #c_v_out = c_v
        
        return c_i,c_v_out
        
    def set_grid_dimensions(self,nodes,time,ktrans=None,D=None,verbose=False,embed_dims=None,grid_dims=None):
        
        if ktrans is None:
            ktrans = self.ktrans
        if D is None:
            D = self.D
        
        self.nr = 10
        self.max_r = 1000. #um #dr*nr
        self.dr = self.max_r / float(self.nr)

        if embed_dims is None:
            pad = self.max_r / 5.
            mnx,mxx = np.min(nodes[:,0]),np.max(nodes[:,0])
            mny,mxy = np.min(nodes[:,1]),np.max(nodes[:,1])
            mnz,mxz = np.min(nodes[:,2]),np.max(nodes[:,2])
            self.embedDims = [ [mnx-pad,mxx+pad],
                          [mny-pad,mxy+pad],
                          [mnz-pad,mxz+pad] ]
        else:
            self.embedDims = embed_dims

        self.dt = time[1] - time[0]
        self.max_time = time[-1]
        self.nt = len(time)
        
        k = self.dt
        h = self.dr
        
        test = k<=np.square(h)/(2.*D)
        
        if verbose or not test:
            print('Max TIME: {} s, {} min'.format(self.max_time,self.max_time/60.))
            print('dt: {} s'.format(self.dt))
            print('nt: {} s'.format(self.nt))
            print('Max RADIUS: {} um'.format(self.max_r))
            print('dr: {} um'.format(self.dr))
            print('nr: {} um'.format(self.nr))
            print('D: {} um2/s'.format(D))
            print('Ktrans: {} /s'.format(ktrans))
            print('Embedding dims: {} um'.format(self.embedDims))
            print('k<=h2/2D   k (dt) = {}, h (dr) = {}, h2/2D = {}'.format(k,h,np.square(h)/(2.*D)))        
        
        assert test
        
        # Regrid
        ntg = self.nt
        ng = [np.int(np.ceil((self.embedDims[i][1]-self.embedDims[i][0])/self.pixSize[i])) for i in range(3)]
        #print('Embedding matrix: {}'.format(ng))
        flatten_z = False
        if flatten_z:
            ng[2] = 1
            self.grid = np.zeros([ntg,ng[0],ng[1]])#,ng[2]])
            self.count = np.zeros([ntg,ng[0],ng[1]])#,ng[2]])
        else:
            self.grid = np.zeros([ntg,ng[0],ng[1],ng[2]])
            self.count = np.zeros([ntg,ng[0],ng[1],ng[2]])
        
        #self.ng = ng
        #xg = np.linspace(embedDims[0][0],embedDims[0][1],num=ng[0],dtype='float')
        #yg = np.linspace(embedDims[1][0],embedDims[1][1],num=ng[1],dtype='float')
        #zg = np.linspace(embedDims[2][0],embedDims[2][1],num=ng[2],dtype='float')
        #tg = np.linspace(0,max_time,num=ntg,dtype='float')
        #tgInd = np.linspace(0,nt-1,num=ntg,dtype='int')

        if grid_dims is None:
            self.ntg = self.nt
            self.ng = [np.int(np.ceil((self.embedDims[i][1]-self.embedDims[i][0])/self.pixSize[i])) for i in range(3)]
            self.grid_dims = [self.ntg,self.ng[0],self.ng[1],self.ng[2]]
        else:
            self.grid_dims = grid_dims
            self.ntg = grid_dims[0]
            self.ng = [grid_dims[1],grid_dims[2],grid_dims[3]]

        dx = (self.embedDims[0][1]-self.embedDims[0][0]) / float(self.ng[0])
        dy = (self.embedDims[1][1]-self.embedDims[1][0]) / float(self.ng[1])
        dz = (self.embedDims[2][1]-self.embedDims[2][0]) / float(self.ng[2])
        self.dx,self.dy,self.dz = dx,dy,dz
        
    def interstitial_diffusion(self,nodes,index,conc,time,plot_conc=False,set_grid_dims=True,progress=True,flatten_z=False,med_reg=False,store_results=True):
        
        """ Simulates vascular exchange and diffusion through interstitium
        """        
        
        # Get size of vessel segment
        nnode = len(nodes)
        nodes = np.asarray(nodes)
 
        if store_results:
            if set_grid_dims:       
                self.set_grid_dimensions(nodes,time,ktrans=self.ktrans,D=self.D)
                
            self.grid = np.zeros(self.grid_dims)
            self.count = np.zeros(self.grid_dims)
        
        # Create lists of model parameters corresponding to number of nodes in current vessel segment
        ktrans = [self.ktrans/60.]*nnode #/min
        ef = [self.ef]*nnode #/min
        D = [self.D]*nnode # um2/s
        
        # Set radial array for finite element calc
        r = np.linspace(0,self.max_r,num=self.nr)
        
        # Set dimensions of embedding lattice
        embedDims = self.embedDims
        xg = np.linspace(embedDims[0][0],embedDims[0][1],num=self.ng[0],dtype='float')
        yg = np.linspace(embedDims[1][0],embedDims[1][1],num=self.ng[1],dtype='float')
        zg = np.linspace(embedDims[2][0],embedDims[2][1],num=self.ng[2],dtype='float')
        #tg = np.linspace(0,self.max_time,num=self.ntg,dtype='float')
        tgInd = np.linspace(0,self.nt-1,num=self.ntg,dtype='int')
        
        if progress: # progress bar
            pbar = tqdm(total=nnode)
            
        # Initialise array for modified vascular component (following leakage into interstitium)
        conc_out = conc * 0.
        
        # Loop through each node in current vessel segment
        for ni,curNode in enumerate(nodes):
            if progress: # progress bar
                pbar.update(1)
            
            if ni<nnode-1 and index[ni]==index[ni+1]: # for cylinder only...
            #if True:

                c_v = conc[ni,:]
                c_i,c_v_out = self.radial_diffusion(curNode,c_v,D[ni],ktrans[ni],ef[ni],self.dt,self.dr,self.nr,self.nt,time)
                conc_out[ni,:] = c_v_out
                
                if store_results: # To bypass finite element bit (speeds up calculation and reduces memory)
                    # Regrid
                    c_i = c_i[:,tgInd] 
    
                    for ri,radius in enumerate(r):
                        #sph = self.sphere_coordinates(radius,curNode,10)
                        sph = self.cylinder_coordinates(radius,curNode,nodes[ni+1],self.feNSample)
                        #sph = self.circle_coordinates(radius,curNode,nodes[ni+1],self.feNSample)
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
    
                                if flatten_z:
                                    self.grid[:,xsc,ysc] += c_i[ri,:]
                                    self.count[:,xsc,ysc] += 1
                                else:
                                    self.grid[:,xsc,ysc,zsc] += c_i[ri,:]
                                    self.count[:,xsc,ysc,zsc] += 1
                                    
            if ni==nnode-1 and index[ni]==index[ni-1]:
                conc_out[ni,:] = conc_out[ni-1,:]

        if med_reg and store_results:
            pbar = tqdm(total=self.nt)
            print 'Median regularisation...'
            for i in xrange(self.nt):
                pbar.update(1)
                self.grid[i] = scipy.ndimage.filters.median_filter(self.grid[i],size=13)
            pbar.close()

        if progress:
            pbar.close()
        #self.normalise()
            
        return conc_out
        
    def normalise(self):
        self.grid[self.count>0] = self.grid[self.count>0] / self.count[self.count>0]
        
    def smooth_grid(self,n,grid=None):
        if grid is None:
            grid = self.grid
        for i in xrange(grid.shape[0]):
            grid[i] = scipy.ndimage.filters.median_filter(grid[i],size=n)
        return grid
        
    def save_grid(self,path,grid=None,pixdim=None,format='nifti'):
                
        if grid is None:
            grid = self.grid
            
        print('Max C_i: {} mM'.format(np.max(grid)))
            
        if format.lower()=='amira':
            pass
        else:
            import nibabel as nib
            ndims = len(grid.shape)
            if ndims==3:
                tmp = np.swapaxes(grid,0,1)
                tmp = np.swapaxes(tmp,1,2)
                img = nib.Nifti1Image(tmp[:,:,None,:],affine=np.eye(4))
            elif ndims==4:
                tmp = np.swapaxes(grid,0,1) # Shift time from first to last dimension
                tmp = np.swapaxes(tmp,1,2)
                tmp = np.swapaxes(tmp,2,3)
                img = nib.Nifti1Image(tmp,affine=np.eye(4))
            hdr = img.header
            if pixdim is None:
                hdr['pixdim'][1] = self.dx
                hdr['pixdim'][2] = self.dy
                hdr['pixdim'][3] = self.dz
                hdr['pixdim'][4] = self.dt
            else:
                hdr['pixdim'][1] = pixdim[0]
                hdr['pixdim'][2] = pixdim[1]
                hdr['pixdim'][3] = pixdim[2]
                hdr['pixdim'][4] = pixdim[3]
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
        
def main():
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
    
if __name__ == "__main__":
    pass
    #main()
