# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 09:37:29 2017

@author: simon

Statistical analysis of Amira SpatialGraph file

"""

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pickle
import nibabel as nib
from tqdm import tqdm, trange # progress bar
import os

class Statistics(object):
    
    def __init__(self, graph, path=None):
        
        self.graph = graph
        self.nodes = None
        self.edges = None
        rad_names = ['Radii','thickness','Radius']
        self.radius_field_name = None
        for rad_name in rad_names:
            radii = self.graph.get_data(rad_name)
            if radii is not None:
                self.radius_field_name = rad_name
                
        self.radii = None
        self.nconn = None
        self.branching_angles = None
        self.node_connections = None
        
        self.edge_intervessel_distance = None
        self.edge_length = None
        self.edge_volume = None
        self.edge_radii = None
        self.edge_euclidean = None      
        self.edge_tortuosity = None                   
        
        #if self.nodes is None:
        #    print('Generating node list...')
        #    self.nodes = self.graph.node_list(path=path)
        #    print('Node list complete')
        #    print('Generating node geometry...')
        #    self.node_geometry(self.nodes)
        #    print('Node geometry complete')
        #if False: #self.edges is None:
        #    print('Generating edge list...')
        #    #self.edges = self.graph.edges_from_node_list(self.nodes)
        #    print('Edge list complete')
            
        print('Generating node geometry...')
        self.node_geometry()
        print('Generating edge geometry...')
        self.edge_geometry()
        print('Edge geometry complete')
        
        
    def coords_in_cube(self,coords,cube):
        
        res = np.zeros(len(coords),dtype='bool')
        for i,c in enumerate(coords):
            if c[0]>=cube[0] and c[0]<=cube[1] and \
               c[1]>=cube[2] and c[1]<=cube[3] and \
               c[2]>=cube[4] and c[2]<=cube[5]:
                res[i] = True
            else:
                res[i] = False
        return res
        
    def summary_image(self,voxel_size=[250,250,250.],parameter='Radii',output_path=''):

        nse = self.graph.edge_spatial_extent()
        dnse = [np.abs(x[1]-x[0]) for x in nse]
        
        nstep = [np.int(np.ceil(dnse[i] / np.float(voxel_size[i]))) for i in range(3)]
        pixel_volume = np.product(voxel_size)
        
        volume = np.zeros(nstep)
        bfrac = np.zeros(nstep)
        radius = np.zeros(nstep)
        length = np.zeros(nstep)
        flow = np.zeros(nstep)
        count = np.zeros(nstep)
        
        x = np.linspace(nse[0][0],nse[0][1],num=nstep[0])
        y = np.linspace(nse[1][0],nse[1][1],num=nstep[1])
        z = np.linspace(nse[2][0],nse[2][1],num=nstep[2])
        
        pbar = tqdm(total=len(self.edges))
                    
        for ei,edge in enumerate(self.edges):
            pbar.update(1)

            radii,lengths,volumes,coords,flows = self.blood_volume([edge],sum_edge=False)
            
            for cInd,coords in enumerate(edge.coordinates):
                xInds = [i1 for i1,xp in enumerate(x) if (xp-coords[0])>=0]
                yInds = [i1 for i1,yp in enumerate(y) if (yp-coords[1])>=0]
                zInds = [i1 for i1,zp in enumerate(z) if (zp-coords[2])>=0]
                try:
                    i = xInds[0]
                    j = yInds[0]
                    k = zInds[0]
                except Exception as e:
                    print(e)
                    import pdb
                    pdb.set_trace()
                    
                #print i,j,k,len(volumes),coords.shape
                
                if cInd<volumes.shape[0]:
                    volume[i,j,k] += volumes[cInd]
                    if flows is not None:
                        flow[i,j,k] += flows[cInd]
                        
                    radius[i,j,k] += radii[cInd]
                    length[i,j,k] += lengths[cInd]
                    count[i,j,k] += 1
                    
        pbar.close()
        
        # Take averages
        radius = radius / count
        radius[~np.isfinite(radius)] = 0.
        length = length / count
        length[~np.isfinite(length)] = 0.
        bfrac = volume / pixel_volume

        img = nib.Nifti1Image(bfrac,affine=np.eye(4))
        ofile = output_path+'blood_volume.nii'
        nib.save(img,ofile)
    
        img = nib.Nifti1Image(radius,affine=np.eye(4))
        ofile = output_path+'radius.nii'
        #ofile = 'C:\\Users\\simon\\Dropbox\\radius.nii'
        nib.save(img,ofile)
        
        img = nib.Nifti1Image(length,affine=np.eye(4))
        #ofile = 'C:\\Users\\simon\\Dropbox\\length.nii'
        ofile = output_path+'length.nii'
        nib.save(img,ofile)
        
        img = nib.Nifti1Image(flow,affine=np.eye(4))
        #ofile = 'C:\\Users\\simon\\Dropbox\\flow.nii'
        ofile = output_path+'flow.nii'
        nib.save(img,ofile)
        
    def edge_geometry(self): #,edges):
        
        #pbar = tqdm(total=len(edges))

        #edgeInds = self.graph.edgepoint_edge_indices()
        #points = self.graph.get_data('EdgePointCoordinates')
        #radii = self.graph.get_data(self.radius_field_name)
        
        graph = self.graph
        nodecoords = graph.get_data('VertexCoordinates')
        edgeconn = graph.get_data('EdgeConnectivity')
        edgepoints = graph.get_data('EdgePointCoordinates')
        nedgepoints = graph.get_data('NumEdgePoints')
        radii = self.graph.get_data(self.radius_field_name)
        
        nedges = edgeconn.shape[0]
        
        edgeInds = np.zeros(edgepoints.shape[0],dtype='int')
        for edge_ind in range(nedges):
            nep = nedgepoints[edge_ind]
            x0 = np.sum(nedgepoints[:edge_ind])
            x1 = x0 + nep
            edgeInds[x0:x1] = edge_ind
 
        self.edge_intervessel_distance = np.zeros(nedges)   
        self.edge_length = np.zeros(nedges)
        self.edge_radii = np.zeros(nedges)
        self.edge_volume = np.zeros(nedges)
        self.edge_euclidean = np.zeros(nedges)      
        self.edge_tortuosity = np.zeros(nedges)   

        for edge_ind in trange(nedges):
            #try:
            if True:
                nep = nedgepoints[edge_ind]
                x0 = np.sum(nedgepoints[:edge_ind])
                x1 = x0 + nep
                pts = edgepoints[x0:x1]
                rads = radii[x0:x1]
            
                #pts = edge.coordinates
                #rads = edge.get_scalar(self.radius_field_name)
                
                # Define search range
                rng = [np.min(pts,axis=0),np.max(pts,axis=0)]
                #if rng[0]<0.:
                #   rng[0] *= 1.2*np.max(rads)
                #else:
                #   rng[0] *= 0.8*np.max(rads)
                #if rng[1]<0.:
                #   rng[1] *= 0.8*np.max(rads)
                #else:
                #   rng[1] *= 1.2*np.max(rads)
                   
                lim = 10. #100. #um
                # Get nearby points
                #inds = [(edgepoints[:,0]>rng[0][0]-lim) & (edgepoints[:,0]<rng[1][0]+lim) & (edgepoints[:,1]>rng[0][1]-lim) & (edgepoints[:,1]<rng[1][1]+lim) & (edgepoints[:,2]>rng[0][2]-lim) & (edgepoints[:,2]<rng[1][2]+lim) & (edgeInds!=edge_ind)]
                inds = []
                if len(inds)>0:
                    curPoints = edgepoints[inds[0]]
                    curRadii = radii[inds[0]]
                else:
                    curPoints,curRadii = [],[]

                dist = np.zeros(pts.shape[0]-1)
                
                length = np.sum([np.linalg.norm(pts[i]-pts[i-1]) for i,x in enumerate(pts[1:])])
                volume = np.sum([np.linalg.norm(pts[i]-pts[i-1])*np.square(rads[i]) for i,x in enumerate(pts[1:])])
                
                for i in range(pts.shape[0]-1):
                    if len(curPoints)>0:
                        dist[i] = np.min([np.linalg.norm(pts[i]-p)-rads[i]-curRadii[j] for j,p in enumerate(curPoints)])
                        
                if len(curPoints)>0:
                    self.edge_intervessel_distance[edge_ind] = np.max(dist)
                else:
                    self.edge_intervessel_distance[edge_ind] = -1.
                self.edge_length[edge_ind] = length
                self.edge_volume[edge_ind] = volume
                self.edge_radii[edge_ind] = np.mean(rads)
                self.edge_euclidean[edge_ind] = np.linalg.norm(pts[-1]-pts[0])

                self.edge_tortuosity[edge_ind] = self.edge_euclidean[edge_ind] / length
                
            #except Exception as e:
            #    print('Error, edge {}: {}'.format(edge,e))

            
    def _branching_angle(self,vec1,vec2):

        if np.linalg.norm(vec1)*np.linalg.norm(vec2)==0.:
            return 0.
        rad = np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
        deg = np.rad2deg(rad)
        return deg
        
    def volume(self,coords):
        from scipy.spatial import ConvexHull
        return ConvexHull(coords).volume
            
    def node_geometry(self): #,nodes):
    
        graph = self.graph
        nodecoords = graph.get_data('VertexCoordinates')
        edgeconn = graph.get_data('EdgeConnectivity')
        edgepoints = graph.get_data('EdgePointCoordinates')
        nedgepoints = graph.get_data('NumEdgePoints')
        radii = self.graph.get_data(self.radius_field_name)
        
        nnodes = nodecoords.shape[0]
        
        self.branching_angles = []
        self.node_connections = []

        for node_ind in trange(nnodes):
            sind = np.where((edgeconn[:,0]==node_ind) | (edgeconn[:,1]==node_ind))

            if len(sind[0])>0:
                if len(sind[0])>1:
                    self.node_connections.append(len(sind[0]))
                    
                for edge_ind in sind[0]:
                
                    # Edge direction
                    if edgeconn[edge_ind,0]==node_ind:
                        direction = 1
                    else:
                        direction = -1
                        
                    nep = nedgepoints[edge_ind]
                    x0 = np.sum(nedgepoints[:edge_ind])
                    x1 = x0 + nep
                    pts = edgepoints[x0:x1]
                    
                    if direction==-1:
                        pts = pts[::-1]
                        
                    for edge_ind2 in sind[0]:
                        if edge_ind!=edge_ind2:
                            
                            # Edge direction
                            if edgeconn[edge_ind2,0]==node_ind:
                                direction2 = 1
                            else:
                                direction2 = -1
                                
                            nep2 = nedgepoints[edge_ind2]
                            x02 = np.sum(nedgepoints[:edge_ind2])
                            x12 = x02 + nep2
                            pts2 = edgepoints[x02:x12]
                            
                            if direction2==-1:
                                pts2 = pts2[::-1]

                            veci = pts[0]-pts[1]
                            vecj = pts2[0]-pts2[1]
                            if not all(x==y for x,y in zip(veci,vecj)):
                                self.branching_angles.append(self._branching_angle(veci,vecj))
        
    def histogram(self,v,range=None,xlabel=None,nbins=50,show=False):
        
        # the histogram of the data
        n, bins, patches = plt.hist(v, nbins, range=range, normed=1, facecolor='green', alpha=0.75)
        
        # add a 'best fit' line
        #y = mlab.normpdf( bins, mu, sigma)
        #l = plt.plot(bins, y, 'r--', linewidth=1)
        
        if xlabel is not None:
            plt.xlabel(xlabel)
        #plt.ylabel('Probability')
        #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
        #plt.axis([40, 160, 0, 0.03])
        #plt.grid(True)
        
        if show:
            plt.show()
        
    def boxplot(self,v):
        
        plt.boxplot(v)
        
    def blood_volume(self,edges,sum_edge=True):
        
        nedge = len(edges)
        
        radii = [None] * nedge #np.zeros(self.graph.nedgepoint)
        lengths = [None] * nedge #np.zeros(self.graph.nedgepoint)
        volumes = [None] * nedge #np.zeros(self.graph.nedgepoint)
        coords = [None] * nedge #np.zeros((self.graph.nedgepoint,3))
        if 'Flow' in edges[0].scalarNames:
            flow = [None] * nedge #np.zeros(self.graph.nedgepoint)
        else:
            flow = None
             
        pbar = tqdm(total=len(edges))
        #import pdb
        #pdb.set_trace()
        #pcount = 0
        for ei,edge in enumerate(edges):
            pbar.update(1)
            npoints = edge.npoints
            rad = edge.get_scalar(self.radius_field_name)
            radii[ei] = rad
            lengths[ei] = np.append(edge.length,0)
            curVol = np.pi*np.square(rad[0:-1])*edge.length
            volumes[ei] = curVol
            if flow is not None:
                curFlow = edge.get_scalar('Flow')
                flow[ei] = curFlow
            
            if sum_edge:
                volumes[ei] = np.sum(volumes[ei])
                lengths[ei] = np.sum(lengths[ei])
            coords[ei] = edge.coordinates
            
            #pcount += npoints
        pbar.close()
            
        return radii,lengths,volumes,coords,flow
        
    def do_stats(self,path=None):
        
        #print('Calculating statistics...')
        #print('Estimating network parameters...')
        #import pdb
        #pdb.set_trace()
        #radii,lengths,volumes,coords,flow = self.blood_volume(self.edges)
        #print('Finished estimating network parameters...')

        print('Calculating stats...')
        nconn = np.asarray(self.node_connections)
        ba = np.asarray(self.branching_angles)
        ba = ba[np.isfinite(ba)]
        ivd = self.edge_intervessel_distance[self.edge_intervessel_distance>0.]

        self.histogram(ba,range=[0,180],nbins=50,xlabel='Vessel branching angle (deg)')
        if path is not None:
            plotfile = os.path.join(path,'branching_angle_histogram.png')
            plt.savefig(plotfile) #,transparent=True)

        self.histogram(self.edge_radii,nbins=30,xlabel='Vessel radius (um)')
        if path is not None:
            plotfile = os.path.join(path,'radius_histogram.png')
            plt.savefig(plotfile) #,transparent=True)

        #lengthsFlat = [item for sublist in lengths for item in sublist]
        self.histogram(self.edge_length,nbins=50,xlabel='Vessel length (um)')
        if path is not None:
            plotfile = os.path.join(path,'vessel_length_histogram.png')
            plt.savefig(plotfile) #,transparent=True)

        self.histogram(self.edge_volume,nbins=50,xlabel='Vessel volume (um3)')
        if path is not None:
            plotfile = os.path.join(path,'vessel_volume_histogram.png')
            plt.savefig(plotfile) #,transparent=True)
            
        self.histogram(self.edge_tortuosity,nbins=50,xlabel='Vessel volume (um3)')
        if path is not None:
            plotfile = os.path.join(path,'vessel_volume_histogram.png')
            plt.savefig(plotfile) #,transparent=True)
            
        self.histogram(ivd,nbins=50,xlabel='Intervessel distance (um)')
        if path is not None:
            plotfile = os.path.join(path,'intervessel_distance_histogram.png')
            plt.savefig(plotfile) #,transparent=True)

        #blood_volume = np.sum(volumes)
        #coords = self.graph.get_data('VertexCoordinates')
        #try:
        #    total_volume = self.volume(coords)
        #except Exception as e:
        #    print(e)
        #    total_volume = -1.
        #blood_fraction = blood_volume / total_volume
        #print(('TUMOUR VOLUME (um3): {}'.format(total_volume)))
        #print(('BLOOD VOLUME (um3): {}'.format(blood_volume)))
        #print(('BLOOD VOLUME FRACTION: {}'.format(blood_fraction)))
        
        if path is not None:
            #with open(os.path.join(path,'volume.txt'),'w') as fo:
            #    fo.write('TUMOUR VOLUME (um3) \t{}\n'.format(total_volume))
            #    fo.write('BLOOD VOLUME (um3) \t{}\n'.format(blood_volume))
            #    fo.write('BLOOD VOLUME FRACTION \t{}\n'.format(blood_fraction))
                
            with open(os.path.join(path,'stats.txt'),'w') as fo:
                
                # Write header
                hdr = ['PARAM','Mean','SD','median','min','max']
                hdr = '\t'.join(hdr)+'\n'
                fo.write(hdr)
                
                params = [self.edge_radii,self.edge_length,ivd,nconn,ba,self.edge_volume]
                paramNames = ['Radius (um)','Vessel length (um)','Intervessel distance( um)','Number of connections','Branching angle (deg)','Vessel volume (um3)']
                
                for v,n in zip(params,paramNames):
                    print(n)
                    try:
                        cur = [n,np.mean(v),np.std(v),np.median(v),np.min(v),np.max(v)]
                    except:
                        cur = [n,-1.,-1.,-1.,-1.,-1.]
                    cur = ['{}'.format(c) for c in cur]
                    cur = '\t'.join(cur)+'\n'
                    fo.write(cur)
                    
            with open(os.path.join(path,'radii.p'),'wb') as fo:
                pickle.dump(self.edge_radii,fo)
            with open(os.path.join(path,'branching_angle.p'),'wb') as fo:
                pickle.dump(ba,fo)
            with open(os.path.join(path,'vessel_length.p'),'wb') as fo:
                pickle.dump(self.edge_length,fo)
            with open(os.path.join(path,'nconn.p'),'wb') as fo:
                pickle.dump(nconn,fo)
            with open(os.path.join(path,'vessel_volume.p'),'wb') as fo:
                pickle.dump(self.edge_volume,fo)
            with open(os.path.join(path,'intervessel_distance.p'),'wb') as fo:
                pickle.dump(ivd,fo)

def main():
    #dir_ = 'C:\\Users\\simon\\Dropbox\\Mesentery\\'
    #f = dir_ + 'Flow2AmiraPressure.am'
    # dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\'
    # #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
    
    # #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\'
    
    # #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\'
    # #f = dir_+'spatialGraph_RIN.am'
    
    #dir_ = r"C:\Users\simon\Dropbox\VDA_1_lectin\Control\SW#1"
    #f = dir_+r'\SW1_spatialGraph_scaled.am'
    # pixsize = 6.98
    # dir_ = r"G:\OPT\2015.11.VDA_1 study\VDA Colorectal cancer\Control\LS\LS#2"
    # f = dir_+r'\LS2_bg_removed_frangi_response_skeletonised_with_radius.SptGraph.am'
    # pixsize = 8.21
    # #dir_ = r"G:\OPT\2015.11.VDA_1 study\VDA Colorectal cancer\Control\LS\LS#4"
    # #pixsize = 8.21
    
    #dirs = [r"C:\Users\simon\Dropbox\VDA_1_lectin\Control\SW#2",
    #        r"C:\Users\simon\Dropbox\VDA_1_lectin\Control\SW#3"]
    #fs = [r'\SW2_spatialGraph_scaled.am',
    #      r'\SW3_spatialGraph_scaled.am']
    #dirs = [r'C:\Users\simon\Dropbox\VDA_1_lectin\Treated\LS#1']
    #fs = [r'\LS1t_vessel_seg_frangi_response_skel_with_radius.am']
    #pixsize = [4.78]
    
    dirs = '/mnt/data2/Sahai/export/nifti'
    #fs = 'HET7_proc.SptGraph.am'
    fs = 'HET7_proc.SptGraph.am'
    odir = os.path.join(dirs,'HET7')
    if not os.path.exists(odir):
        os.mkdir(odir)
    pixsize = 1.
          
    #for i,dir_ in enumerate(dirs):
    i = 0
    if True:
        
        f = fs #[i]
        pix = pixsize #[i]
        dir_ = dirs
          
        from pymira import spatialgraph
        graph = spatialgraph.SpatialGraph()
        print('Reading graph...')
        graph.read(os.path.join(dir_,f))
        print('Graph read')
       
        if pix is not None and pix!=1.:
            ofile = os.path.join(dir_,'spatialGraph_scaled.am')
            graph.rescale_coordinates(pix,pix,pix)
            graph.rescale_radius(pix,ofile=ofile)
        
        stats = Statistics(graph,path=dir_)  
        stats.do_stats(path=odir)
    
    # ofile = dir_+'\spatialGraph_scaled.am'
    # graph.rescale_coordinates(pixsize,pixsize,pixsize)
    # graph.rescale_radius(pixsize,ofile=ofile)
    
    # #epi = graph.edge_point_index()
    # #edgeCoords = graph.get_data('EdgePointCoordinates')
    # #nodes = graph.node_list()
    # #edges = graph.edges_from_node_list(nodes)
    # #testInd = 1000
    # #testCoords = edgeCoords[testInd]
    # #testEdge = [e for e in edges if e.index==epi[testInd]]
    # #testCoords in testEdge[0].coordinates
    # #import pdb
    # #pdb.set_trace()
    
    # #import pdb
    # #pdb.set_trace()
    #stats = Statistics(graph,path=dir_)
    
    #stats.do_stats(path=dir_)
    # #stats.do_stats(output_directory=None)
    # #stats.summary_image(voxel_size=[125.,125.,125.],output_path=dir_)
    # stats.do_stats(output_directory=dir_)
    # #stats.do_stats(output_directory=None)
    # #stats.summary_image(voxel_size=[250.,250.,250.])
    # #stats.do_stats(output_directory=dir_)
    # stats.do_stats(output_directory=None)
    # stats.summary_image(voxel_size=[125.,125.,125.],output_path=dir_)
    
    # #stats.do_stats(output_directory=dir_)
    # stats.do_stats(output_directory=None)
    # #stats.summary_image(voxel_size=[250.,250.,250.])
    # #stats.do_stats(output_directory=dir_)
    # #stats.do_stats(output_directory=None)
    # stats.summary_image(voxel_size=[125.,125.,125.],output_path=dir_)

if __name__=='__main__':
#    try:
#        #import cProfile
#        #cProfile.run('main()')
#    except:
#        pass
#    #import pdb
#    #pdb.set_trace()
    main()
