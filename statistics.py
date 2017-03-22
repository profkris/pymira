# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 09:37:29 2017

@author: simon
"""

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pickle
import nibabel as nib

class Statistics(object):
    
    def __init__(self, graph):
        
        self.graph = graph
        self.nodes = None
        self.edges = None
        
        if self.nodes is None:
            print('Generating node list...')
            self.nodes = self.graph.node_list()
            print('Node list complete')
            print('Generating node geometry...')
            self.node_geometry(self.nodes)
            print('Node geometry complete')
        if self.edges is None:
            print('Generating edge list...')
            self.edges = self.graph.edges_from_node_list(self.nodes)
            print('Edge list complete')
            print('Generating edge geometry...')
            self.edge_geometry(self.edges)
            print('Edge geometry complete')
        
        self.radii = None
        self.nconn = None
        self.branching_angle = None
        self.lengths = None
        
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
        
    def summary_image(self,voxel_size=[250,250,250.],parameter='Radii'):
        
        #if parameter=='Radii':
        #    values = self.radii
        
        nodeCoords = self.graph.get_data('VertexCoordinates')
        edgeCoords = self.graph.get_data('EdgePointCoordinates')
         
        nse = self.graph.edge_spatial_extent()
        dnse = [np.abs(x[1]-x[0]) for x in nse]
        
        epi = self.graph.edge_point_index() # maps edge points to edge indicess
        
        nstep = [np.int(np.ceil(dnse[i] / np.float(voxel_size[i]))) for i in range(3)]
        
        blood_volume = np.zeros(nstep)
        bfrac = np.zeros(nstep)
        radius = np.zeros(nstep)
        length = np.zeros(nstep)
        flows = np.zeros(nstep)
        
        for i in range(nstep[0]):
            for j in range(nstep[1]):
                for k in range(nstep[2]):
                    x0 = i*voxel_size[0] + nse[0][0]
                    x1 = x0 + voxel_size[0] - 1
                    y0 = j*voxel_size[1] + nse[1][0]
                    y1 = y0 + voxel_size[1] - 1
                    z0 = k*voxel_size[2] + nse[2][0]
                    z1 = z0 + voxel_size[2] - 1

                    cin = self.coords_in_cube(edgeCoords,[x0,x1,y0,y1,z0,z1])

                    if np.any(cin):
                        edgePointInds = np.asarray([ind for ind,v in enumerate(cin) if v])
                        curEdges = np.unique([e for e in self.edges if e.index in epi[edgePointInds]])
                        
                        # Filter any points outside voxel
                        bvol = 0.
                        for ed in curEdges:
                            try:
                                cin2 = self.coords_in_cube(ed.coordinates,[x0,x1,y0,y1,z0,z1])
                                radii,lengths,volumes,coords,flow = self.blood_volume([ed],sum_edge=False)
                                volumes = np.append([0.],[volumes])
                                bvol += np.sum(volumes[cin2])
                            except Exception,e:
                                print(e)
                                import pdb
                                pdb.set_trace()
                        blood_volume[i,j,k] = bvol
                        if flow is not None:
                            flows[i,j,k] = np.mean(flow)
                        pixel_volume = np.product(voxel_size)
                        bfrac[i,j,k] = bvol / pixel_volume
                        radius[i,j,k] = np.mean(radii)
                        length[i,j,k] = np.mean(lengths)
                        print i,j,k,bfrac[i,j,k],radius[i,j,k],length[i,j,k]

        img = nib.Nifti1Image(bfrac,affine=np.eye(4))
        ofile = 'C:\\Users\\simon\\Dropbox\\blood_volume.nii'
        nib.save(img,ofile)
    
        img = nib.Nifti1Image(radius,affine=np.eye(4))
        ofile = 'C:\\Users\\simon\\Dropbox\\radius.nii'
        nib.save(img,ofile)
        
        img = nib.Nifti1Image(length,affine=np.eye(4))
        ofile = 'C:\\Users\\simon\\Dropbox\\length.nii'
        nib.save(img,ofile)
        
        img = nib.Nifti1Image(flow,affine=np.eye(4))
        ofile = 'C:\\Users\\simon\\Dropbox\\flow.nii'
        nib.save(img,ofile)
        
    def edge_geometry(self,edges):
        
        for edge in edges:        
            length = np.zeros(edge.npoints-1)
            pts = edge.coordinates
            for i in range(0,edge.npoints-1):
                length[i] = np.linalg.norm(pts[i+1]-pts[i]) # Length (um)
            edge.length = length
            edge.euclidean = np.linalg.norm(pts[-1]-pts[0])
            
    def _branching_angle(self,vec1,vec2):

        if np.linalg.norm(vec1)*np.linalg.norm(vec2)==0.:
            return 0.
        rad = np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
        deg = np.rad2deg(rad)
        return deg
        
    def volume(self,coords):
        from scipy.spatial import ConvexHull
        return ConvexHull(coords).volume
            
    def node_geometry(self,nodes):

        for node in nodes:
            ba = []
            for i in range(0,len(node.edges)):
                for j in range(0,len(node.edges)):
                    if i!=j:
                        ci = node.edges[i].coordinates
                        if node.edges[i].at_start_node(node.index) is False:
                            ci = ci[::-1]
                        veci = ci[0]-ci[1]
                        cj = node.edges[j].coordinates
                        if node.edges[j].at_start_node(node.index) is False:
                            cj = cj[::-1]
                        vecj = cj[0]-cj[1]
                        if not all(x==y for x,y in zip(veci,vecj)):
                            ba.append(self._branching_angle(veci,vecj))
            node.branching_angle = np.unique(ba)
        
    def histogram(self,v,range=None,xlabel=None,nbins=50):
        
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
        
        plt.show()
        
    def boxplot(self,v):
        
        plt.boxplot(v)
        
    def blood_volume(self,edges,sum_edge=True):
        
        radii = np.zeros(0)
        lengths = np.zeros(0)
        volumes = np.zeros(0)
        coords = []
        if 'Flow' in edges[0].scalarNames:
            flow = np.zeros(0)
        else:
            flow = None
            
        for edge in edges:
            rad = edge.get_scalar('Radii')
            radii = np.append(radii,rad)
            lengths = np.append(lengths,edge.length)
            curVol = np.pi*np.square(rad[1:])*edge.length
            volumes = np.append(volumes,curVol)
            if flow is not None:
                curFlow = edge.get_scalar('Flow')
                flow = np.append(flow,curFlow)
            
            if sum_edge:
                curVol = np.sum(curVol)
                lengths = np.sum(lengths)
            coords.append(edge.coordinates)
            
        return radii,lengths,volumes,coords,flow
        
    def do_stats(self,output_directory=None):
        
        print('Calculating statistics...')
        radii,lengths,volumes,_,_ = self.blood_volume(self.edges)

        nconn = np.asarray([node.nconn for node in self.nodes])
        ba = np.zeros(0)
        for node in self.nodes:
            ba = np.append(ba,node.branching_angle)
        ba = ba[np.isfinite(ba)]
            
        self.radii = radii
        self.nconn = nconn
        self.branching_angle = ba
        self.lengths = lengths
        self.volumes = volumes

        plt.figure()
        self.histogram(ba,range=[0,180],nbins=50,xlabel='Vessel branching angle (deg)')
        if output_directory is not None:
            plotfile = output_directory+'branching_angle_histogram.pdf'
            plt.savefig(plotfile,transparent=True)
        plt.figure()
        self.histogram(radii,range=[0,30],nbins=30,xlabel='Vessel radius (um)')
        if output_directory is not None:
            plotfile = output_directory+'radius_histogram.pdf'
            plt.savefig(plotfile,transparent=True)
        plt.figure()
        self.histogram(lengths,range=[0,80],nbins=50,xlabel='Vessel length (um)')
        if output_directory is not None:
            plotfile = output_directory+'vessel_length_histogram.pdf'
            plt.savefig(plotfile,transparent=True)
        self.histogram(volumes,range=[0,80],nbins=50,xlabel='Vessel volume (um3)')
        if output_directory is not None:
            plotfile = output_directory+'vessel_volume_histogram.pdf'
            plt.savefig(plotfile,transparent=True)

        blood_volume = np.sum(volumes)
        coords = self.graph.get_data('VertexCoordinates')
        total_volume = self.volume(coords)
        blood_fraction = blood_volume / total_volume
        print('TUMOUR VOLUME (um3): {}'.format(total_volume))
        print('BLOOD VOLUME (um3): {}'.format(blood_volume))
        print('BLOOD VOLUME FRACTION: {}'.format(blood_fraction))
        
        if output_directory is not None:
            with open(output_directory+'volume.txt','wb') as fo:
                fo.write('TUMOUR VOLUME (um3) \t{}\n'.format(total_volume))
                fo.write('BLOOD VOLUME (um3) \t{}\n'.format(blood_volume))
                fo.write('BLOOD VOLUME FRACTION \t{}\n'.format(blood_fraction))
                
            with open(output_directory+'stats.txt','wb') as fo:
                
                # Write header
                hdr = ['PARAM','Mean','SD','median','min','max']
                hdr = '\t'.join(hdr)+'\n'
                fo.write(hdr)
                
                params = [radii,lengths,nconn,ba,volumes]
                paramNames = ['Radius (um)','Vessel length (um)','Number of connections','Branching angle (deg)','Vessel volume (um3)']
                
                for v,n in zip(params,paramNames):
                    cur = [n,np.mean(v),np.std(v),np.median(v),np.min(v),np.max(v)]
                    cur = ['{}'.format(c) for c in cur]
                    cur = '\t'.join(cur)+'\n'
                    fo.write(cur)
                    
            with open(output_directory+'radii.p','wb') as fo:
                pickle.dump(radii,fo)
            with open(output_directory+'branching_angle.p','wb') as fo:
                pickle.dump(ba,fo)
            with open(output_directory+'vessel_length.p','wb') as fo:
                pickle.dump(lengths,fo)
            with open(output_directory+'nconn.p','wb') as fo:
                pickle.dump(nconn,fo)
            with open(output_directory+'vessel_volume.p','wb') as fo:
                pickle.dump(volumes,fo)
        
#dir_ = 'C:\\Users\\simon\\Dropbox\\Mesentery\\'
#f = dir_ + 'Flow2AmiraPressure.am'
#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\'
#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\'

f = dir_+'spatialGraph_RIN.am'
from pymira import spatialgraph
graph = spatialgraph.SpatialGraph()
print('Reading graph...')
graph.read(f)
print('Graph read')

#epi = graph.edge_point_index()
#edgeCoords = graph.get_data('EdgePointCoordinates')
#nodes = graph.node_list()
#edges = graph.edges_from_node_list(nodes)
#testInd = 1000
#testCoords = edgeCoords[testInd]
#testEdge = [e for e in edges if e.index==epi[testInd]]
#testCoords in testEdge[0].coordinates
#import pdb
#pdb.set_trace()

stats = Statistics(graph)
stats.do_stats(output_directory=dir_)
#stats.do_stats(output_directory=None)
#stats.summary_image(voxel_size=[250.,250.,250.])