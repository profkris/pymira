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
import progressbar

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
        
        #bar = progressbar.ProgressBar(maxval=len(self.edges),redirect_stdout=True)#,max_value=len(self.edges))
                    
        for ei,edge in enumerate(self.edges):
            #bar.update(ei)

            radii,lengths,volumes,coords,flows = self.blood_volume([edge],sum_edge=False)
            
            for cInd,coords in enumerate(edge.coordinates):
                xInds = [i1 for i1,xp in enumerate(x) if (xp-coords[0])>=0]
                yInds = [i1 for i1,yp in enumerate(y) if (yp-coords[1])>=0]
                zInds = [i1 for i1,zp in enumerate(z) if (zp-coords[2])>=0]
                try:
                    i = xInds[0]
                    j = yInds[0]
                    k = zInds[0]
                except Exception,e:
                    print e
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
        
        img = nib.Nifti1Image(flows,affine=np.eye(4))
        #ofile = 'C:\\Users\\simon\\Dropbox\\flow.nii'
        ofile = output_path+'flow.nii'
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

#dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\SW1222\\1\\'
f = dir_+'spatialGraph_RIN.am'

#dir_ = r"G:\OPT\2015.11.VDA_1 study\VDA Colorectal cancer\Control\LS\LS#1"
#f = dir_+r'\ls1_vessel_seg_skel_with_radius.SptGraph.am'

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

#import pdb
#pdb.set_trace()
stats = Statistics(graph)
stats.do_stats(output_directory=dir_)
#stats.do_stats(output_directory=None)
#stats.summary_image(voxel_size=[250.,250.,250.])
#stats.do_stats(output_directory=dir_)
stats.do_stats(output_directory=None)
stats.summary_image(voxel_size=[125.,125.,125.],output_path=dir_)
