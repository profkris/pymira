# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:30:01 2017

@author: simon
"""

import numpy as np
from pymira import spatialgraph, inject_agent
import pickle

def main():         
    #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T - Post-VDA\\1\\'
    #dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\'
    #f = os.path.join(dir_,'spatialGraph_RIN.am')
    dir_ = 'C:\\Users\\simon\\Dropbox\\Mesentery\\'
    f = dir_ + 'ct_output.am'
    
    editor = spatialgraph.Editor()

    graph = spatialgraph.SpatialGraph()
    print('Reading graph...')
    graph.read(f)
    print('Graph read')
    
#    nodeFile = dir_+'nodeList.dill'
#    with open(nodeFile ,'rb') as fo:
#        nodeList = pickle.load(fo)
#    edges = graph.edges_from_node_list(nodeList)
    
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
    
    radii = graph.get_data('Radii')
    epi = graph.edgepoint_edge_indices()
    
    dt = np.ediff1d(time)
    dt = np.append(0.,dt)
    exposure = np.cumsum(conc,axis=0) / (2.*np.pi*radii)
    removedEdgePoint = np.zeros(conc.shape,dtype='int')
    #thr = np.linspace(np.max(exposure),np.min(exposure),num=10)
    #thr = np.max(exposure) / 100.
    #thr = np.median(exposure[-1])
    thr = 1e-12
    removedEdgePoint[exposure>thr] = 1
    
    nodeFile = dir_+'nodeList.dill'
    with open(nodeFile ,'rb') as fo:
        nodeList = pickle.load(fo)
    edges = graph.edges_from_node_list(nodeList)
    
    #import pdb
    #pdb.set_trace()
    inj = inject_agent.InjectAgent()
    inletNodes = []
    for node in nodeList:
        inj.vertex_flow_ordering(node)
        if node.is_inlet:
            inletNodes.append(node.index)
    
    #def add_concentration(self,edges,time,conc_time=0.):
    #out_time = np.asarray([time[0],time[int(len(time)/2.)],time[-1]])
    out_time = [time[-1]] #np.asarray(time[-1])
    for idx,t in enumerate(out_time):
        #idx = (np.abs(time-conc_time)).argmin()
        #conc = np.zeros([len(edges),edges[0].concentration.shape[0]])
        
        count = 0
        remEdgeInds = []
        for i,edge in enumerate(edges):
            curVals = exposure[idx,epi==edge.index]
            if idx==0:
                edge.add_scalar('Exposure',curVals)
            else:
                edge.set_scalar('Exposure',curVals)
            sind = np.where(curVals>thr)
            if len(sind[0])>0:
                count += 1
                #curVals[:] = 0
                remEdgeInds.append(edge.index)
                #new_graph.remove_edge(edge.index)
            else:
                #curVals[:] = 1
                pass
            #edge.add_scalar('ExposureThreshold_{}'.format(t),curVals)
            #edge.add_scalar('ExposureThreshold'.format(t),curVals)
                        
        new_graph = graph.node_list_to_graph(nodeList)
        remEdges = True
        if len(remEdgeInds)>0 and remEdges is True:
            new_graph = editor.delete_edges(new_graph,remEdgeInds)
            nl = new_graph.node_list()
            graphList = new_graph.identify_graphs()
            graphInds = np.unique(graphList)
            # See if any graphs have been isolated from an inlet
            for graphInd in graphInds:
                graphNodes = [nodeList[i] for i,g in enumerate(graphList) if g==graphInd]
                hasInlet = np.any([n.nconn==1 and n.flow_direction[0]>0 for n in graphNodes])
                hasOutlet = np.any([n.nconn==1 and n.flow_direction[0]<0 for n in graphNodes])
                if not hasInlet or not hasOutlet:
                    remEdgeInds = [[e.index for e in n.edges] for n in graphNodes]
                    remEdgeInds = [item for sublist in remEdgeInds for item in sublist]
                    new_graph = editor.delete_edges(new_graph,remEdgeInds)

        print 't={}, thr={}: {} {}%'.format(t,thr,count,count*100./float(len(edges)))
        new_graph.write(dir_+'exposure\exposure{}.am'.format(t))
            
    #new_graph = graph.node_list_to_graph(nodeList)
    #new_graph.write(dir_+'ct_output_with_exposure.am')
        
if __name__ == "__main__":
    main()