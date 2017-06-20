# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:30:01 2017

@author: simon

Estimate exposure to a drug (uses injection simulaiton results from injectagent.py)

"""

import numpy as np
from pymira import spatialgraph, inject_agent
import pickle
import os

def load_timeseries_graph(dir_):
    
    from pymira import spatialgraph
    
    graph = []
    time = []
    
    for file in os.listdir(dir_):
        if file.endswith('.am'):
            #print(os.path.join("/mydir", file))
            g = spatialgraph.SpatialGraph()
            try:
                g.read(os.path.join(dir_,file))
                graph.append(g)
                t = np.float(g.get_parameter_value('Time'))
                time.append(t)
            except Exception,e:
                print e
                
    time = np.asarray(time)
    graph = np.asarray(graph)
    srtInds = time.argsort()
    time = time[srtInds]
    graph = graph[srtInds]
    
    return graph,time
    
def concentration_from_timeseries(graphs,log=True):
    
    nt = len(graphs)
    conc = None
    
    for gi,graph in enumerate(graphs):
        concFieldInd = [i for i,x in enumerate(graph.fieldNames) if 'Concentration' in x]
        
        field = graph.fields[concFieldInd[0]]
        
        if 'data' in field:
            if conc is None:
                npoint = len(field['data'])
                conc = np.zeros((nt,npoint),dtype='float')
            if not log:
                conc[gi,:] = field['data']
            else:
                conc[gi,:] = np.exp(field['data'])
        else:
            conc[gi,:] = conc[gi-1,:]
            print('Data missing: {}'.format(gi))
            
    return conc

def main():
    
    dir_ = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\ca1_kt0p0001'
    
    if True:
        
        editor = spatialgraph.Editor()
    
        print('Reading graph...')
        graphs,time = load_timeseries_graph(os.path.join(dir_,'vascular_recon'))
        print('Graph read')
        
        #nt = len(graphs)
        points = graphs[0].get_data('EdgePointCoordinates')
        npoints = points.shape[0]
        nEdgePoint = graphs[0].get_data('NumEdgePoints')
        edgePointIndex = np.zeros(npoints,dtype='int')
        
        offset = 0
        edgeCount = 0
        for npCur in nEdgePoint:
            edgePointIndex[offset:offset+npCur] = edgeCount
            edgeCount += 1
            offset += npCur
            
        assert offset==npoints

        #npoint = points.shape[0]
        
        conc = concentration_from_timeseries(graphs)
        
        radii = graphs[0].get_data('Radii')
        epi = graphs[0].edgepoint_edge_indices()
        
        #dt = np.ediff1d(time)
        #dt = np.append(0.,dt)
        #import pdb
        #pdb.set_trace()
        exposure = np.cumsum(conc,axis=0) / (2.*np.pi*radii)
        removedEdgePoint = np.zeros(conc.shape,dtype='int')
        #thr = np.linspace(np.max(exposure),np.min(exposure),num=10)
        #thr = np.max(exposure) / 100.
        #thr = np.median(exposure[-1])
        thr = 1e-14
        import pdb
        pdb.set_trace()
        removedEdgePoint[exposure>thr] = 1
        
        nodeFile = os.path.join(dir_,'nodeList.dill')
        if os.path.isfile(nodefile):
            with open(nodeFile ,'rb') as fo:
                nodeList = pickle.load(fo)
        else:
            import pdb
            pdb.set_trace()
        edges = graphs[0].edges_from_node_list(nodeList)
        
        #import pdb
        #pdb.set_trace()
        inj = inject_agent.InjectAgent()
        
        #def add_concentration(self,edges,time,conc_time=0.):
        #out_time = np.asarray([time[0],time[int(len(time)/2.)],time[-1]])
        out_time = np.asarray([time[-1]])
        #out_time = [time[-1]] #np.asarray(time[-1])
        #out_time = [time[-1]]
        
        #import dill
        #filename = dir_+'globalsave.pkl'
        #dill.dump_session(filename)
    else:
        import dill
        filename = dir_+'globalsave.pkl'
        dill.load_session(filename)
        import pdb
        pdb.set_trace()
        
    #import pdb
    #pdb.set_trace()
      
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
                        
        new_graph = graphs[0].node_list_to_graph(nodeList)
        remEdges = True
        if len(remEdgeInds)>0 and remEdges is True:
            #new_graph = editor.delete_edges(new_graph,remEdgeInds)
            nl = new_graph.node_list()
            inletNodes = []
            for node in nl:
                inj.vertex_flow_ordering(node)                    
                if node.nconn>0 and node.is_inlet:
                    inletNodes.append(node.index)
                    
            #import pdb
            #pdb.set_trace()
            
            new_graph.set_graph_sizes() # shouldn't have to do this...
            graphList, graphSize = new_graph.identify_graphs()
            graphInds = np.unique(graphList)
            # See if any graphs have been isolated from an inlet
            
            #remEdgeInds = []
                            
            largestNet = True
            if not largestNet:
                for graphInd in graphInds:
                    graphNodes = [nl[i] for i,g in enumerate(graphList) if g==graphInd]

                    hasInlet = np.any([n.nconn==1 and n.flow_direction[0]>0 for n in graphNodes])
                    hasOutlet = np.any([n.nconn==1 and n.flow_direction[0]<0 for n in graphNodes])
                    if not hasInlet or not hasOutlet:
                        remEdgeIndsCur = [[e.index for e in n.edges] for n in graphNodes]
                        remEdgeIndsCur = [[x] if type(x) is not list else x for x in remEdgeIndsCur]
                        remEdgeIndsCur = [item for sublist in remEdgeIndsCur for item in sublist] # flatten
                        remEdgeInds.append(remEdgeIndsCur)
            else:
                #import pdb
                #pdb.set_trace()
                if len(graphSize)>0:
                    #graph_size = np.histogram(graphList,bins=np.linspace(0,np.max(graphList),1))[0]
                    mxGraphInd = np.argmax(graphSize)
                    graphNodes = [nl[i] for i,g in enumerate(graphList) if g!=mxGraphInd]
                    remEdgeInds.append([[e.index for e in n.edges] for n in graphNodes])

            remEdgeInds = [[x] if type(x) is not list else x for x in remEdgeInds] # make sure all elemenets are lists
            remEdgeInds = [item for sublist in remEdgeInds for item in sublist] # flatten
            count += len(remEdgeInds)
            remEdgeInds = np.unique(remEdgeInds)
            new_graph = editor.delete_edges(new_graph,remEdgeInds)

        print 't={}, thr={}: {} {}%'.format(t,thr,count,count*100./float(len(edges)))
        new_graph.write(os.path.join(dir_,'exposure\exposure{}.am'.format(t)))
            
    #new_graph = graph.node_list_to_graph(nodeList)
    #new_graph.write(dir_+'ct_output_with_exposure.am')
        
if __name__ == "__main__":
    main()