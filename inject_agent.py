# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 09:23:48 2017

@author: simon
"""

from pymira import spatialgraph
import pymira.front as frontPKG
import numpy as np
import os
import pickle
import warnings

np.seterr(all='warn')
warnings.filterwarnings('error')

class InjectAgent(object):
    
    def __init__(self):
        self.dt = 0.1
        self.nt = 100
        self.time = np.arange(0,self.nt,self.dt)
        # To convert from (nL/min) to (um^3/s) use conversion below
        self.fluxConversion = 1e6/60.
        
        self.Q_limit = 1e-7
        
        self.nodeList = None
        self.graph = None
        
    def vertex_flow_ordering(self,node):
        
        if node.nconn==0:
            return
    
        node.inFlow = 0.
        node.outFlow = 0.
        node.flow = []
        node.flow_direction = []
        node.delta_pressure = []
        
        for edge in node.edges:        
            reverse = True
            if edge.at_start_node(node.index):
                reverse = False

            pressure = edge.get_scalar('Pressure',reverse=reverse)
            flow = edge.get_scalar('Flow',reverse=reverse)
            radius = edge.get_scalar('Radii',reverse=reverse)
            
            delta_pressure = (pressure[-1]-pressure[0])
            if delta_pressure>0.:
                # Pressure increasing along vessel - inflow
                flowDir = -1
                eFlow = flow[0] #-np.abs(flow[0])
                node.inFlow += np.abs(flow[0])
            elif delta_pressure<0.:
                # Pressure decreasing along vessel - outflow
                flowDir = 1
                eFlow = flow[0] #np.abs(flow[0])
                node.outFlow += np.abs(flow[0])
            elif delta_pressure==0.:
                flowDir = 0
                eFlow = flow[0] #np.abs(flow[0])
                
            velocity = np.asarray([(f*self.fluxConversion)/(np.square(r)*np.pi) for f,r in zip(flow,radius)])
            length = np.zeros(edge.npoints-1)
            pts = edge.coordinates
            for i in range(edge.npoints-1):
                length[i] = np.linalg.norm(pts[i+1]-pts[i]) # Length (um)
    
            delay = np.asarray([np.sum(l/np.abs(v)) if v!=0. else 1.e6 for l,v in zip(length,velocity)]) # seconds
                
            node.flow_direction.append(flowDir)
            node.flow.append(eFlow)
            node.delta_pressure.append(delta_pressure)
            
            edge.velocity = velocity
            edge.length = length
            edge.delay = delay
            edge.concentration = np.zeros([edge.npoints,len(self.time)])
            
        if 0 in node.flow_direction:
            # Find combination that minimises inflow/outflow difference
            inFlow = [node.flow[i] for i,dir in enumerate(node.flow_direction) if dir>0]
            outFlow = [node.flow[i] for i,dir in enumerate(node.flow_direction) if dir<0]
            noPD = [node.flow[i] for i,dir in enumerate(node.flow_direction) if dir==0]
            noPDInds = [i for i,dir in enumerate(node.flow_direction) if dir==0]
            if len(noPD)>2:
                #print('Number of no pressure-drop branches > 2! {}'.format(len(noPD)))
                pass
            elif len(noPD)==2 and noPD[0]==noPD[1] and node.outFlow==node.inFlow:
                #print('Number of no pressure-drop branches = 2 and are equal flow {}'.format(len(noPD)))
                node.flow_direction[noPDInds[0]] = -1
                node.flow_direction[noPDInds[1]] = 1
                node.inFlow = np.sum(inFlow)+np.sum(noPD[0])
                node.outFlow = np.sum(outFlow)+np.sum(noPD[1])
            else:
                if np.abs(np.sum(inFlow)-np.sum(outFlow)-np.sum(noPD))<0.01:
                    # Zero pressure-drop edges are outflows (+1 direction)
                    for ind in noPDInds:
                        node.flow_direction[ind] = 1
                    node.outFlow = np.sum(outFlow)+np.sum(noPD)
                elif np.abs(np.sum(inFlow)-np.sum(outFlow)+np.sum(noPD))<0.01:
                    # Zero pressure-drop edges are inflows (-1 direction)
                    for ind in noPDInds:
                        node.flow_direction[ind] = -1
                    node.inFlow = np.sum(inFlow)+np.sum(noPD)
                else:
                    #print('No pressure drop edges cannot be balanced!')
                    pass
            
        if node.inFlow==0. and node.outFlow>0:
            node.is_inlet = True
        else:
            node.is_inlet = False
            
#        if node.nconn==1:
        if node.inFlow==0:
            node.inFlow = node.outFlow
        elif node.outFlow==0:
            node.outFlow = node.inFlow
            
        node.Q = np.zeros(node.nconn,dtype='float') + 0.
        for i,edge in enumerate(node.edges):
            if node.flow_direction[i]>=0:
                if node.inFlow!=0.:
                    #edge.Q = node.flow[i] / node.inFlow
                    node.Q[i] = node.flow[i] / node.inFlow
                    if node.Q[i]>1.001 and 0 in node.flow_direction:
                        
                        #import pdb
                        #pdb.set_trace()
                        pass
                else:
                    node.Q[i] = 0.
                    #edge.Q = 0.
            else:
                #edge.Q = 0.
                node.Q[i] = 0.
                
    def parker(self,t,delay):
        # PARKER (average arterial input function from human population)
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
            
            s = np.where(tMin<0.)
            if len(s[0])>0:
                conc[s[0]] = 0.
                
            if not np.all(np.isfinite(conc)):
                conc[:] = 0.
        except:
            import pdb
            pdb.set_trace()
            

        return conc
        
    def auc(self,edges):
        
        for edge in edges:
            auc = np.sum(edge.concentration,axis=1)*self.dt
            edge.add_scalar('AUC',auc)
            
    def inflow_ratio(self,nodeList):
        
        for node in nodeList:
            if node.outFlow>0.:
                r = node.inFlow / node.outFlow
            else:
                r = 0.
            node.add_scalar('Flow ratio',r)
            
    def get_concentration(self,edges):
        
        conc = np.zeros([len(edges),edges[0].concentration.shape[0]])        
        for i,edge in enumerate(edges):
            conc[i,:] = edge.concentration
            
    def save_graph(self,output_directory=None):
        
        if self.graph is None or self.nodeList is None:
            return

        # Calculate AUC        
        edges = self.graph.edges_from_node_list(self.nodeList)
        self.auc(edges)
        
        new_graph = self.graph.node_list_to_graph(self.nodeList)
        if output_directory is not None:
            ofile = os.path.join(output_directory,'')+'ct_output.am'
            print('Writing graph to file ({})...'.format(ofile))
            new_graph.write(ofile)
            print('Writing complete')
            
            #ofile = os.path.join(output_directory,'')+'ct.p'            
            #print('Saving concentration data ({})...'.format(ofile))
            #conc = self.get_concentration(edges)
            #with open(ofile, 'wb') as handle:
            #    pickle.dump(ofile, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def inject(self, graph, output_directory=None):

        self.graph = graph        
        nodeCoords = graph.get_data('VertexCoordinates')
        edgeConn = graph.get_data('EdgeConnectivity')
        nedgepoints = graph.get_data('NumEdgePoints')
        edgeCoords = graph.get_data('EdgePointCoordinates')
        radius = graph.get_data('Radii')
        flow = graph.get_data('Flow')
        pressure = graph.get_data('Pressure')
        
        # Generate node list
        nodeList = graph.node_list()
        self.nodeList = nodeList
        for node in nodeList:
            self.vertex_flow_ordering(node)
        
        nconn = np.asarray([x.nconn for x in nodeList])
        
        sTerminal = np.where(nconn==1)[0]
        nTerminal = len(sTerminal)
        
        sInlet = [node.index for node in nodeList if node.is_inlet==True]# and node.flow_direction[0]>=0]
        nInlet = len(sInlet)
        total_inflow = np.abs(np.sum([nodeList[i].flow[0] for i in sInlet]))
        inletNodes = [nodeList[x] for x in sInlet]
        for inletNode in inletNodes:
            inletNode.inletQ = inletNode.flow[0] / total_inflow
            inletNode.inletDelay = 0.
            
        inletQ = [n.outFlow for n in inletNodes]
        mxQind = np.argmax(inletQ)
        inletNodes = [inletNodes[mxQind]]
        
        #self.inflow_ratio(nodeList)
        #new_graph = graph.node_list_to_graph(nodeList)
        #new_graph.write(output_directory+'flow_ratio.am')

        for inletCount,inletNode in enumerate(inletNodes):
            #import pdb
            #pdb.set_trace()
            print('Inlet: {}'.format(inletNode.index))
            visited = []
            visited_edges = []
            curNode = inletNode
            
            front = frontPKG.Front(inletNode,delay=inletNode.inletDelay,Q=inletNode.inletQ)
            
            # START WALK...
            endloop = False
            while endloop is False:
                
                if len(front.Q)>0:
                    mnQ = np.min(front.Q)
                    print('front size: {}, min Q: {}'.format(front.front_size,mnQ))

                if front.front_size>0:              
                    (current_nodes,delay,Q) = front.get_current_front()
                    
                    for n,curNode in enumerate(current_nodes):
                        
                        visited.append(curNode.index)
                        
                        res = [(nodeList[nodeIndex],curNode.edges[i],curNode.Q[i]) for i,nodeIndex in enumerate(curNode.connecting_node) if curNode.flow[i]>0.]
                                #if curNode.edges[i].index not in visited_edges and curNode[i].flow>=0.]
                        if len(res)>0:
                            # Choose which branch
                            next_nodes = [r[0] for r in res if r[2]!=0.]
                            via_edges = [r[1] for r in res if r[2]!=0.]
                            edge_Q = [r[2] for r in res if r[2]!=0.]
                            
                            delay_from = []
                            Q_from = []
                            for ve,via_edge in enumerate(via_edges):
                                visited_edges.append(via_edge.index)
                                try:
                                    conc = Q[n]*edge_Q[ve]*self.parker(self.time,delay[n])
                                    if np.all(conc==0):
                                        import pdb
                                        pdb.set_trace()
                                    via_edge.concentration += np.repeat([conc],via_edge.npoints,axis=0)
                                except:
                                    import pdb
                                    pdb.set_trace()
                        
                                delay_from.append(np.sum(via_edge.delay))
                                Q_from.append(Q[n]*edge_Q[ve])
                                
                            # Eliminate nodes that have a Q lower than limit
                            inds = [i for i,q in enumerate(Q_from) if q>self.Q_limit]
                            if len(inds)>0:
                                front.step_front([next_nodes[i] for i in inds],
                                                 delay=[delay_from[i] for i in inds],
                                                 Q=[Q_from[i] for i in inds])
                            else:
                                pass
                                #print('Exit Q: {}'.format(Q_from))
                        else:
                            #endloop = True
                            pass

                    front.complete_step()
                else:
                    endloop = True
                    print('Exit step {} - front size 0'.format(front.nstep))
                    break
                
            self.save_graph(output_directory=output_directory)
            
dir_ = 'C:\\Users\\simon\\Dropbox\\160113_paul_simulation_results\\LS147T\\1\\'
#pDir = dir+'paths\'
f = dir_+'spatialGraph_RIN.am'
graph = spatialgraph.SpatialGraph()
graph.read(f)
sfile = dir_+'ct_var_output.sav'
ia = InjectAgent()
try:
    ia.inject(graph,output_directory=dir_)
except KeyboardInterrupt:
    print('Ctrl-C interrupt! Saving graph')
    ia.save_graph(output_directory=dir_)
#except Exception,e:
#    import pdb
#    pdb.set_trace()