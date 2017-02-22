# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 11:49:52 2016

@author: simon
"""

from pymira import amiramesh
import numpy as np

class SpatialGraph(amiramesh.AmiraMesh): #
    
    def __init__(self,header_from=None,initialise=False,scalars=[]):
        amiramesh.AmiraMesh.__init__(self)
        if header_from is not None:
            import copy
            self.parameters = copy.deepcopy(header_from.parameters)
            self.definitions = copy.deepcopy(header_from.definitions)
            self.header = copy.deepcopy(header_from.header)
            self.fieldNames = copy.deepcopy(header_from.fieldNames)
        if initialise:
            self.initialise(scalars=scalars)
            
    def initialise(self,scalars=[]):
        self.fileType = '3D ASCII 2.0'
        self.filename = ''
        
        self.add_definition('VERTEX',[0])
        self.add_definition('EDGE',[0])
        self.add_definition('POINT',[0])
        
        self.add_parameter('ContentType','HxSpatialGraph')

        self.add_field(name='VertexCoordinates',marker='@1',
                              definition='VERTEX',type='float',
                              nelements=3,nentries=[0])
        self.add_field(name='EdgeConnectivity',marker='@2',
                              definition='EDGE',type='int',
                              nelements=2,nentries=[0])
        self.add_field(name='NumEdgePoints',marker='@3',
                              definition='EDGE',type='int',
                              nelements=1,nentries=[0])
        self.add_field(name='EdgePointCoordinates',marker='@4',
                              definition='POINT',type='float',
                              nelements=3,nentries=[0])
                              
        offset = len(self.fields) + 1
        if type(scalars) is not list:
            scalars = [scalars]
        for i,sc in enumerate(scalars):
            self.add_field(name=sc,marker='@{}'.format(i+offset),
                              definition='POINT',type='float',
                              nelements=1,nentries=[0])
                              
        self.fieldNames = [x['name'] for x in self.fields]
        
    def read(self,*args,**kwargs):
        if not amiramesh.AmiraMesh.read(self,*args,**kwargs):
            return False
        if self.get_parameter_value("ContentType")!="HxSpatialGraph":
            print 'Warning: File is not an Amira SpatialGraph!'

        self.set_graph_sizes()
                
        return True
        
    def set_graph_sizes(self):
        try:
            self.nnode = self.get_definition('VERTEX')['size'][0]
        except:
            pass
        try:
            self.nedge = self.get_definition('EDGE')['size'][0]
        except:
            pass
        try:
            self.nedgepoint = self.get_definition('POINT')['size'][0]
        except:
            pass
    
    def reset_data(self):
        for x in self.fields:
            x['data'] = None
        for x in self.definitions:
            x['size'] = [0]
        for x in self.fields:
            x['shape'] = [0L,x['nelements']]
            x['nentries'] = [0L]
            
    def add_node(self,node=None,index=0,coordinates=[0.,0.,0.]):
        nodeCoords = self.get_field('VertexCoordinates')['data']
        if node is not None:
            coordinates = node.coords
            index = node.index
        if nodeCoords is not None:
            newData = np.vstack([nodeCoords, np.asarray(coordinates)])
            self.set_definition_size('VERTEX',newData.shape)
        else:
            newData = np.asarray(coordinates)
            self.set_definition_size('VERTEX',[1,newData.shape])
        self.set_data(newData,name='VertexCoordinates')
        
        nodeCoords = None # free memory (necessary?)
    
    def add_node_connection(self,startNode,endNode,edge):
        edgeConn = self.get_field('EdgeConnectivity')['data']
        nedgepoints = self.get_field('NumEdgePoints')['data']
        edgeCoords = self.get_field('EdgePointCoordinates')['data']
        
        # Add connection
        if edgeConn is not None:
            newData = np.vstack([edgeConn, np.asarray([startNode.index,endNode.index])])
            self.set_definition_size('EDGE',[1,newData.shape[0]])
        else:
            newData = np.asarray([startNode.index,endNode.index])
            self.set_definition_size('EDGE',newData.shape[0])
        self.set_data(newData,name='EdgeConnectivity')
        
        # Add number of edge points
        npoints = edge.coordinates.shape[0]
        if nedgepoints is not None:
            newData = np.vstack([nedgepoints, np.asarray(npoints)])
        else:
            newData = np.asarray(npoints)
        self.set_data(newData,name='NumEdgePoints')
        
        # Add edge point coordinates
        if edgeCoords is not None:
            newData = np.vstack([edgeCoords, np.asarray(edge.coordinates)])
        else:
            newData = np.asarray(edge.coordinates)
        self.set_definition_size('POINTS',newData.shape[0])
        self.set_data(newData,name='EdgePointCoordinates')

    def number_of_node_connections(self,file=None):
    
       #Identify terminal nodes
       nodeCoords = amdata.fields[0]['data']
       nnode = nodeCoords.shape[0]
       nConn = np.asarray([0]*nnode)
       conn = amdata.fields[1]['data']
       
       for i in range(nnode):
           ntmp1 = len(np.where(conn[:,0]==i)[0])
           ntmp2 = len(np.where(conn[:,1]==i)[0])
           nConn[i] = ntmp1 + ntmp2
           
       return nConn
       
    def node_list(self):
        # Convert graph to a list of node (and edge) objects
        nodeCoords = self.get_field('VertexCoordinates')['data']
        nnode = nodeCoords.shape[0]
        return [Node(graph=self,index=nodeIndex) for nodeIndex in range(nnode)]
        
    def clone(self):
        import copy
        return copy.deepcopy(self)
        
    def node_spatial_extent(self):
        
        nodecoords = self.get_data('VertexCoordinates')
        rx = [np.min(nodecoords[:,0]),np.max(nodecoords[:,0])]
        ry = [np.min(nodecoords[:,1]),np.max(nodecoords[:,1])]
        rz = [np.min(nodecoords[:,2]),np.max(nodecoords[:,2])]
        return [rx,ry,rz]
        
    def constrain_nodes(self,xrange=[None,None],yrange=[None,None],zrange=[None,None],no_copy=True):
        nodeCoords = self.get_data('VertexCoordinates')
        edgeConn = self.get_data('EdgeConnectivity')
        nedgepoints = self.get_data('NumEdgePoints')
        edgeCoords = self.get_data('EdgePointCoordinates')
        
        assert len(xrange)==2
        assert len(yrange)==2
        assert len(zrange)==2

        if not no_copy:        
            graph = self.clone()
        else:
            graph = self

        nnode = len(nodeCoords[:,0])
        nedge = len(edgeConn[:,0])
        nedgepoint = len(edgeCoords[:,0])
        
        r = self.node_spatial_extent()

        if xrange[1] is None:
            xrange[1] = r[0][1]
        if yrange[1] is None:
            yrange[1] = r[1][1]
        if zrange[1] is None:
            zrange[1] = r[2][1]
        xrange = [np.max([r[0][0],xrange[0]]),np.min([r[0][1],xrange[1]])]
        yrange = [np.max([r[1][0],yrange[0]]),np.min([r[1][1],yrange[1]])]
        zrange = [np.max([r[2][0],zrange[0]]),np.min([r[2][1],zrange[1]])]
        
        keepNode = np.zeros(nnode) + 1
        for i in range(nnode):
            x,y,z = nodeCoords[i,:]
            if x<xrange[0] or x>xrange[1] or y<yrange[0] or y>yrange[1] or z<zrange[0] or z>zrange[1]:
                keepNode[i] = 0
        
        # Delete nodes
        nodes_to_delete = np.where(keepNode==0)
        nodes_to_keep = np.where(keepNode==1)
        if len(nodes_to_keep[0])==0:
            print('No nodes left!')
            return
        #graph.remove_nodes(nodes_to_delete[0])
        nodeCoords_ed = np.delete(nodeCoords,nodes_to_delete[0],axis=0)
        # Update VERTEX definition
        vertex_def = graph.get_definition('VERTEX')
        vertex_def['size'] = [len(nodes_to_keep[0])]
        
        # Delete connected edges
        keepEdge = np.zeros(nedge) + 1
        for index,e in enumerate(nodes_to_delete[0]):
            s0 = np.where(edgeConn[:,0]==e)
            s1 = np.where(edgeConn[:,1]==e)
            if len(s0[0])>0:
                keepEdge[s0[0]] = 0
            if len(s1[0])>0:
                keepEdge[s1[0]] = 0
        edges_to_delete = np.where(keepEdge==0)
        edges_to_keep = np.where(keepEdge==1)
        edgeConn_ed = np.delete(edgeConn,edges_to_delete[0],axis=0)
        #import pdb
        #pdb.set_trace()
        # Offset edge indices to 0
        unq = np.unique(edgeConn_ed)
        nunq = len(unq)
        newInds = np.arange(nunq)
        # Update edge indices
        for i,u in enumerate(unq):
            sInds = np.where(edgeConn_ed==u)
            edgeConn_ed[sInds[0][:],sInds[1][:]] = newInds[i]
            #s0 = np.where(edgeConn_ed[:,0]==u)
            #edge_conn_ed_newInds[s0[0],0] = newInds[u]
        #edge_conn_ed_newInds = newInds[unq[edgeConn_ed]]
        nedgepoints_ed = np.delete(nedgepoints,edges_to_delete[0],axis=0)
        # Update EDGE definition
        edge_def = graph.get_definition('EDGE')
        edge_def['size'] = [len(edges_to_keep[0])]

        keepEdgePoint = np.zeros(nedgepoint) + 1
        for edgeIndex in edges_to_delete[0]:
            [strt,fin] = self.edgepoint_indices(edgeIndex)
            pts = edgeCoords[strt:fin,:]
            e = self.get_edge(edgeIndex)
            if e.end_node_index in nodes_to_keep[0] and \
               e.start_node_index in nodes_to_keep[0]:
                import pdb
                pdb.set_trace()
            keepEdgePoint[strt:fin] = 0
        edgepoints_to_delete = np.where(keepEdgePoint==0)
        edgepoints_to_keep = np.where(keepEdgePoint==1)
        edgeCoords_ed = np.delete(edgeCoords,edgepoints_to_delete[0],axis=0)
        # Update POINT definition
        edgepoint_def = graph.get_definition('POINT')
        edgepoint_def['size'] = [len(edgepoints_to_keep[0])]
        
        graph.set_data(nodeCoords_ed,name='VertexCoordinates')
        graph.set_data(edgeConn_ed,name='EdgeConnectivity')
        graph.set_data(nedgepoints_ed,name='NumEdgePoints')
        graph.set_data(edgeCoords_ed,name='EdgePointCoordinates')
        
        #Check for any other scalar fields
        scalars = [f for f in self.fields if f['definition'].lower()=='point' and len(f['shape'])==1]
        for sc in scalars:
            data_ed = np.delete(sc['data'],edgepoints_to_delete[0],axis=0)
            graph.set_data(data_ed,name=sc['name'])

        if not no_copy:                
            return graph
            
#    def remove_nodes(self,index):
#        nodeCoords = self.get_data('VertexCoordinates')
#        edgeConn = self.get_data('EdgeConnectivity')
#        nedgepoints = self.get_data('NumEdgePoints')
#        edgeCoords = self.get_data('EdgePointCoordinates')
#        
#        nnode = len(nodeCoords[:,0])
#        assert index>=0
#        assert index<nnode
#        
#        nodeCoords_ed = np.delete(nodeCoords,index,axis=0)
#        
#        s0 = np.where(edgeConn[:,0]==index)
#        s1 = np.where(edgeConn[:,1]==index)
#        if len(s0[0])>0:
#            for i,e in enumerate(s0[0]):
#                self.remove_edge(i)
#                
#    def remove_edge(self,index):
#        edgeConn = self.get_data('EdgeConnectivity')
#        nedgepoints = self.get_data('NumEdgePoints')
#        edgeCoords = self.get_data('EdgePointCoordinates')
        
    def remove_field(self,fieldName):
        f = [(i,f) for (i,f) in enumerate(self.fields) if f['name']==fieldName]
        if f[0][1] is None:
            print('Could not locate requested field: {}'.format(fieldName))
            return
        _  = self.fields.pop(f[0][0])
    
    def remove_intermediate_nodes(self,file=None,nodeList=None):
        
        import pickle
        import os
        import spatialgraph
        
        if nodeList is None:
            pFile = self.dir+'\\node_list.p'
        #if not os.path.isfile(pFile):
        #if True:
            print 'Generating node list...'
            nodeList = self.node_list()
            print 'Node list complete.'
            #print 'Saving node list to {}'.format(pFile)
            #with open(pFile,'wb') as f:
            #    pickle.dump(nodeList,f,protocol=pickle.HIGHEST_PROTOCOL)
            #print 'Node list saved.'
        #else:
        #    print 'Restoring node list...'
        #    with open(pFile,'rb') as f:
        #        nodeList = pickle.load(f)
        
        nnode = self.nnode
        nedge = self.nedge        
        nconn = np.array([node.nconn for node in nodeList])
        new_nodeList = []
        
        # Initialise list for mapping old to new node indices
        node_index_lookup = np.zeros(nnode) - 1.
        # Mark if a node has become an edge point
        node_now_edge = np.zeros(nnode)
        node_converted = np.zeros(nnode)
        edge_converted = np.zeros(nedge)
        
#        newNodeCoord = None
#        newEdgeConn = None
        
        newNodeCount = 0
        newEdgeCount = 0
        for cntr,node in enumerate(nodeList):
            print(' {} of {}. Now edge? {}'.format(cntr,len(nodeList),node_now_edge[node.index]))
            # Is the current node branching?
            if (node.nconn==1 or node.nconn>2) and node_now_edge[node.index]==0 and node_converted[node.index]==0:
                # If so, make a new node object
                print('NODE {}'.format(newNodeCount))
                newNode = spatialgraph.Node(index=newNodeCount,coordinates=node.coords,connecting_node=[])

                node_converted[node.index] = 1
                new_nodeList.append(newNode)
                
                newNodeIndex = newNodeCount
                node_index_lookup[node.index] = newNodeIndex
                newNodeCount += 1
  
                # Loop through each branch
                for connecting_node_index in node.connecting_node:
                    # Create an edge                    
                    curNodeIndex = connecting_node_index
                    endNode = None
                    visited = [node.index]
                    
                    connecting_edge = [e for x,e in zip(node.connecting_node,node.edges) if x==connecting_node_index]
                    if connecting_edge[0].start_node_index!=node.index:
                        reverse_edge_indices = True
                        ecoords = connecting_edge[0].coordinates
                        ecoords = ecoords[::-1,:]
                        #ecoords = ecoords[0:-1,:]
                        scalars = [s[::-1] for s in connecting_edge[0].scalars]
                        #scalars = [s[0:-1] for s in scalars]
                    else:
                        reverse_edge_indices = False
                        #ecoords = connecting_edge[0].coordinates[0:-1,:]
                        #scalars = [s[0:-1] for s in connecting_edge[0].scalars]
                        ecoords = connecting_edge[0].coordinates
                        scalars = connecting_edge[0].scalars
                    print('EDGE {}'.format(newEdgeCount))
                        
                    newEdge = Edge(index=newEdgeCount,start_node_index=node.index,
                                       start_node_coords=node.coords,
                                       coordinates=ecoords,
                                       npoints=ecoords.shape[0],
                                       scalars=scalars,
                                       scalarNames=connecting_edge[0].scalarNames)
                    newEdgeCount += 1
                    
                    #if newEdgeCount==2:
                    if True:
                        import pdb
                        pdb.set_trace()
                    
                    # Start walking - complete when a branching node is encountered
                    endFlag = False
                    #newEdge = None
                        
                    while endFlag is False:
                        #print curNodeIndex
                        curNode = nodeList[curNodeIndex]
                        visited.append(curNodeIndex)
                        
                        if len(connecting_edge)>1:
                            import pdb
                            pdb.set_trace()
                                       
                        # If it's an intermediate (connecting) node
                        if curNode.nconn==2:
                            next_node_index = [x for x in curNode.connecting_node if x not in visited]
                            connecting_edge = [e for x,e in zip(curNode.connecting_node,curNode.edges) 
                                               if x not in visited and edge_converted[e.index]==0 ]
                            if connecting_edge[0].start_node_index!=curNode.index:
                                reverse_edge_indices = True
                            else:
                                reverse_edge_indices = False

                            # Add in connecting edge points
                            if newEdgeCount==2:
                                import pdb
                                pdb.set_trace()
                            #if connecting_edge[0].npoints>2:
                            if True:
                                # Reverse edge coordinates if necessary
                                if reverse_edge_indices:
                                    scalars = [s[::-1] for s in connecting_edge[0].scalars]
                                    #scalars = [s[1:-1] for s in scalars]
                                    coords = connecting_edge[0].coordinates
                                    coords = coords[::-1,:]
                                    newEdge.add_point(coords,scalars=scalars,remove_last=True)
                                else:
                                    #scalars = [s[1:-1] for s in connecting_edge[0].scalars]
                                    scalars = connecting_edge[0].scalars
                                    newEdge.add_point(connecting_edge[0].coordinates,
                                                  scalars=scalars,remove_last=True)
                            elif connecting_edge[0].npoints==2:
                                # Reverse edge coordinates if necessary
                                if reverse_edge_indices:
                                    #scalars = [s[-1] for s in connecting_edge[0].scalars]
                                    scalars = [s[::-1] for s in connecting_edge[0].scalars]
                                    coords = connecting_edge[0].coordinates
                                    coords = coords[::-1,:]
                                    newEdge.add_point(connecting_edge[0].coordinates[-1],
                                                  scalars=scalars,remove_last=True)
                                else:
                                    scalars = [s[0] for s in connecting_edge[0].scalars]
                                    newEdge.add_point(connecting_edge[0].coordinates[0],
                                                  scalars=scalars,remove_last=True)
                            else:
                                import pdb
                                pdb.set_trace()

                            # Mark that node is now an edge point
                            node_now_edge[curNodeIndex] = 1

                            # If we've run out of nodes, then quit;
                            # Otherwise, walk to the next node
                            if len(next_node_index)==0:                                
                                endFlag = True
                            else:
                                curNodeIndex = next_node_index[0]
                                edge_converted[connecting_edge[0].index] = 1
                                
                        elif curNode.nconn==1:
                            endFlag = True # Terminal node
                            end_node_index = curNode.index
                        else: # Branch point
                            endFlag = True
                            end_node_index = curNode.index
                            
                        # Sort out end nodes and edges
                        if endFlag:
                            # Find end node
                            if newEdge is None:
                                import pdb
                                pdb.set_trace()
                            if node_converted[curNodeIndex]==1 and node_now_edge[curNodeIndex]==0:
                                end_node_index_new = int(node_index_lookup[end_node_index])
                                endNode = new_nodeList[end_node_index_new]
                            else:
                                print('NODE {}'.format(newNodeCount))
                                end_node_index_new = newNodeCount
                                endNode = spatialgraph.Node(index=end_node_index_new,coordinates=curNode.coords,connecting_node=[])
                                node_converted[curNodeIndex] = 1
                                new_nodeList.append(endNode) #[newNodeCount] = endNode
                                node_index_lookup[end_node_index] = newNodeCount
                                newNodeCount += 1
                                
#                            if connecting_edge[0].start_node_index!=curNode.index:
#                                reverse_edge_indices = True
#                            else:
#                                reverse_edge_indices = False
#
#                            if connecting_edge[0].start_node_index!=curNode.index:
#                                reverse_edge_indices = True
#                                scalars = [s[::-1] for s in connecting_edge[0].scalars]
#                                scalars = [s[0] for s in scalars]
#                                coords = connecting_edge[0].coordinates
#                                coords = coords[::-1,:]
#                                newEdge.add_point(coords[0,:],scalars=scalars,is_end=True,end_index=end_node_index_new)
#                            else:
#                                reverse_edge_indices = False
#                                scalars = [s[0] for s in connecting_edge[0].scalars]
#                                newEdge.add_point(connecting_edge[0].coordinates[0,:],
#                                                  scalars=scalars,is_end=True,end_index=end_node_index_new)

                            #scalars = [s[-1] for s in connecting_edge[0].scalars]
                            #newEdge.add_point(endNode.coords,is_end=True,end_index=end_node_index_new,
                            #                  scalars=scalars)
                            newEdge.complete_edge(endNode.coords,end_node_index_new)
                            newNode.add_edge(newEdge)
                            endNode.add_edge(newEdge,reverse=True)
                            
                            break
                        
                    assert endNode is not None
             
        #return new_nodeList
        new_nedge = newEdgeCount
        new_nnode = newNodeCount
        
        nodeCoords = np.asarray([n.coords for n in new_nodeList])
        edges = np.unique([f for f in n.edges for n in new_nodeList])
        edgeConn = np.asarray([[x.start_node_index,x.end_node_index] for x in edges])
        edgeCoords = np.concatenate([x.coordinates for x in edges])
        nedgepoint = np.array([x.npoints for x in edges])
        scalarNames = edges[0].scalarNames
        scalarData = [x.scalars for x in edges]
        scalars = []
        nscalar = len(scalarNames)
        for i in range(nscalar): 
            scalars.append(np.concatenate([s[i] for s in scalarData]))
        
        import spatialgraph
        new_graph = SpatialGraph(initialise=True,scalars=scalarNames)
        
        new_graph.set_definition_size('VERTEX',new_nnode)
        new_graph.set_definition_size('EDGE',new_nedge)
        new_graph.set_definition_size('POINT',edgeCoords.shape[0])

        new_graph.set_data(nodeCoords,name='VertexCoordinates')
        new_graph.set_data(edgeConn,name='EdgeConnectivity')
        new_graph.set_data(nedgepoint,name='NumEdgePoints')
        new_graph.set_data(edgeCoords,name='EdgePointCoordinates')
        for i,s in enumerate(scalars):
            new_graph.set_data(s,name=scalarNames[i])
        
        return new_graph
        
    def get_node(self,index):
        return Node(graph=self,index=index)
        
    def get_edge(self,index):
        return Edge(graph=self,index=index)
        
    def get_scalars(self):
        return [f for f in self.fields if f['definition'].lower()=='point' and len(f['shape'])==1]
        
    def edgepoint_indices(self,edgeIndex):
        nedgepoints = self.get_data('NumEdgePoints')
        edgeCoords = self.get_data('EdgePointCoordinates')
        nedge = len(nedgepoints)
        
        assert edgeIndex>=0
        assert edgeIndex<nedge
        
        npoints = nedgepoints[edgeIndex]
        if edgeIndex==0:
            start_index = 0
        else:
            start_index = np.sum(nedgepoints[0:edgeIndex])
        end_index = start_index + npoints
        
        return [start_index,end_index]
        
    def _print(self):
        print('GRAPH')
        print('Fields: {}'.format(self.fieldNames))
        for f in self.fields:
            print(f)
            
    def sanity_check(self,deep=False):
        
        self.set_graph_sizes()
        
        for d in self.definitions:
            defName = d['name']
            defSize = d['size'][0]
            fields = [f for f in self.fields if f['definition']==defName]
            for f in fields:
                if f['nentries'][0]!=defSize:
                    print('{} field size does not match {} definition size!'.format(f['name'],defName))
                if f['shape'][0]!=defSize:
                    print('{} shape size does not match {} definition size!'.format(f['name'],defName))
                if not all(x==y for x,y in zip(f['data'].shape,f['shape'])):
                    print('{} data shape does not match shape field!'.format(f['name']))

        if deep:
            for nodeInd in range(self.nnode):
                node = self.get_node(nodeInd)
                for i,e in enumerate(node.edges):
                    if not node.edge_indices_rev[i]:
                        if not all(x==y for x,y in zip(e.start_node_coords,node.coords)):
                            print('Node coordinates ({}) do not match start of edge ({}) coordinates'.format(node.index,e.index))
                        if not all(x==y for x,y in zip(e.coordinates[0,:],e.start_node_coords)):
                            print('Edge start point does not match edge/node start ({}) coordinates'.format(e.index))
                        if not all(x==y for x,y in zip(e.coordinates[-1,:],e.end_node_coords)):
                            print('Edge end point does not match edge/node end ({}) coordinates'.format(e.index))
                    else:
                        if not all(x==y for x,y in zip(e.end_node_coords,node.coords)):
                            print('Node coordinates ({}) do not match end of edge ({}) coordinates'.format(node.index,e.index))
                        if not all(x==y for x,y in zip(e.coordinates[0,:],e.start_node_coords)):
                            print('Edge end point does not match edge start (REVERSE) ({}) coordinates'.format(e.index))
                        if not all(x==y for x,y in zip(e.coordinates[-1,:],e.end_node_coords)):
                            print('Edge start point does not match edge end (REVERSE) ({}) coordinates'.format(e.index))
    
class Node(object):
    
    def __init__(self, graph=None, index=0, edge_indices=None, edge_indices_rev=None,
                 connecting_node=None, edges=None, coordinates=None ):
                     
        self.index = index
        self.nconn = 0
        self.edge_indices = edge_indices
        self.edge_indices_rev = edge_indices_rev
        self.connecting_node = connecting_node
        self.edges = edges
        self.coords = coordinates
        
        if graph is not None:
            vertCoords = graph.get_field('VertexCoordinates')['data']
            edgeConn = graph.get_field('EdgeConnectivity')['data']
            
            s0 = np.where(edgeConn[:,0]==index)
            ns0 = len(s0[0])
            s1 = np.where(edgeConn[:,1]==index)
            ns1 = len(s1[0])
            
            self.coords = vertCoords[index,:]
            self.nconn = ns0 + ns1
            
            self.edge_indices = edge_indices
            if self.edge_indices is None:
                self.edge_indices = []
            self.edge_indices_rev = []
            self.connecting_node = []
            self.edges = []
    
            for e in s0[0]:
                self.edge_indices.append(e)
                self.edges.append(Edge(graph=graph,index=e))
                self.edge_indices_rev.append(False)
                self.connecting_node.append(edgeConn[e,1])
            for e in s1[0]:
                self.edge_indices.append(e)
                self.edge_indices_rev.append(True)
                self.edges.append(Edge(graph=graph,index=e))
                self.connecting_node.append(edgeConn[e,0])
                
    def add_edge(self,edge,reverse=False):
        if self.edges is None:
            self.edges = []
        self.edges.append(edge)
        if self.edge_indices_rev is None:
            self.edge_indices_rev = []
        self.edge_indices_rev.append(reverse)
        if self.connecting_node is None:
            self.connecting_node = []
        if not reverse:
            self.connecting_node.append(edge.end_node_index)
        else:
            self.connecting_node.append(edge.start_node_index)
            
    def _print(self):
        print('NODE ({}):'.format(self.index))
        print('Coordinate: {}'.format(self.coords))
        print('Connected to: {}'.format(self.connecting_node))
        if len(self.connecting_node)>0:
            edgeInd = [e.index for e in self.edges]
            print('Connected via edges: {}'.format(edgeInd))
            
class Edge(object):
    
    def __init__(self, index=0, graph=None, 
                 start_node_index=None, start_node_coords=None,
                 end_node_index=None, end_node_coords=None,
                 npoints=0, coordinates=None, scalars=None,
                 scalarNames=None):
        self.index = index
        self.start_node_index = start_node_index
        self.start_node_coords = start_node_coords # numpy array
        self.end_node_index = end_node_index
        self.end_node_coords = end_node_coords # numpy array
        self.npoints = npoints
        self.coordinates = coordinates # numpy array
        self.complete = False
        self.scalars = scalars
        self.scalarNames = scalarNames
        
        if graph is not None:
            nodeCoords = graph.get_field('VertexCoordinates')['data']
            edgeConn = graph.get_field('EdgeConnectivity')['data']
            nedgepoints = graph.get_field('NumEdgePoints')['data']
            self.coordinates = self.get_coordinates_from_graph(graph,index)
            self.start_node_index = edgeConn[index,0]
            self.start_node_coords = nodeCoords[self.start_node_index,:]
            self.end_node_index = edgeConn[index,1]
            self.end_node_coords = nodeCoords[self.end_node_index,:]
            self.npoints = nedgepoints[index]
            self.coordinates = self.get_coordinates_from_graph(graph,index)
            self.complete = True
            self.scalars,self.scalarNames = self.get_scalars_from_graph(graph,index)
            
        #if self.coordinates is not None:
        #    assert self.npoints==len(self.coordinates[:,0])
        #    # Make sure start coordinates match
        #    assert all([x==y for x,y in zip(self.coordinates[0,:],self.start_node_coords)])
        #    # Make sure end coordinates match
        #    assert all([x==y for x,y in zip(self.coordinates[-1,:],self.end_node_coords)])
        
    def get_coordinates_from_graph(self,graph,index):
        nedgepoints = graph.get_field('NumEdgePoints')['data']
        coords = graph.get_field('EdgePointCoordinates')['data']
        nprev = np.sum(nedgepoints[0:index])
        ncur = nedgepoints[index]
        e_coords = coords[nprev:nprev+ncur,:]
        return e_coords
        
    def get_scalars_from_graph(self,graph,index):
        scalars = graph.get_scalars()
        if len(scalars)==0:
            return None,None
        nedgepoints = graph.get_field('NumEdgePoints')['data']
        nprev = np.sum(nedgepoints[0:index])
        ncur = nedgepoints[index]
        scalarData = []
        scalarNames = []
        for s in scalars:
            scalarData.append(s['data'][nprev:nprev+ncur])
            scalarNames.append(s['name'])
        return scalarData,scalarNames
        
    def add_point(self,coords,is_end=False,end_index=None,scalars=None,remove_last=False):
        coords = np.asarray(coords)
        #assert len(coords)==3
        if len(coords.shape)==2:
            npoints = coords.shape[0]
            p0 = coords[0]
        else:
            npoints = 1
            p0 = coords
        
        if self.coordinates is None:
            self.coordinates = []
            self.scalars = []
            self.scalarNames = []
            if self.start_node_coords is None:
                self.start_node_coords = np.asarray(p0)
            self.coordinates = np.asarray(coords)
            self.npoints = npoints
        else:
            if remove_last:
                if self.npoints>1:
                    self.coordinates = self.coordinates[0:-1,:]
                else:
                    self.coordinates = []
                self.npoints -= 1
            self.coordinates = np.vstack([self.coordinates, np.asarray(coords)])
            self.npoints += npoints
            if scalars is not None:
                if remove_last:
                    if self.npoints==0:
                        self.scalars = []
                    else:
                        self.scalars = [s[0:-1] for s in self.scalars]
                for i,sc in enumerate(scalars):
                    self.scalars[i] = np.append(self.scalars[i],scalars[i])
        if is_end:
            self.complete(np.asarray(coords),end_index)
            
    def complete_edge(self,end_node_coords,end_node_index):
        self.end_node_coords = np.asarray(end_node_coords)
        self.end_node_index = end_node_index
        self.complete = True
            
    def _print(self):
        print('EDGE ({})'.format(self.index))
        print('Number of points: {}'.format(self.npoints))
        print('Start node (index,coords): {} {}'.format(self.start_node_index,self.start_node_coords))
        print('End node (index,coords): {} {}'.format(self.end_node_index,self.end_node_coords))
        if self.scalarNames is not None:
            print('Scalar fields: {}'.format(self.scalarNames))
        if not self.complete:
            print('Incomplete...')