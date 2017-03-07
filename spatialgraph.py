# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 11:49:52 2016

@author: simon
"""

from pymira import amiramesh
import numpy as np

class SpatialGraph(amiramesh.AmiraMesh): #
    
    def __init__(self,header_from=None,initialise=False,scalars=[],node_scalars=[]):
        amiramesh.AmiraMesh.__init__(self)
        if header_from is not None:
            import copy
            self.parameters = copy.deepcopy(header_from.parameters)
            self.definitions = copy.deepcopy(header_from.definitions)
            self.header = copy.deepcopy(header_from.header)
            self.fieldNames = copy.deepcopy(header_from.fieldNames)
        if initialise:
            self.initialise(scalars=scalars,node_scalars=node_scalars)
            
    def initialise(self,scalars=None,node_scalars=None):
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
        
        if len(scalars)>0:
            if type(scalars) is not list:
                scalars = [scalars]
            for i,sc in enumerate(scalars):
                self.add_field(name=sc,marker='@{}'.format(offset),
                                  definition='POINT',type='float',
                                  nelements=1,nentries=[0])
                offset = len(self.fields) + 1
                                  
        if len(node_scalars)>0:
            if type(node_scalars) is not list:
                node_scalars = [node_scalars]
            for i,sc in enumerate(node_scalars):
                self.add_field(name=sc,marker='@{}'.format(i+offset),
                                  definition='VERTEX',type='float',
                                  nelements=1,nentries=[0])
                offset = len(self.fields) + 1
                              
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
           #ntmp1 = len(np.where(conn[:,0]==i)[0])
           #ntmp2 = len(np.where(conn[:,1]==i)[0])
           ntmp1 = len([j for j,x in enumerate(conn) if i in x])
           nConn[i] = ntmp1 #+ ntmp2
           
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
        
        assert len(xrange)==2
        assert len(yrange)==2
        assert len(zrange)==2

        if not no_copy:        
            graph = self.clone()
        else:
            graph = self

        nodeCoords = self.get_data('VertexCoordinates')
        nnode = len(nodeCoords)
        #nedge = len(edgeConn)
        #nedgepoint = len(edgeCoords)
        
        # Spatial extent of nodes
        r = self.node_spatial_extent()

        # Locate nodes outside of ranges
        if xrange[1] is None:
            xrange[1] = r[0][1]
        if yrange[1] is None:
            yrange[1] = r[1][1]
        if zrange[1] is None:
            zrange[1] = r[2][1]
        xrange = [np.max([r[0][0],xrange[0]]),np.min([r[0][1],xrange[1]])]
        yrange = [np.max([r[1][0],yrange[0]]),np.min([r[1][1],yrange[1]])]
        zrange = [np.max([r[2][0],zrange[0]]),np.min([r[2][1],zrange[1]])]
        
        # Mark which nodes to keep / delete
        keepNode = np.zeros(nnode,dtype='bool') + True
        for i in range(nnode):
            x,y,z = nodeCoords[i,:]
            if x<xrange[0] or x>xrange[1] or y<yrange[0] or y>yrange[1] or z<zrange[0] or z>zrange[1]:
                keepNode[i] = False
                
        nodes_to_delete = np.where(keepNode==False)
        nodes_to_keep = np.where(keepNode==True)
        if len(nodes_to_keep[0])==0:
            print('No nodes left!')
            return
        
        editor = Editor()
        return editor.delete_nodes(nodes_to_delete[0])
        
    def remove_field(self,fieldName):
        f = [(i,f) for (i,f) in enumerate(self.fields) if f['name']==fieldName]
        if f[0][1] is None:
            print('Could not locate requested field: {}'.format(fieldName))
            return
        _  = self.fields.pop(f[0][0])
        
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
                            print('Node coordinates ({}) do not match start of edge ({}) coordinates: {} {}'.format(node.index,e.index,e.start_node_coords,node.coords))
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

    def edges_from_node_list(self,nodeList):
        
        edges = [n.edges for n in nodeList]
        import itertools
        edges = np.asarray(list(set(itertools.chain(*edges)))) # Flatten and get unique entries
        return edges        

    def node_list_to_graph(self,nodeList):
        
        nodeCoords = np.asarray([n.coords for n in nodeList])
        nnode = nodeCoords.shape[0]
        
        edges = self.edges_from_node_list(nodeList)

        edgeConn = np.asarray([[x.start_node_index,x.end_node_index] for x in edges])
        edgeCoords = np.concatenate([x.coordinates for x in edges])
        nedgepoint = np.array([x.npoints for x in edges])
        
        scalarNames = edges[0].scalarNames
        scalarData = [x.scalars for x in edges]        
        scalars = []
        nscalar = len(scalarNames)
        for i in range(nscalar): 
            scalars.append(np.concatenate([s[i] for s in scalarData]))
        
        nodeScalarNames = nodeList[0].scalarNames
        nodeScalarData = np.asarray([x.scalars for x in nodeList])
        nnodescalar = len(nodeScalarNames)
        nodeScalars = np.zeros([nnodescalar,nnode])
        for i in range(nnodescalar):
            #nodeScalars.append(np.concatenate([s[i] for s in nodeScalarData]))
            nodeScalars[i,:] = nodeScalarData[i::nnodescalar,0]
        
        #import spatialgraph
        graph = SpatialGraph(initialise=True,scalars=scalarNames,node_scalars=nodeScalarNames)
        
        graph.set_definition_size('VERTEX',nodeCoords.shape[0])
        graph.set_definition_size('EDGE',edgeConn.shape[0])
        graph.set_definition_size('POINT',edgeCoords.shape[0])

        graph.set_data(nodeCoords,name='VertexCoordinates')
        graph.set_data(edgeConn,name='EdgeConnectivity')
        graph.set_data(nedgepoint,name='NumEdgePoints')
        graph.set_data(edgeCoords,name='EdgePointCoordinates')
        for i,s in enumerate(scalars):
            graph.set_data(s,name=scalarNames[i])
        for i,s in enumerate(nodeScalars):
            graph.set_data(s,name=nodeScalarNames[i])
        
        return graph

class Editor(object):
    
    def delete_nodes(self,graph,nodes_to_delete):
        
        nodeCoords = graph.get_data('VertexCoordinates')
        edgeConn = graph.get_data('EdgeConnectivity')
        nedgepoints = graph.get_data('NumEdgePoints')
        edgeCoords = graph.get_data('EdgePointCoordinates')
        
        nnode = len(nodeCoords)
        nedge = len(edgeConn)
        nedgepoint = len(edgeCoords)
        
        # Modify node coordinates
        #nodeCoords_ed = np.delete(nodeCoords,nodes_to_delete[0],axis=0)
        nodes_to_keep = [x for x in range(nnode) if x not in nodes_to_delete]
        nodeCoords_ed = np.asarray([nodeCoords[x] for x in nodes_to_keep])
    
        # Modify connected edges
        keepEdge = np.in1d(edgeConn, nodes_to_keep).reshape(edgeConn.shape)
        keepEdge = np.asarray([all(x) for x in keepEdge])
        edges_to_delete = np.where(keepEdge==False)
        edges_to_keep = np.where(keepEdge==True)
        #edgeConn_ed = np.delete(edgeConn,edges_to_delete[0],axis=0)
        edgeConn_ed = np.asarray([edgeConn[x] for x in edges_to_keep][0])

        # Offset edge indices to 0
        #unqNodeIndices = np.unique(edgeConn_ed)
        unqNodeIndices = nodes_to_keep
        nunq = len(unqNodeIndices)
        newInds = np.arange(nunq)            
        edgeConn_ed_ref = np.zeros(edgeConn_ed.shape,dtype='int') - 1
        edgeConn_was = np.zeros(edgeConn_ed.shape,dtype='int') - 1
        # Update edge indices
        for i,u in enumerate(unqNodeIndices):
            sInds = np.where(edgeConn_ed==u)
            #newIndex = newInds[np.where(nodes_to_keep==u)[0]][0]
            newIndex = newInds[i]
            edgeConn_ed_ref[sInds[0][:],sInds[1][:]] = newIndex #newInds[i]
            edgeConn_was[sInds[0][:],sInds[1][:]] = u
        edgeConn_ed = edgeConn_ed_ref

        # Modify edgepoint number array
        #nedgepoints_ed = np.delete(nedgepoints,edges_to_delete[0],axis=0)
        nedgepoints_ed = np.asarray([nedgepoints[x] for x in edges_to_keep][0])

        # Mark which edgepoints to keep / delete
        keepEdgePoint = np.zeros(nedgepoint,dtype='bool') + True
        for edgeIndex in edges_to_delete[0]:
            npoints = nedgepoints[edgeIndex]
            strt = np.sum(nedgepoints[0:edgeIndex])
            fin = strt + npoints
            keepEdgePoint[strt:fin] = False

        # Modify edgepoint coordinates
        edgepoints_to_delete = np.where(keepEdgePoint==False)
        edgepoints_to_keep = np.where(keepEdgePoint==True)
        #edgeCoords_ed = np.delete(edgeCoords,edgepoints_to_delete[0],axis=0)
        edgeCoords_ed = np.asarray([edgeCoords[x] for x in edgepoints_to_keep[0]])
        
        # Sanity check
        check_sanity = True
        if check_sanity:
            for i,ed in enumerate(edgeConn_ed):                
                try:
                    s_nc = nodeCoords_ed[ed[0]]
                    f_nc = nodeCoords_ed[ed[1]]
                    strt = np.sum(nedgepoints_ed[0:i])
                    fin = strt + nedgepoints_ed[i]
                    pts = edgeCoords_ed[strt:fin]
                    if any([x!=y for x,y in zip(pts[0],s_nc)]):
                        edOld = edgeConn_was[i]
                        s_nc_old = nodeCoords[edOld[0]]
                        f_nc_old = nodeCoords[edOld[1]]
                        print('Node index {}'.format(i))
                        print('Old start node: {}'.format(s_nc_old))
                        print('New start node: {}'.format(s_nc))
                        print('Start point: {}'.format(pts[0]))
                        print('Old end node: {}'.format(f_nc_old))
                        print('New end node: {}'.format(f_nc))
                        print('End point: {}'.format(pts[-1]))
                        import pdb
                        pdb.set_trace()
                    if any([x!=y for x,y in zip(pts[-1],f_nc)]):
                        edOld = edgeConn_was[i]
                        s_nc_old = nodeCoords[edOld[0]]
                        f_nc_old = nodeCoords[edOld[1]]
                        print('Node index {}'.format(i))
                        print('Old start node: {}'.format(s_nc_old))
                        print('New start node: {}'.format(s_nc))
                        print('Start point: {}'.format(pts[0]))
                        print('Old end node: {}'.format(f_nc_old))
                        print('New end node: {}'.format(f_nc))
                        print('End point: {}'.format(pts[-1]))
                        import pdb
                        pdb.set_trace()
                except Exception,e:
                    import pdb
                    pdb.set_trace()
        
        # Update VERTEX definition
        vertex_def = graph.get_definition('VERTEX')
        vertex_def['size'] = [len(nodes_to_keep)]
        # Update EDGE definition
        edge_def = graph.get_definition('EDGE')
        edge_def['size'] = [len(edges_to_keep[0])]
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
            
        graph.set_graph_sizes()
        return graph
    
    def remove_intermediate_nodes(self,graph,file=None,nodeList=None):
        
        import pickle
        import os

        print 'Generating node list...'
        nodeList = graph.node_list()
        print 'Node list complete.'
        
        nnode = graph.nnode
        nedge = graph.nedge        
        nconn = np.array([node.nconn for node in nodeList])
        new_nodeList = []
        new_edgeList = []
        
        # Initialise list for mapping old to new node indices
        node_index_lookup = np.zeros(nnode,dtype='int') - 1
        edge_index_lookup = np.zeros(nedge,dtype='int') - 1
        # Mark if a node has become an edge point
        node_now_edge = np.zeros(nnode,dtype='int')
        node_converted = np.zeros(nnode,dtype='int')
        node_edges_checked = np.zeros(nnode,dtype='int')
        edge_converted = np.zeros(nedge,dtype='int')
        
        newNodeCount = 0
        newEdgeCount = 0
        
        for cntr,node in enumerate(nodeList):
            
            if node.index==3754:
                import pdb
                pdb.set_trace()
            
            # Is the current node branching (or terminal)?
            if (node.nconn==1 or node.nconn>2) and node_now_edge[node.index]==0 and node_edges_checked[node.index]==0:
                # If so, make a new node object
                if node_converted[node.index]==0:
                    print('NODE (START) {} {} {}'.format(newNodeCount,node.index,node.nconn))
                    newNode = Node(index=newNodeCount,coordinates=node.coords,connecting_node=[],old_index=node.index)
                    # Mark node as having been converted to a new node (rather than an edge)
                    node_converted[node.index] = 1
                    new_nodeList.append(newNode)
                    newNodeIndex = newNodeCount
                    node_index_lookup[node.index] = newNodeIndex
                    newNodeCount += 1
                else:
                    print('NODE (START, REVISITED) {} {} {}'.format(newNodeCount,node.index,node.nconn))
                    ind = node_index_lookup[node.index]
                    if ind<0:
                        import pdb
                        pdb.set_trace()
                    newNode = new_nodeList[ind]
                    newNodeIndex = newNode.index
                    
                node_edges_checked[node.index] = 1
                
                edges_complete = np.zeros(node.nconn,dtype='bool') + False
                
                # Loop through each branch
                for node_counter,connecting_node_index in enumerate(node.connecting_node):

                    # Initialise variables                    
                    curNodeIndex = connecting_node_index
                    endNode = None
                    visited = [node.index]
                    visited_edges = []

                    # Compile connecting edges -
                    connecting_edge = [e for x,e in zip(node.connecting_node,node.edges) if x==connecting_node_index]
                    
                    if connecting_edge[0].index==4098:
                        import pdb
                        pdb.set_trace()
                        
                    for connEdge in connecting_edge:

                        # Check if edge has already been converted (e.g. the return of a loop)
                        if edge_converted[connEdge.index]==0:
                            # Check whether to reverse coordinates in edge
                            if connEdge.end_node_index==node.index:
                                reverse_edge_indices = True
                                ecoords = connEdge.coordinates
                                ecoords = ecoords[::-1,:]
                                scalars = [s[::-1] for s in connEdge.scalars]
                            elif connEdge.start_node_index==node.index:
                                reverse_edge_indices = False
                                ecoords = connEdge.coordinates
                                scalars = connEdge.scalars
                            else:
                                import pdb
                                pdb.set_trace()
                                
                            # Create edge object to add points to during walk
                            print('EDGE {}'.format(newEdgeCount))
                            newEdge = Edge(index=newEdgeCount,start_node_index=newNode.index,
                                               start_node_coords=newNode.coords,
                                               coordinates=ecoords,
                                               npoints=ecoords.shape[0],
                                               scalars=scalars,
                                               scalarNames=connEdge.scalarNames)
                                               
                            new_edgeList.append(newEdge)
                            assert len(new_edgeList)==newEdgeCount+1
                            
                            edge_index_lookup[connEdge.index] = newEdgeCount
                            visited_edges.append(connEdge.index)
                            edge_converted[connEdge.index] = 1
                            
                            newEdgeCount += 1
                        
                            # Start walking - complete when a branching node is encountered
                            endFlag = False
                            
                            while endFlag is False:
                                curNode = nodeList[curNodeIndex]
                                visited.append(curNodeIndex)
                                
                                # If it's an intermediate (connecting) node
                                if curNode.nconn==2:
                                    # Check which connecting nodes have been visited already
                                    next_node_index = [x for x in curNode.connecting_node if x not in visited]
                                    # Get connecting edge (connected to unvisited, unconverted node)
                                    connecting_edge_walk = [e for x,e in zip(curNode.connecting_node,curNode.edges) if x not in visited and edge_converted[e.index]==0 ]
                                    # If no unvisited nodes have been identified...
                                    if len(connecting_edge_walk)==0:
                                        # Look for branching nodes that have been visited (i.e. loops)
                                        connecting_edge_walk = [e for x,e in zip(curNode.connecting_node,curNode.edges) if edge_converted[e.index]==0 ]
                                        if len(connecting_edge_walk)==1:
                                            foundConn = False
                                            # Check both start and end node indices
                                            for i,j in enumerate([connecting_edge_walk[0].start_node_index,connecting_edge_walk[0].end_node_index]):
                                                if nodeList[j].nconn > 2:
                                                    #Loop!
                                                    # Look for a connecting branch point
                                                    next_node_index = [j]
                                                    foundConn = True
                                            # If still nothing found...
                                            if not foundConn:
                                                import pdb
                                                pdb.set_trace()
                                                
                                    # If a connected edge has been found...
                                    if len(connecting_edge_walk)>0:
                                        # Check whether to reverse edge points
                                        if connecting_edge_walk[0].end_node_index==curNode.index:
                                            reverse_edge_indices = True
                                        elif connecting_edge_walk[0].start_node_index==curNode.index:
                                            reverse_edge_indices = False
                                        else:
                                            import pdb
                                            pdb.set_trace()
            
                                        # Add in connecting edge points
                                        # Reverse edge coordinates if necessary
                                        if reverse_edge_indices:
                                            scalars = [s[::-1] for s in connecting_edge_walk[0].scalars]
                                            #scalars = [s[1:-1] for s in scalars]
                                            coords = connecting_edge_walk[0].coordinates
                                            coords = coords[::-1,:]
                                            newEdge.add_point(coords,scalars=scalars,remove_last=True)
                                        else:
                                            scalars = connecting_edge_walk[0].scalars
                                            newEdge.add_point(connecting_edge_walk[0].coordinates,
                                                          scalars=scalars,remove_last=True)
        
                                        # Mark that node is now an edge point
                                        node_now_edge[curNodeIndex] = 1
            
                                        # If we've run out of nodes, then quit;
                                        # Otherwise, walk to the next node
                                        if len(next_node_index)==0:                                
                                            endFlag = True
                                        else:
                                            curNodeIndex = next_node_index[0]
                                            edge_converted[connecting_edge_walk[0].index] = 1
                                            edge_index_lookup[connecting_edge_walk[0].index] = newEdge.index
                                    else: # No connected edges found
                                        print('No connected edges...')
                                        endFlag = True
                                        
                                else: # Branch or terminal point
                                    endFlag = True
                                    end_node_index = curNode.index
                                    # Add in final edge coordinates, if necessary
                                    if not all([x==y for x,y in zip(newEdge.coordinates[-1,:],curNode.coords)]):
                                        # Reverse edge coordinates if necessary
                                        if connEdge.start_node_index!=curNode.index:
                                            scalars = [s[::-1] for s in connEdge.scalars]
                                            coords = connEdge.coordinates
                                            coords = coords[::-1,:]
                                            newEdge.add_point(coords,scalars=scalars,remove_last=True)
                                        else:
                                            scalars = connEdge.scalars
                                            newEdge.add_point(connEdge.coordinates,
                                                          scalars=scalars,remove_last=True)
                                    
                                # Sort out end nodes and edges
                                if endFlag:
                                    # Find end node
                                    if newEdge is None:
                                        import pdb
                                        pdb.set_trace()
                                    # If node has already been visited
                                    if node_converted[curNodeIndex]==1 and node_now_edge[curNodeIndex]==0:
                                        end_node_index_new = int(node_index_lookup[end_node_index])
                                        if end_node_index_new<0:
                                            import pdb
                                            pdb.set_trace()
                                        endNode = new_nodeList[end_node_index_new]
                                        #print('REVISITED NODE {} (END)'.format(endNode.index))
                                    # If node hasn't been converted, and isn't an edge
                                    elif node_now_edge[curNodeIndex]==0:
                                        print('NODE (END) {} {}'.format(newNodeCount,curNode.index))
                                        end_node_index_new = newNodeCount
                                        endNode = Node(index=end_node_index_new,coordinates=curNode.coords,connecting_node=[],old_index=curNode.index)
                                        node_converted[curNodeIndex] = 1
                                        new_nodeList.append(endNode) #[newNodeCount] = endNode
                                        node_index_lookup[end_node_index] = newNodeCount
                                        newNodeCount += 1
                                    else:
                                        import pdb
                                        pdb.set_trace()
                                        
                                    try:
                                        stat = newEdge.complete_edge(endNode.coords,end_node_index_new)
                                        if stat!=0:
                                            import pdb
                                            pdb.set_trace()
                                        print('EDGE COMPLETE: end node {}'.format(endNode.index))
                                    except Exception,e:
                                        print(e)
                                        import pdb
                                        pdb.set_trace()
                                        
                                    res = newNode.add_edge(newEdge)
                                    if not res:
                                        import pdb
                                        pdb.set_trace()
                                    if endNode.index!=newNode.index:
                                        res = endNode.add_edge(newEdge,reverse=True)
                                        if not res:
                                            import pdb
                                            pdb.set_trace()
                                    
                                    edges_complete[node_counter] = True
                                    
                                    break
                        else: # Edge has already been converted
                            newEdgeIndex = edge_index_lookup[connEdge.index]
                            if newEdgeIndex<0:
                                import pdb
                                pdb.set_trace()
                            newEdge = new_edgeList[newEdgeIndex]
                            if newEdge.start_node_index==newNode.index:
                                res = newNode.add_edge(newEdge)
                                if not res:
                                    print('Error: Edge {} is already attached to node {}'.format(newEdge.index,newNode.index))
                            elif newEdge.end_node_index==newNode.index:
                                res = newNode.add_edge(newEdge,reverse=True)
                                if not res:
                                    print('Error: Edge {} is already attached to node {}'.format(newEdge.index,newNode.index))
                            else:
                                import pdb
                                pdb.set_trace()
                            edges_complete[node_counter] = True
                        
#                    if edges_complete[node_counter]==False:
#                        import pdb
#                        pdb.set_trace()
                        
#                if newNode.nconn==2:
#                    import pdb
#                    pdb.set_trace()
#                if newNode.nconn!=node.nconn:
#                    import pdb
#                    pdb.set_trace()
#                    #assert endNode is not None
#                if not all(edges_complete):
#                    import pdb
#                    pdb.set_trace()

        #return new_nodeList
        se = np.where(edge_converted==0)
        elu = np.where(edge_index_lookup<0)
        incomplete_edges = [e for e in new_edgeList if e.complete is False]
        incomp = np.where(edge_converted==0)
        node2 = [n for n in new_nodeList if n.nconn==2]
        import pdb
        pdb.set_trace()
        new_nedge = newEdgeCount
        new_nnode = newNodeCount
        
        new_graph = graph.node_list_to_graph(new_nodeList)
        return new_graph

class Node(object):
    
    def __init__(self, graph=None, index=0, edge_indices=None, edge_indices_rev=None,
                 connecting_node=None, edges=None, coordinates=None, old_index=None,
                 scalars=None, scalarNames=None ):
                     
        self.index = index
        self.nconn = 0
        self.edge_indices = edge_indices
        self.edge_indices_rev = edge_indices_rev
        self.connecting_node = connecting_node
        self.edges = edges
        self.coords = coordinates
        self.old_index = old_index
        self.scalars = []
        self.scalarNames = []
        
        if graph is not None:
            vertCoords = graph.get_field('VertexCoordinates')['data']
            edgeConn = graph.get_field('EdgeConnectivity')['data']
            
            #s0 = [j for j,x in enumerate(edgeConn) if index in x]
            #s0 = np.where(edgeConn==index)
            #ns0 = len(s0)
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
            
        current_edge_indices = [e.index for e in self.edges]  
        if edge.index in current_edge_indices:
            return False
            
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
        self.nconn += 1
        return True
        
    def add_scalar(self,name,values):
        
        if name in self.scalarNames:
            print('Error: Node scalar field {} already exists!'.format(name))
            return
            
        if len(self.scalars)==0:
            self.scalars = [values]
            self.scalarNames = [name]
        else:
            self.scalars.append([values])
            self.scalarNames.append(name)
            
    def get_scalar(self,name):
        scalar = [x for i,x in enumerate(self.scalars) if self.scalarNames[i]==name]
        if len(scalar)==0:
            return None
        scalar[0]
            
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
            #self.end_node_index = edgeConn[index,1]
            #self.end_node_coords = nodeCoords[self.end_node_index,:]
            self.npoints = nedgepoints[index]
            self.coordinates = self.get_coordinates_from_graph(graph,index)
            self.scalars,self.scalarNames = self.get_scalars_from_graph(graph,index)
            self.complete_edge(nodeCoords[edgeConn[index,1],:],edgeConn[index,1])
            
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
            self.complete_edge(np.asarray(coords),end_index)
            
    def complete_edge(self,end_node_coords,end_node_index):
        stat = 0
        self.end_node_coords = np.asarray(end_node_coords)
        self.end_node_index = end_node_index
        self.complete = True
        
        if not all([x==y for x,y in zip(self.end_node_coords,self.coordinates[-1,:])]):
            print('Warning: End node coordinates do not match last edge coordiate!')
            stat = -1
        if not all([x==y for x,y in zip(self.start_node_coords,self.coordinates[0,:])]):
            print('Warning: Start node coordinates do not match first edge coordiate!')
            stat = -2
            
        return stat
            
    def at_start_node(self,index):
        if index==self.start_node_index:
            return True
        else:
            return False
            
    def add_scalar(self,name,values):
        
        # TO DO: add support for repeated scalars
        
        if len(values)!=self.npoints:
            print('Error: Scalar field has incorrect number of points')
            return
        if name in self.scalarNames:
            print('Error: Scalar field {} already exists!'.format(name))
            return
            
        if len(self.scalars)==0:
            self.scalars = values
            self.scalarNames = [name]
        else:
            self.scalars.append(values)
            self.scalarNames.append(name)
            
    def get_scalar(self,name,reverse=False):
        scalar = [x for i,x in enumerate(self.scalars) if self.scalarNames[i]==name]
        if len(scalar)==0:
            return None
        if reverse:
            return scalar[0][::-1]
        else:
            return scalar[0]
            
    def _print(self):
        print('EDGE ({})'.format(self.index))
        print('Number of points: {}'.format(self.npoints))
        print('Start node (index,coords): {} {}'.format(self.start_node_index,self.start_node_coords))
        print('End node (index,coords): {} {}'.format(self.end_node_index,self.end_node_coords))
        if self.scalarNames is not None:
            print('Scalar fields: {}'.format(self.scalarNames))
        if not self.complete:
            print('Incomplete...')