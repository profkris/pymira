# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 09:37:47 2017

@author: simon

Merge two spatial grphs
Required for converting Paul Sweeney's files

"""

import pymira.spatialgraph as sp
import os
import numpy as np
join = os.path.join
from tqdm import tqdm, trange

def merge_graphs(graph1,graph2):

    # Add fields from graph2 that aren't in graph1

    dif1  = list(set(graph1.fieldNames) - set(graph2.fieldNames))
    dif2  = list(set(graph2.fieldNames) - set(graph1.fieldNames))
    
    for fName in dif2:
        f = graph2.get_field(fName)
        marker = graph1.generate_next_marker()
        f['marker'] = marker
        #print(('Adding {} {}...'.format(marker,fName)))
        graph1.fields.append(f)
        
def combine_graphs(graph1,graph2):

    # Combine fields common to both graphs
    
    req_fields = ['VertexCoordinates', 'EdgePointCoordinates', 'EdgeConnectivity', 'NumEdgePoints']

    # Common fields
    fields = list(set(graph1.fieldNames).intersection(graph2.fieldNames))
    
    add_fields = list(set(fields) - set(req_fields))
    #breakpoint()
    
    for fName in req_fields:
        f1 = graph1.get_field(fName)
        f2 = graph2.get_field(fName)
        if fName=='VertexCoordinates':
            nnode1 = f1['data'].shape[0]
            nnode2 = f2['data'].shape[0]
        elif fName=='EdgePointCoordinates':
            npoints1 = f1['data'].shape[0]
            npoints2 = f2['data'].shape[0]
        elif fName=='EdgeConnectivity':
            nconn1 = f1['data'].shape[0]
            nconn2 = f2['data'].shape[0]
    
    for fName in fields:
        f1 = graph1.get_field(fName)
        f2 = graph2.get_field(fName)
        
        if fName=='EdgeConnectivity':
            # Offset all values by the number of nodes in graph 1
            data = np.concatenate([f1['data'],f2['data']+nnode1])
        else:
            data = np.concatenate([f1['data'],f2['data']])
        
        #print('Combining {}'.format(fName))
        graph1.set_data(data,name=fName)
    
    #print(nnode1,nnode2)
    graph1.set_definition_size('VERTEX',nnode1+nnode2)
    graph1.nnode = nnode1+nnode2
    graph1.set_definition_size('EDGE',nconn1+nconn2)
    graph1.nedge = nconn1+nconn2
    graph1.set_definition_size('POINT',npoints1+npoints2)
    graph1.nedgepoints = npoints1+npoints2
    return graph1
    
def combine_cco(path,mFiles,ofile):
    
    # Flag inlet / outlet nodes that shouldn't be connected
    ignore_node = np.zeros(len(mFiles),dtype='int')
    
    for i,f in enumerate(tqdm(mFiles)):
    
        if i==0:
            ignore_node[i] = 0
        else:
            ignore_node[i] = graph.nnode
 
        graph_to_add = sp.SpatialGraph()
        #print('Merging {}'.format(f))
        graph_to_add.read(join(path,f))
        
        marker = graph_to_add.generate_next_marker()
        if 'artery' in f:
            vesselType = np.zeros(graph_to_add.nedgepoint)
        elif 'vein' in f:
            vesselType = np.zeros(graph_to_add.nedgepoint) + 1
        if 'upper' in f:
            midLinePos = np.zeros(graph_to_add.nedgepoint)
        elif 'lower' in f:
            midLinePos = np.zeros(graph_to_add.nedgepoint) + 1
        else:
            midLinePos = np.zeros(graph_to_add.nedgepoint)
        marker = graph_to_add.generate_next_marker()
        graph_to_add.add_field(name='VesselType',marker=marker,definition='POINT',type='float',nelements=1,data=vesselType)
        marker = graph_to_add.generate_next_marker()
        graph_to_add.add_field(name='midLinePos',marker=marker,definition='POINT',type='float',nelements=1,data=midLinePos)
        
        if i>0:
            graph = combine_graphs(graph,graph_to_add)
        else:
            graph = graph_to_add
               
    # Connect artery / vein endpoints
    if False:
        verts = graph.get_data('VertexCoordinates')
        conn = graph.get_data('EdgeConnectivity')
        points = graph.get_data('EdgePointCoordinates')
        npoints = graph.get_data('NumEdgePoints')
        vType = graph.get_data('VesselType')
        mlp = graph.get_data('midLinePos') 
        radii = graph.get_data('Radii') 
        nnodeconn = graph.number_of_node_connections()
        a_endnodes = np.where(nnodeconn==1)
        def node_vtype(verts,nodeIndex,conns,npoints,vtype_points,edgeId=None):
            if edgeId is None:
                edgeIds = np.where((conns[:,0]==nodeIndex) | (conns[:,1]==nodeIndex))
                edgeId = edgeIds[0][0]
            npts = int(npoints[edgeId])
            x0 = int(np.sum(npoints[0:edgeId]))
            vtype = vtype_points[x0:x0+npts]
            pts = points[x0:x0+npts,:]
            node = verts[nodeIndex]
            if np.all(pts[0,:]==node):
                return vtype[0]
            else:
                return vtype[-1]
        vType_nodes = np.asarray([node_vtype(verts,i,conn,npoints,vType) for i in range(graph.nnode)])
        #mlp_nodes = np.asarray([node_vtype(verts,i,conn,npoints,mlp) for i in range(graph.nnode)])
        a_endpoints = np.asarray([i for i in range(graph.nnode) if vType_nodes[i]==0 and nnodeconn[i]==1 and i not in ignore_node])
        v_endpoints = np.asarray([i for i in range(graph.nnode) if vType_nodes[i]==1 and nnodeconn[i]==1 and i not in ignore_node])
        
        # Find closest endpoints and join them together
        from scipy.spatial.distance import cdist
        a_v = verts[a_endpoints]
        v_v = verts[v_endpoints]
        dists = cdist(a_v,v_v)
        
        # Find closest venous endpoint to each artery endpoint
        mn = np.argmin(dists,axis=1)
        mn_v = np.argmin(dists,axis=0)

        # Record where veins have been connected
        v_conn = np.zeros(v_v.shape[0])
        for ai in range(a_v.shape[0]):
            new_conn = np.asarray([a_endpoints[ai],v_endpoints[mn[ai]]])
            v_conn[mn[ai]] += 1
            conn = np.concatenate([conn,np.reshape(new_conn,[1,2])])
            new_edge = np.asarray([a_v[ai],v_v[mn[ai]]])
            points = np.concatenate([points,new_edge])
            npoints = np.concatenate([npoints,[2]])
            vType = np.concatenate([vType,[2,2]])
            radii = np.concatenate([radii,[2.5,2.5]])
            mlp = np.concatenate([mlp,[-1,-1]])
            
        for vi in range(v_v.shape[0]):
            if v_conn[vi]==0:
                new_conn = np.asarray([v_endpoints[vi],a_endpoints[mn_v[vi]]])
                v_conn[mn_v[vi]] += 1
                conn = np.concatenate([conn,np.reshape(new_conn,[1,2])])
                new_edge = np.asarray([v_v[vi],a_v[mn_v[vi]]])
                points = np.concatenate([points,new_edge])
                npoints = np.concatenate([npoints,[2]])
                vType = np.concatenate([vType,[2,2]])
                radii = np.concatenate([radii,[2.5,2.5]])
                mlp = np.concatenate([mlp,[-1,-1]])
            
        # Add new connections to graph
        graph.set_data(points,name='EdgePointCoordinates')
        graph.set_data(npoints,name='NumEdgePoints')
        graph.set_data(conn,name='EdgeConnectivity')
        graph.set_data(vType,name='VesselType')
        graph.set_data(radii,name='Radii')
        graph.set_data(mlp,name='midLinePos')
        graph.set_definition_size('POINT',points.shape[0])
        graph.set_definition_size('EDGE',conn.shape[0])

    graph.sanity_check()
    print('Combined graph written to {}'.format(join(path,ofile)))
    graph.write(join(path,ofile))
    
    return graph
    
def test_plot_cco():
    path = '/mnt/data2/retinasim/cco/graph'
    f = 'retina_cco.am'
    graph = sp.SpatialGraph()
    graph.read(join(path,f))
    graph.plot_graph()

if __name__=='__main__':
    path = '/mnt/data2/retinasim/cco/graph'
    mFiles = [  'retina_artery_upper_cco.csv.am',
                'retina_vein_upper_cco.csv.am',
                'retina_artery_lower_cco.csv.am',
                'retina_vein_lower_cco.csv.am',
             ]
    ofile = 'retina_cco.am'
    combine_cco(path,mFiles,ofile)
    test_plot_cco()
