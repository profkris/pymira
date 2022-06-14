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

def merge_graphs(graph1,graph2):

    # Add fields from graph2 that aren't in graph1

    dif1  = list(set(graph1.fieldNames) - set(graph2.fieldNames))
    dif2  = list(set(graph2.fieldNames) - set(graph1.fieldNames))
    
    for fName in dif2:
        f = graph2.get_field(fName)
        marker = graph1.generate_next_marker()
        f['marker'] = marker
        print(('Adding {} {}...'.format(marker,fName)))
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
        
        print('Combining {}'.format(fName))
        graph1.set_data(data,name=fName)
    
    print(nnode1,nnode2)
    graph1.set_definition_size('VERTEX',nnode1+nnode2)
    graph1.nnode = nnode1+nnode2
    graph1.set_definition_size('EDGE',nconn1+nconn2)
    graph1.nedge = nconn1+nconn2
    graph1.set_definition_size('POINT',npoints1+npoints2)
    graph1.nedgepoints = npoints1+npoints2
    return graph1
    
def combine_cco(path,mFiles,ofile):
    
    for i,f in enumerate(mFiles):
 
        graph_to_add = sp.SpatialGraph()
        print('Merging {}'.format(f))
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
        marker = graph_to_add.generate_next_marker()
        graph_to_add.add_field(name='VesselType',marker=marker,definition='POINT',type='float',nelements=1,data=vesselType)
        marker = graph_to_add.generate_next_marker()
        graph_to_add.add_field(name='midLinePos',marker=marker,definition='POINT',type='float',nelements=1,data=midLinePos)
        
        if i>0:
            graph = combine_graphs(graph,graph_to_add)
        else:
            graph = graph_to_add
        print(graph.nnode)
               
    # Connect artery / vein endpoints
    verts = graph.get_data('VertexCoordinates')
    conn = graph.get_data('EdgeConnectivity')
    points = graph.get_data('EdgePointCoordinates')
    npoints = graph.get_data('NumEdgePoints')
    vType = graph.get_data('VesselType')
    mlp = graph.get_data('midLinePos') 
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
    a_endpoints = np.asarray([i for i in range(graph.nnode) if vType_nodes[i]==0 and nnodeconn[i]==1])
    v_endpoints = np.asarray([i for i in range(graph.nnode) if vType_nodes[i]==1 and nnodeconn[i]==1])
    
    # Find closest endpoints and join them together
    from scipy.spatial.distance import cdist
    a_v = verts[a_endpoints]
    v_v = verts[v_endpoints]
    dists = cdist(a_v,v_v)

    breakpoint()   

    graph.sanity_check()
    print('Combined graph written to {}'.format(join(path,ofile)))
    graph.write(join(path,ofile))

if __name__=='__main__':
    path = '/mnt/data2/retinasim/cco/graph'
    mFiles = [  'retina_artery_upper_cco.csv.am',
                'retina_vein_upper_cco.csv.am',
                'retina_artery_lower_cco.csv.am',
                'retina_vein_lower_cco.csv.am',
             ]
    ofile = 'retina_cco.am'
    #combine_graphs.
    combine_cco(path,mFiles,ofile)
    #combine_cco()
