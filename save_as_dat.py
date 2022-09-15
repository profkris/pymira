from pymira import spatialgraph as sp
import numpy as np
from tqdm import tqdm, trange
arr = np.asarray

"""
Converting to Paul Sweeney's data format for REANIMATE processing
Using variable names from his C++ code...
"""

def amira_to_dat(graph,ofile,thickness=None,network_name='anon'):

    vertexCoordinates = graph.get_data('VertexCoordinates')
    edgeConnectivity = graph.get_data('EdgeConnectivity')
    nedgePoints = graph.get_data('NumEdgePoints')
    edgePointCoordinates = graph.get_data('EdgePointCoordinates')
    
    if thickness is None:
        thickness = graph.get_radius_field()['data']
    
    nvertex = vertexCoordinates.shape[0]
    nedge = edgeConnectivity.shape[0]
    npoint = edgePointCoordinates.shape[0]
    nseg = npoint - nedge
    nnod = nvertex + npoint - 2*nedge

    maxEdge = np.max(edgePointCoordinates,axis=0)
    maxVertex = np.max(vertexCoordinates,axis=0)
    mx = np.max(np.vstack([maxEdge,maxVertex]),axis=0)
    minEdge = np.min(edgePointCoordinates,axis=0)
    minVertex = np.min(vertexCoordinates,axis=0)
    mn = np.min(np.vstack([minEdge,minVertex]),axis=0)
    
    alx,aly,alz = mx - mn

    with open(ofile,'w') as handle:
        handle.write(f"{network_name} network derived from Amira Spatial Graph\n")
        handle.write(f"{alx} {aly} {alz} box dimensions in microns - adjust as needed\n")
        handle.write(f"32 32 32 number of tissue points in x,y,z directions - adjust as needed\n")
        handle.write(f"10	outer bound distance - adjust as needed\n")
        handle.write(f"100	max. segment length - adjust as needed\n")
        handle.write(f"30	maximum number of segments per node - adjust as needed\n")
        handle.write(f"{nseg}	total number of segments\n")
        handle.write(f"SegName Type StartNode EndNode Diam   Flow[nl/min]    Hd\n");

        segnodname = np.zeros([2,nseg])
        k,cnt = 0,0
        diammax = -1.
        for i in range(nedge):
            segnodname1 = edgeConnectivity[i,0]
            segnodname2 = edgeConnectivity[i,1]
            diam = np.max([2.*thickness[k],4.])
            if diam > diammax:
                diammax = diam
                idiammax = k

            handle.write(f"{k} {3} {segnodname1} {segnodname2} {diam} {0.0, 0.45}\n")
            segnodname[0,cnt] = segnodname1
            segnodname[1,cnt] = segnodname2
            cnt += 1
            k += 1
            
        return
            
        segnodname = segnodname[:,0:cnt] #.submat(0,0,1,cnt-1);
        print(f"Max diameter = {diammax}, idx = {idiammax}")

        if cnt!=nseg:
            print("*** Error: incorrect number of segments ***")

        handle.write(f"{nnod}   number of nodes\n")
        handle.write(f"Name    x       y       z\n")
        cnt = 0
        
        for i in range(nvertex):   # nodes from nvertex
            handle.write(f"{i+1} {vertexCoordinates[i,0]} {vertexCoordinates[i,1]} {vertexCoordinates[i,2]}\n")
            cnt += 1

        k = -1
        inod = nvertex
        for i in range(nedge):   # nodes from npoint
            for j in range(1,nedgePoints[i]):
                k += 1
                inod += 1
                if j>1 and j<nedgePoints[i]:
                    handle.write(f"{inod} {edgePointCoordinates[k,0]} {edgePointCoordinates[k,1]} {edgePointCoordinates[k,2]}\n")
                    cnt += 1
                    
        if (cnt != nnod):
            print("*** Error: incorrect number of nodes")


if __name__=='__main__':

    fname = '/mnt/data2/retinasim/data/cco_circ_domain/graph/retina_cco.am'
    graph = sp.SpatialGraph()
    graph.read(fname)
    
    if True:
        ed = sp.Editor()
        graph = ed.remove_intermediate_nodes(graph)

        # Calculate point scalars from edge scalars (i.e. one value per edge)
        rad_field = graph.get_radius_field()
        sc = graph.point_scalars_to_edge_scalars(name=rad_field['name'])
    
    ofile = fname.replace('.am','.dat')
    breakpoint()
    amira_to_dat(graph,ofile,thickness=sc[0])


