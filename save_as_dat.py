from pymira import spatialgraph as sp
import numpy as np
from tqdm import tqdm, trange
arr = np.asarray

"""
Converting to Paul Sweeney's data format for REANIMATE processing
Using variable names from his C++ code...
"""

def amira_to_dat(graph,ofile,network_name='anon'):

    vertexCoordinates = graph.get_data('VertexCoordinates')
    edgeConnectivity = graph.get_data('EdgeConnectivity')
    nedgePoints = graph.get_data('NumEdgePoints')
    edgePointCoordinates = graph.get_data('EdgePointCoordinates')
    
    rad_field = graph.get_radius_field()
    thickness = graph.point_scalars_to_edge_scalars(name=rad_field['name'])

    vessType = graph.point_scalars_to_edge_scalars(name='VesselType')
    if vessType is None:
        vessType = np.zeros(thickness.shape[0])

    nvertex = vertexCoordinates.shape[0]
    nedge = edgeConnectivity.shape[0]
    npoint = edgePointCoordinates.shape[0]
    #nseg = npoint - nedge
    nseg = nedge
    nnod = nvertex #+ npoint - 2*nedge

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

        flow = 1.
        diammax = -1.
        for i in range(nedge):
            segnodname1 = edgeConnectivity[i,0] + 1
            segnodname2 = edgeConnectivity[i,1] + 1
            diam = np.max([2.*thickness[i],4.])
            if diam > diammax:
                diammax = diam
                idiammax = i
                
            # Convert vessel type codes
            if vessType[i]==0: # artery
                vt = 1
            elif vessType[i]==2: # capillary
                vt = 2
            elif vessType[i]==1: # vein
                vt = 3

            handle.write(f"{i+1} {vt} {segnodname1} {segnodname2} {diam} {flow} {0.45}\n")
            
        print(f"Max diameter = {diammax}, idx = {idiammax}")

        handle.write(f"{nnod}   number of nodes\n")
        handle.write(f"Name    x       y       z\n")
        for i in range(nvertex):   # nodes from nvertex
            handle.write(f"{i+1} {vertexCoordinates[i,0]} {vertexCoordinates[i,1]} {vertexCoordinates[i,2]} {1.0}\n")

        if False:
            inod = nvertex
            for i in range(nedge):   # nodes from npoint
                x0 = np.sum(nedgePoints[:i])
                x1 = x0 + nedgePoints[i]
                pts = edgePointCoordinates[x0:x1]
                
                for j,p in ennumerate(pts): # Maybe need to limit this to exclude start and end points?
                    handle.write(f"{inod} {p[0]} {p[1]} {p[2]}\n")
                    inod += 1
                        
            if (cnt != nnod):
                print("*** Error: incorrect number of nodes")

        nodtype = graph.get_node_count()

        endnodes = np.where(nodtype==1)
        nnodbc = len(endnodes[0])
        handle.write(f"{nnodbc}   number of boundary nodes\n")
        handle.write(f"Name    bctyp     bcprfl     bchd\n")
        inod = 0
        consthd = 0.45
        #pressure = [80.,80.,10.,10.]
        pressure = [20.,18.] #,10.,10.,10.]
        #breakpoint()
        for inod,pr in zip(endnodes[0],pressure):
            handle.write(f"{inod+1} {0} {pr} {consthd}\n")

if __name__=='__main__':

    #fname = '/mnt/data2/retinasim/data/cco_circ_domain/graph/retina_cco_a2v.am'
    #fname = '/mnt/data2/retinasim/data/cco_circ_domain/graph/retina_cco_a2v_datprep.am'
    fname ='/mnt/data2/retinasim/data/cco_circ_domain/graph/test_network.am'
    graph = sp.SpatialGraph()
    graph.read(fname)
    
    if True:
        ed = sp.Editor()
        graph = ed.remove_intermediate_nodes(graph)
        graph.write(fname.replace('.am','_datprep.am'))
    
    #ofile = '/mnt/data2/retinasim/data/cco_circ_domain/graph/reanimate/retina_cco_datprep.dat'
    ofile = '/mnt/data2/retinasim/data/cco_circ_domain/graph/reanimate/test_network.dat'
    amira_to_dat(graph,ofile)
    print(f'Written {ofile}')


