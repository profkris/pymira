import csv
import numpy as np
import spatialgraph

def csv2amira(filepath, ofile=''):
    start_coords, end_coords, radii = [], [], []
    with open(filepath, 'r') as fh:
        reader = csv.reader(fh)
        for row in reader:
            start_coords.append(np.asarray([float(x) for x in row[0:3]]))
            end_coords.append(np.asarray([float(x) for x in row[3:6]]))
            radii.append(float(row[-1]))
            
    start_coords = np.asarray(start_coords)
    end_coords = np.asarray(end_coords)
    radii = np.asarray(radii)
    
    nseg = start_coords.shape[0]
    
    # Infer connectivity from uniqueness of nodes
    all_coords = np.vstack([start_coords,end_coords])
    node = np.zeros(all_coords.shape[0],dtype='int') - 1
    node_coords = []
    node_count = 0
    for i,crd in enumerate(all_coords):
        if node[i]==-1:
            inds = np.where((all_coords[:,0]==crd[0]) & (all_coords[:,1]==crd[1]) & (all_coords[:,2]==crd[2]))
            if len(inds[0])>0:
                node[inds[0]] = node_count
                node_count += 1
                node_coords.append(crd)
                
    node_coords = np.asarray(node_coords)
    start_node_inds = node[:nseg]
    end_node_inds = node[nseg:]
    
    conns = np.zeros([nseg,2],dtype='int')
    points = np.zeros([nseg*2,3],dtype=start_coords.dtype)
    point_radii = np.zeros(nseg*2)
    npoints = np.zeros(nseg,dtype='int') + 2
    for i in range(nseg):
        conns[i] = [start_node_inds[i],end_node_inds[i]]
        points[2*i] = node_coords[start_node_inds[i]]
        points[2*i+1] = node_coords[end_node_inds[i]]
    
    graph = spatialgraph.SpatialGraph(initialise=True,scalars=['Radii'])
    graph.set_definition_size('VERTEX',node_count)
    graph.set_definition_size('EDGE',nseg)
    graph.set_definition_size('POINT',nseg*2)
    graph.set_data(node_coords,name='VertexCoordinates')
    graph.set_data(conns,name='EdgeConnectivity')
    graph.set_data(npoints,name='NumEdgePoints')
    graph.set_data(points,name='EdgePointCoordinates')
    graph.set_data(point_radii,name='Radii')

    if True:
        graph.sanity_check(deep=True)

    graph.write(ofile)
            
if __name__=='__main__':
    filepath = 'tests/test.csv'
    opath = 'tests/test.am'
    csv2amira(filepath, ofile=opath)
