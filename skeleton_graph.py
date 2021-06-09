import numpy as np
arr = np.asarray
import os
from os.path import join
from pymira import spatialgraph
import nibabel as nib
from matplotlib import pyplot as plt
from scipy import ndimage
from skimage.morphology import skeletonize
from tqdm import trange

opath = '/home/simon/Desktop/Share/'

def shift_volume(X, dx, dy, dz):
    pad = [[np.abs(dx),np.abs(dx)],[np.abs(dy),np.abs(dy)],[np.abs(dz),np.abs(dz)]]
    X = np.pad(X.copy(),pad)
    X = np.roll(X, dy, axis=0)
    X = np.roll(X, dx, axis=1)
    X = np.roll(X, dz, axis=2)
    return X[pad[0][0]:-pad[0][1],pad[1][0]:-pad[1][1],pad[2][0]:-pad[2][1]]
    
def nearest_neighbours(vol,x,y,z):
    neigh = np.zeros(9+8+9) * np.nan
    dirs = []
    count = 0
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                if [i,j,k]!=[0,0,0]:
                    dirs.append([i,j,k])
                    xp,yp,zp = x+i, y+j, z+k
                    if xp>=0 and yp>=0 and zp>0 and xp<vol.shape[0] and yp<vol.shape[1] and zp<vol.shape[2]:
                        neigh[count] = vol[xp,yp,zp]
                    count += 1
    return neigh, np.asarray(dirs)

path = '/mnt/data2/Sahai/export/nifti'
fname = join(path,'HET7_vessel_ROI.nii')

#graph = spatialgraph.SpatialGraph()
#graph.read(fname)

print('Reading image...')
img = nib.load(fname)
data = img.get_fdata()
data[data>0] = 1
data = data.astype('uint8')

print('Skeletonising...')
skeleton = skeletonize(data) #, method='lee')
skel_conn = skeleton.copy()
for i in [-1,0,1]:
    for j in [-1,0,1]:
        for k in [-1,0,1]:
            if i!=0 and j!=0 and k!=0:
                skel_conn += shift_volume(skeleton,i,j,k)
skel_conn *= skeleton
print('Distance transform...')
dist_metric = ndimage.generate_binary_structure(3, 2)
#distmap,edges = ndimage.distance_transform_cdt(data,metric='cityblock',return_indices=True) # cityblock
resolution = arr([8.0,0.95,0.95])
distmap,edges = ndimage.distance_transform_edt(data,sampling=resolution,return_indices=True) # cityblock

tr = np.eye(4)
tr[0,0],tr[1,1],tr[2,2] = resolution
dimg = nib.Nifti1Image(distmap, tr)
ofile = join(opath,'distance_map.nii')
nib.save(dimg, ofile)
print('Wrote distance map to {}'.format(ofile))
dimg = nib.Nifti1Image(skel_conn, tr)
ofile = join(opath,'connectivity.nii')
nib.save(dimg, ofile)
print('Wrote connectivity map to {}'.format(ofile))

# Loop through all skeleton pixels looking for connections
inds = np.where(skel_conn>0)
edge_points = []
nedge_points = []
cur_edge_points = []
edge_points = []
conns = []
radius = []
visited = skeleton*0
# Get node coordinates
node_inds = np.where((skel_conn>1) | (skel_conn==1))
nodes = arr([arr([node_inds[0][i],node_inds[1][i],node_inds[2][i]]) for i in range(len(node_inds[0]))])

nnode = len(node_inds[0])
loop_count = 0
#while True:
if True:
    n_edge_start = len(edge_points)
    # Loop through each node (pixels with >2 or =1 connections)
    for i in trange(nnode):
        start_node_index = i
        # Node position
        x,y,z = node_inds[0][i],node_inds[1][i],node_inds[2][i]
        # Node index
        start_node_coords = [x,y,z]
        
        endloop = False
        # Temporary copy of visited pixel array
        visited_tmp = visited.copy()
        count = 0

        # Neibourhood of starting node
        neigh,dirs = nearest_neighbours(skel_conn,x,y,z)
        vis,dirs = nearest_neighbours(visited_tmp,x,y,z)
        # Locate neighbouring pixels that haven't been visited
        ninds = np.where((neigh>0) & (vis<=neigh))

        # If there are valid neighbouring pixels
        if len(ninds[0])>0:

            # Loop through each pixel neighbouring the start node
            for j in range(len(ninds[0])):
                # Starting poisition (the start node)
                xc,yc,zc = x,y,z
                # Next point location
                nind = ninds[0][j]
                xp,yp,zp = xc+dirs[nind][0], yc+dirs[nind][1], zc+dirs[nind][2]
                visited_tmp[xp,yp,zp] += 1
                # Initialise the edge array
                cur_edge_points = [[xc,yc,zc]]
                cur_radius = [distmap[xc,yc,zc]]

                # Start walking
                while True:
                    # Find valid neighbouring pixels 
                    cur_neigh,cur_dirs = nearest_neighbours(skel_conn,xp,yp,zp)
                    cur_vis,cur_dirs = nearest_neighbours(visited_tmp,xp,yp,zp)
                    ninds_cur = np.where((cur_neigh>0) & (cur_vis<=cur_neigh))
                    
                    if len(ninds_cur[0])>0:
                        # Only take the first index (valid?)
                        nind_cur = ninds_cur[0][0]
                        # Next point location
                        xp,yp,zp = xc+cur_dirs[nind_cur][0], yc+cur_dirs[nind_cur][1], zc+cur_dirs[nind_cur][2]
                        visited_tmp[xp,yp,zp] += 1
                        neigh_cat = skel_conn[xp,yp,zp]
                        # Add edge point to current (putative) edge
                        cur_edge_points.append([xp,yp,zp])
                        cur_radius.append(distmap[xp,yp,zp])
                        if neigh_cat==1 or neigh_cat>2: # terminal node or branch point
                            end_node_index = np.where((nodes[:,0]==xp) & (nodes[:,1]==yp) & (nodes[:,2]==zp))
                            if len(end_node_index[0])==0:
                                import pdb
                                pdb.set_trace()
                            conns.append([start_node_index,end_node_index[0][0]])
                            endloop = True
                            edge_points.extend(cur_edge_points)
                            radius.extend(cur_radius)
                            nedge_points.append(len(cur_edge_points))
                            for ep in cur_edge_points:
                                visited[ep[0],ep[1],ep[2]] = visited_tmp[ep[0],ep[1],ep[2]]
                                
                        # Set next origin point to current point
                        xc,yc,zc = xp,yp,zp
                    else:
                        # Dead end...
                        endloop = True
                        
                    if endloop:
                        break            
                
    n_edge_end = len(edge_points)
    #if True: #n_edge_start==n_edge_end:
    #    break
    #else:
    #    loop_count += 1
    #    if loop_count>10:
    #        import pdb
    #        pdb.set_trace()
    #        break
            
import pdb
pdb.set_trace() 

conns = arr(conns)
nedge_points = arr(nedge_points)
edge_points = arr(edge_points).astype('float')
nodes = arr(nodes).astype('float')
edge_points[:,0] *= resolution[0]
edge_points[:,1] *= resolution[1]
edge_points[:,2] *= resolution[2]
nodes[:,0] *= resolution[0]
nodes[:,1] *= resolution[1]
nodes[:,2] *= resolution[2]
radius = arr(radius)

graph = spatialgraph.SpatialGraph(initialise=True,scalars=['Radius'])      
graph.set_definition_size('VERTEX',nodes.shape[0])
graph.set_definition_size('EDGE',conns.shape[0])
graph.set_definition_size('POINT',edge_points.shape[0])
graph.set_data(nodes,name='VertexCoordinates')
graph.set_data(conns,name='EdgeConnectivity')
graph.set_data(nedge_points,name='NumEdgePoints')
graph.set_data(edge_points,name='EdgePointCoordinates')
graph.set_data(radius,name='Radius')  

ofile = join(opath,'graph_output_nolee.am')
graph.write(ofile)     
print('Saved to {}'.format(ofile))
