# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 07:20:06 2017

@author: simon
"""

from pymira import spatialgraph
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

file = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\spatialGraph_RIN.am'
cfile = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\crawl_recon\crawl_recon.am'
efile = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1\ca1_kt0p0001\exposure\exposure60.0.am'

file = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\SW1222\1\spatialGraph_RIN.am'
cfile = None
efile = None

graph = spatialgraph.SpatialGraph()
print 'Reading graph...'
graph.read(file)

if cfile is not None:
    cgraph = spatialgraph.SpatialGraph()
    print 'Reading graph...'
    cgraph.read(cfile)

if efile is not None:
    egraph = spatialgraph.SpatialGraph()
    print 'Reading graph...'
    egraph.read(efile)

coords = graph.get_data('VertexCoordinates')
pts = graph.get_data('EdgePointCoordinates')
flow = graph.get_data('Flow')
radii = graph.get_data('Radii')
distance = cgraph.get_data('Distance')
delay = cgraph.get_data('Delay')
exposure = egraph.get_data('Exposure')
#coords = np.asarray([[inds[0][i],inds[1][i],inds[2][i]] for i,tmp in enumerate(inds[0])])

from scipy.spatial import ConvexHull
hull = ConvexHull(coords)
edges = zip(*coords)
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

hullcoords = coords[hull.vertices,:]
ax.scatter(hullcoords[:,0],hullcoords[:,1],hullcoords[:,2])

sc = [1.,0.75,0.5,0.25,0.1,0.]
nsc = len(sc-1)

#com = np.mean(coords,axis=0)
#dist = np.asarray([np.linalg.norm(c-com) for c in coords])

dist_from_hull = np.abs(np.max(np.dot(hull.equations[:, :-1], coords.T).T + hull.equations[:, -1], axis=-1))
pts_dist_from_hull = np.abs(np.max(np.dot(hull.equations[:, :-1], pts.T).T + hull.equations[:, -1], axis=-1))

path2 = r'C:\Users\simon\Dropbox\160113_paul_simulation_results\LS147T\1'
mfile = os.path.join(path2,'match_exposure_flow_radii.dill')
import dill as pickle
with open(mfile,'rb') as fo:
    match = pickle.load(fo)
    
match,params = match[0],match[1]
distance = match[:,5]
radii = match[:,2]

fig = plt.figure()
plt.scatter(pts_dist_from_hull,flow)
plt.ylim(-200,2000)
plt.xlim(-200,2500)
plt.xlabel('Distance from perimeter ($\mu$m)')
plt.ylabel('Intravascular flow (nL /min)')

fig = plt.figure()
plt.scatter(pts_dist_from_hull,radii)
#plt.ylim(-200,2000)
#plt.xlim(-200,2500)
plt.xlabel('Distance from perimeter ($\mu$m)')
plt.ylabel('Vessel radius ($\mu$m)')

fig = plt.figure()
plt.scatter(pts_dist_from_hull[exposure>0],exposure[exposure>0])
#plt.ylim(-200,2000)
plt.ylim(0, 0.002)
plt.xlabel('Distance from perimeter ($\mu$m)')
plt.ylabel('Exposure')

flagInd = 3
flag = match[:,flagInd]
fig = plt.figure()
inds0 = np.where((flag==0) & (exposure>0.))[0]
inds1 = np.where((flag==1) & (exposure>0.))[0]
plt.scatter(pts_dist_from_hull[inds0],exposure[inds0],color='green',marker='.')
plt.scatter(pts_dist_from_hull[inds1],exposure[inds1],color='blue',marker='.')
ax = plt.gca()
#ax.set_yscale('log')
plt.ylim(1e-5,1e-2)
plt.xlim(-200,2500)
plt.xlabel('Intravascular distance from inlet ($\mu$m)')
plt.ylabel('Exposure (mM min /$um^{2}$)')

fig = plt.figure()
plt.scatter(distance,delay)
plt.ylim(-200,2000)
plt.xlim(-1000,10000)
plt.xlabel('Intravascular distance from inlet ($\mu$m)')
plt.ylabel('Intravascular flow (nL /min)')

fig = plt.figure()
#inds0 = np.where((flag==0) & (pts_dist_from_hull>=0.))[0]
#inds1 = np.where((flag==1) & (pts_dist_from_hull>=0.))[0]
inds0 = np.where((flag[0::10]==1) & (pts_dist_from_hull>=0.))[0]
inds1 = np.where((flag[0::10]==2) & (pts_dist_from_hull>=0.))[0]
figsize = (4,8)
plt.figure(figsize=figsize)
d1 = pts_dist_from_hull[inds0]
d2 = pts_dist_from_hull[inds1]
bplot = plt.boxplot([d1,d2],widths=(0.6,0.6),notch=False,vert=True,sym='ko',patch_artist=True,showfliers=False)
plt.ylabel('Distance from perimeter $\mu$m')
plt.xticks([1,2],['Non-labelled','Labelled'])
colors = ['lightgreen','lightblue']
for patch, color in zip(bplot['boxes'], colors): 
    patch.set_facecolor(color)
#plt.scatter(distance[inds0],pts_dist_from_hull[inds0],color='green',marker='.')
#plt.scatter(distance[inds1],pts_dist_from_hull[inds1],color='blue',marker='.')

fig = plt.figure()
inds0 = np.where((flag==0) & (distance>=0.))[0]
inds1 = np.where((flag==1) & (distance>=0.))[0]
plt.scatter(distance[inds0],radii[inds0],color='green',marker='.')
plt.scatter(distance[inds1],radii[inds1],color='blue',marker='.')
#plt.ylim(-200,2000)
#plt.xlim(-200,2500)
plt.xlabel('Intravascular distance from inlet ($\mu$m)')
plt.ylabel('Vessel radius ($\mu$m)')