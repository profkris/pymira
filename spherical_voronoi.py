# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:41:45 2017

@author: simon
"""

import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from scipy.spatial import SphericalVoronoi

points = np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0],
                    [0, 1, 0], [0, -1, 0], [-1, 0, 0], ])
center = np.array([0, 0, 0])
radius = 1
# calculate spherical Voronoi diagram
sv = SphericalVoronoi(points, radius, center)
# sort vertices (optional, helpful for plotting)
sv.sort_vertices_of_regions()
# generate plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# plot the unit sphere for reference (optional)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='y', alpha=0.1)