# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 11:24:58 2021

@author: larsk
"""

import numpy as np
import scipy
import scipy.spatial
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def Weighted_Voronoi_2D(x, y, W):
    Voronoi = scipy.spatial.Voronoi(np.array([x, y, W]).transpose(),
                                    furthest_site=True)
    Voronoi._points = Voronoi._points[:, 0:2]
    Voronoi.vertices = Voronoi.vertices[:, 0:2]
    Voronoi.ndim = 2
    Voronoi.min_bound = Voronoi.min_bound[0:2]
    Voronoi.max_bound = Voronoi.max_bound[0:2]
    return Voronoi


tracked = pd.read_pickle(r'C:\Users\larsk\OneDrive\Documenten\PhD\Projects\20210312\particles.pkl')
z = np.sqrt(10 - tracked.r)
Delaunay = scipy.spatial.Delaunay(np.array([tracked.x, tracked.y]).transpose())
Voronoi = scipy.spatial.Voronoi(np.array([tracked.x, tracked.y]).transpose())
# Voronoi = Weighted_Voronoi_2D(tracked.x, tracked.y, z)

color = [[] for i in range(2500)]
color[:1500] = ['blue' for i in range(1500)]
color[1500:] = ['red' for i in range(1000)]

fig, ax = plt.subplots(dpi=500)
plt.triplot(tracked.x, tracked.y, Delaunay.simplices, linewidth=0.2)
patches = []
for i in range(len(tracked)):
    circle = mpl.patches.Circle((tracked.x[i], tracked.y[i]), tracked.r[i],
                                facecolor=color[i],
                                alpha=0.4)
    patches.append(circle)
p = mpl.collections.PatchCollection(patches, match_original=True)
ax.add_collection(p)
ax.axis('equal')
plt.xlim([-0.1, 10.1])
plt.ylim([-0.1, 5.1])
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
plt.show()

fig, ax = plt.subplots(dpi=500)
patches = []
for i in range(len(tracked)):
    circle = mpl.patches.Circle((tracked.x[i], tracked.y[i]), tracked.r[i],
                                facecolor=color[i],
                                alpha=0.4)
    patches.append(circle)
p = mpl.collections.PatchCollection(patches, match_original=True)
ax.add_collection(p)
scipy.spatial.voronoi_plot_2d(Voronoi,
                              ax=ax,
                              line_width=0.2,
                              point_size=0.4,
                              show_vertices=False,
                              show_points=False)
# ax.axis('equal')
plt.xlim([-0.1, 10.1])
plt.ylim([-0.1, 5.1])
plt.show()