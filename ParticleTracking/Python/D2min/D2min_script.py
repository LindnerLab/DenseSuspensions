# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 14:53:33 2020

@author: Lars Kool
"""

import os
import numpy as np
import pandas as pd
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt # For convenience
from D2min_functions import neighboring_particles, d2min_particle, deviatoric_strain

start_frame = 56
interval = 120
neighbor_cutoff = 70
avg_distance = 45
vmin = 0
vmax = 0.1
verbose = False
save_files = False
path = r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25'
output = 'Pmax'

tracked_complete = pd.read_pickle(os.path.join(path, r'Preprocessed\V1\Complete_tracks_renumbered.pkl'))

nParticles = len(tracked_complete.particle.unique())
nFrames = len(tracked_complete.frame.unique())
nCycles = m.floor((nFrames-start_frame)/interval)
d2min = np.zeros([nParticles, nCycles])
M = np.zeros([nParticles, nCycles, 3, 3])
strain_deviatoric = np.zeros([nParticles, nCycles])

for i in range(nCycles):
    frame = start_frame + i*interval
    xy_frame = np.array([tracked_complete[tracked_complete.frame == frame].x,
                         tracked_complete[tracked_complete.frame == frame].y]
                        ).transpose()
    xy_next = np.array([tracked_complete[tracked_complete.frame == frame+interval].x,
                         tracked_complete[tracked_complete.frame == frame+interval].y]
                        ).transpose()
    r_frame = np.array([tracked_complete[tracked_complete.frame == frame].r])
    neighbors = neighboring_particles(xy_frame, neighbor_cutoff)
    
    for j in range(nParticles):
        xy_particle = xy_frame[j, :].reshape([1, 2])
        xy_particle_next = xy_next[j, :].reshape([1, 2])
        xy_neighbors = xy_frame[neighbors[j], :]
        xy_neighbors_next = xy_next[neighbors[j], :]
        d2min[j, i], M[j, i, :, :] = d2min_particle(xy_particle, xy_particle_next, xy_neighbors, xy_neighbors_next)
        strain_deviatoric[j, i] = deviatoric_strain(M[j, i, :, :])
        
    if verbose:
        fig, ax = plt.subplots(1)
        patches = []
        img = plt.imread(os.path.join(r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25\RawData\V1', '_'+str(frame).zfill(5)+'.tif'))
        crop = [0, 2303, 158, 2178]
        img = img[crop[0]:crop[1]][crop[2]:crop[3]]
        ax.imshow(img,cmap='gray')
        norm = mpl.colors.Normalize(vmin=0, vmax=vmax, clip=True)
        cmap = mpl.cm.get_cmap('viridis')
        for x, y, r, c in zip(xy_frame[:, 0], xy_frame[:, 1], r_frame[:].transpose(), cmap(norm(d2min[:, i]))):
            circle = mpl.patches.Circle((x,y), r, facecolor=c)
            patches.append(circle)
        p = mpl.collections.PatchCollection(patches, match_original=True)
        ax.add_collection(p)
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
        cbar.set_label('$D^{2}_{min}/d_{avg}^{2} (-)$')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.show()

if save_files:
    np.savetxt(os.path.join(path, r'Processed\D2min', output, 'D2min.csv'),
               d2min, delimiter = ',')
    
bins = 10**np.arange(-2,3,0.1)
d2min_binned = np.zeros((len(bins)-1,3))
for i in range(len(bins)-1):
    idx = np.logical_and(strain_deviatoric >= bins[i],
                                                     strain_deviatoric < bins[i+1])
    d2min_binned[i,2] = np.sum(idx)
    d2min_binned[i,0] = 10**np.mean(np.log10(d2min[idx]))
    d2min_binned[i,1] = np.std(np.log10(d2min[idx]))/np.sqrt(d2min_binned[i,2])
    

plt.figure(dpi=500)
plt.plot(bins[:-1],d2min_binned[:,0],
             linestyle='none',
             marker='o',
             markersize=5,
             markerfacecolor='none',
             markeredgecolor='red')
plt.scatter(strain_deviatoric, d2min,
          marker='o',
          s=5,
          facecolor='gray',
          edgecolors='none',
          alpha=0.2)
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.01, 1000])
plt.xlabel('Deviatoric strain (-)')
plt.ylabel('$D^{2}_{min}/D^{2} (-)$')
plt.show()
