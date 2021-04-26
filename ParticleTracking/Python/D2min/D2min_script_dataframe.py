# -*- coding: utf-8 -*-
"""
Contributors : Lars Kool
Affiliations : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France

This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162
More info on this programme: https://caliper-itn.org/
"""

import os
import pickle5 as pickle
import numpy as np
import pandas as pd
import ipyparallel as ipp
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from D2min_functions import neighboring_particles, d2min_particle, deviatoric_strain

def d2min_serial(tracked_complete, i, interval, neighbor_cutoff, nParticles, path, nFill, verbose, save_files):
    d2min = np.zeros([nParticles])
    M = [[] for i in range(nParticles)]
    strain_deviatoric = np.zeros([nParticles])
    
    x = tracked_complete[tracked_complete.frame == i].x
    y = tracked_complete[tracked_complete.frame == i].y
    frames = tracked_complete[tracked_complete.frame == i].frame
    xy_frame = np.array([x, y]).transpose()
    xy_next = np.array([tracked_complete[tracked_complete.frame == i+interval].x,
                         tracked_complete[tracked_complete.frame == i+interval].y]
                        ).transpose()
    r_frame = np.array([tracked_complete[tracked_complete.frame == i].r])
    neighbors = neighboring_particles(xy_frame, neighbor_cutoff)
    
    for j in range(nParticles):
        xy_particle = xy_frame[j, :].reshape([1, 2])
        xy_particle_next = xy_next[j, :].reshape([1, 2])
        xy_neighbors = xy_frame[neighbors[j], :]
        xy_neighbors_next = xy_next[neighbors[j], :]
        d2min[j], M[j] = d2min_particle(xy_particle, xy_particle_next, xy_neighbors, xy_neighbors_next)
    strain_deviatoric = deviatoric_strain(M)
    _dict = {'x': x,
             'y': y,
             'r': r_frame.reshape((nParticles, )),
             'frame': frames,
             'Affine': M,
             'D2min': d2min,
             'Deviatoric': strain_deviatoric}
    D2min = pd.DataFrame(_dict)
    if save_files:
        D2min.to_pickle(os.path.join(path, 'Processed', 'D2min', str(interval),
                                     str(i).zfill(nFill) + '.pkl'))
        
def d2min_dummy(i):
    import D2min_script_dataframe
    d2min_serial(tracked_complete, i, interval, neighbor_cutoff, nParticles, path, nFill, verbose, save_files)
    return

    # if verbose:
    #     fig, ax = plt.subplots(1)
    #     patches = []
    #     img = plt.imread(os.path.join(r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25\RawData\V1', '_'+str(i).zfill(5)+'.tif'))
    #     crop = [0, 2303, 158, 2178]
    #     img = img[crop[0]:crop[1]][crop[2]:crop[3]]
    #     ax.imshow(img,cmap='gray')
    #     norm = mpl.colors.Normalize(vmin=0, vmax=vmax, clip=True)
    #     cmap = mpl.cm.get_cmap('viridis')
    #     for x, y, r, c in zip(xy_frame[:, 0], xy_frame[:, 1], r_frame[:].transpose(), cmap(norm(d2min[:, i]))):
    #         circle = mpl.patches.Circle((x,y), r, facecolor=c)
    #         patches.append(circle)
    #     p = mpl.collections.PatchCollection(patches, match_original=True)
    #     ax.add_collection(p)
    #     cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    #     cbar.set_label('$D^{2}_{min}/d_{avg}^{2} (-)$')
    #     ax.get_xaxis().set_visible(False)
    #     ax.get_yaxis().set_visible(False)
    #     plt.show()

def d2min_parallel(tracked_complete, interval, neighbor_cutoff, nParticles, path, nFill, verbose, save_files):
    files = next(os.walk(os.path.join(path, 'Processed', 'D2min', str(interval))))[2]
    idx = [i for i in range(nFrames - interval) if str(i).zfill(nFill) + '.pkl' not in files]
    # constants to push to the cores:
    workers.push({'tracked_complete': tracked_complete,
                  'interval': interval,
                  'nParticles': nParticles,
                  'neighbor_cutoff': neighbor_cutoff,
                  'path': path,
                  'nFill': nFill,
                  'verbose': verbose,
                  'save_files': save_files})
    # execute:
    _task = workers.map(d2min_dummy, idx)
    _task.wait_interactive()
    _temp = _task.get()
    return

start_frame = 56
interval = 5
neighbor_cutoff = 70
avg_distance = 45
vmin = 0
vmax = 0.1
verbose = False
save_files = True
nFill = 5
path = r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50'
pickle.HIGHEST_PROTOCOL = 5

with open(os.path.join(path, r'Preprocessed\V1\tracked.pkl'), "rb") as fh:
  tracked_complete = pickle.load(fh)

client = ipp.Client()
workers = client[:]

nParticles = len(tracked_complete.particle.unique())
nFrames = len(tracked_complete.frame.unique())
nCycles = m.floor((nFrames-start_frame)/interval)

d2min_parallel(tracked_complete, interval, neighbor_cutoff, nParticles, path, nFill, verbose, save_files)
