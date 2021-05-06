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
import json
import seaborn as sns
import pickle5 as pickle
import numpy as np
import pandas as pd
import ipyparallel as ipp
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt
from D2min.D2min_functions import visualize_D2min
        
def d2min_parallel(settings):
    nFrames = settings['nFrames']
    path = settings['path']
    interval = settings['interval']
    
    files = next(os.walk(os.path.join(path, 'Processed', 'V1', 'D2min', str(interval))))[2]
    nFill = int(np.ceil(np.log10(nFrames)))                                  
    idx = [i for i in range(nFrames - interval)
            if str(i).zfill(nFill) + '.pkl' not in files]
    # constants to push to the cores:
    workers.push({'settings': settings})

    # execute:
    task = workers.map(d2min_dummy, idx)
    task.wait_interactive()
    _temp = task.get()
    return

def d2min_dummy(i):
    import D2min.D2min_functions as f
    _temp = f.d2min_serial(i, settings)
    return _temp

def visualize_particles(frames, directory, version):
    for idx in frames:
        x_frame = []
        y_frame = []
        r_frame = []
        c_frame = []
        
        types = next(os.walk(os.path.join(directory, 'Preprocessed', version, 'PKL')))[1]
        N = len(types)
        palette = sns.color_palette('tab10', N)
        
        for i, _type in enumerate(types):
            path = os.path.join(directory, 'Preprocessed', version, 'PKL', _type)
            files = next(os.walk(path))[2]
            with open(os.path.join(path, files[idx]), "rb") as fh:
                P = pickle.load(fh)
            
            x_frame.append(P.x)
            y_frame.append(P.y)
            r_frame.append(P.r)
            c_frame.append(pd.Series([palette[i] for x in range(len(P.x))]))
            
        x_frame = pd.concat(x_frame)
        y_frame = pd.concat(y_frame)
        r_frame = pd.concat(r_frame)
        c_frame = pd.concat(c_frame)
        
        with open(os.path.join(directory, 'Preprocessed', version, 'settings.json')) as settings:
            crop = json.load(settings)['crop']
        fig, ax = plt.subplots(1)
        patches = []
        path = os.path.join(directory, 'RawData', version)
        files = next(os.walk(path))[2]
        img = plt.imread(os.path.join(path, files[idx]))
        img = img[crop[0]:crop[1]][crop[2]:crop[3]]
        ax.imshow(img,cmap='gray')
        for x, y, r, c in zip(x_frame, y_frame, r_frame, c_frame):
            circle = mpl.patches.Circle((y, x), r,
                                        facecolor=c,
                                        alpha=0.7)
            patches.append(circle)
        p = mpl.collections.PatchCollection(patches, match_original=True)
        ax.add_collection(p)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.show()
    
path = r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50'
interval = 120
save_files = True
nFrames = 1200
neighbor_cutoff = 70

client = ipp.Client()
workers = client[:]

settings = {'path': path,
            'interval': interval,
            'nFrames': nFrames,
            'neighbor_cutoff': neighbor_cutoff,
            'save_files': save_files}
    
d2min_parallel(settings)

# visualize_D2min([3], r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50', 'V1', 120, vmax=0.2, save_files=True)
# visualize_particles(range(241, 360), r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50', r'V1')