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
        
def d2min_parallel(settings):
    nFrames = settings['nFrames']
    path = settings['path']
    interval = settings['interval']
    
    files = next(os.walk(os.path.join(path, 'Processed', 'D2min', str(interval))))[2]    
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

path = r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50'
interval = 5
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

