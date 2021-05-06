# -*- coding: utf-8 -*-
"""
Created on Thu May  6 12:05:06 2021

@author: Charly
"""

import os
import numpy as np
import pandas as pd
import ipyparallel as ipp
import matplotlib as mpl
import matplotlib.pyplot as plt

def visualize_d2min_parallel(settings):
    workers.push({'settings': settings})
    nFrames = settings['nFrames']
    files = next(os.walk(os.path.join(directory, 'Processed', version, 'D2min', str(interval), 'Figures')))[2]
    nFill = int(np.ceil(np.log10(nFrames)))                                  
    idx = [i for i in range(nFrames)
           if str(i).zfill(nFill) + '.jpg' not in files]
    
    task = workers.map(visualize_d2min_dummy, idx)
    task.wait_interactive()
    _temp = task.get()
    return

def visualize_d2min_dummy(i):
    import D2min.D2min_functions as f
    directory = settings['directory']
    version = settings['version']
    interval = settings['interval']
    save_files = settings['save_files']
    _temp = f.visualize_D2min(i, directory, version, interval, vmax=0.2, save_files=save_files)
    return _temp

directory =  r'G:\Lars\Oscillatory Compression\20210331\Avg125_Amp100_Back25_Per600_C60'
version = 'V1'
interval = 5
nFrames = 1200

client = ipp.Client()
workers = client[:]

settings = {'directory': directory,
            'version': version,
            'interval': interval,
            'nFrames': nFrames,
            'save_files': True}

visualize_d2min_parallel(settings)