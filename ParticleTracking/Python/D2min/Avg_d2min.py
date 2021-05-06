# -*- coding: utf-8 -*-
"""
Created on Tue May  4 10:49:08 2021

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
    
    files = next(os.walk(os.path.join(directory, 'Processed', version, 'D2min', str(interval), 'Figures')))[2]
    nFill = int(np.ceil(np.log10(files)))                                  
    idx = [i for i in range(files)
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
    _temp = f.visualize_D2min(i, directory, version, interval)
    return _temp

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

directory =  r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50'
version = 'V1'
interval = 5
nFrames = 1200

client = ipp.Client()
workers = client[:]

path = os.path.join(directory, 'Processed', version, 'D2min', str(interval))
files = next(os.walk(path))[2]

D2min_avg = np.zeros((nFrames - interval,))
for i, file in enumerate(files):
    tracked = pd.read_pickle(os.path.join(path, file))
    D2min_avg[i] = np.mean(tracked.D2min)
    
plt.figure(dpi=500)
plt.plot(np.arange(0, (nFrames-interval)/120, 1/120), D2min_avg,
         linestyle='None',
         marker='o',
         markersize=1)
plt.yscale('log')
plt.ylim([0.001, 100])
plt.xlabel('Cycle')
plt.ylabel('$\\leftangle D^2_{min} \\rightangle$')