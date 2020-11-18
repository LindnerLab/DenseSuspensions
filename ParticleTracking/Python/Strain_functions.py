# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:08:35 2020

@author: Lars Kool
"""

import os
import numpy as np
import pandas as pd
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

def calc_strain(tracked, Cauchy=False, Hencky=False, direction='x', x_max=np.nan):
    nParticles = len(tracked.particle.unique())
    nFrames = len(tracked.frame.unique())
    
    # Set the expected sorting
    if all(tracked[:nParticles].particle == tracked.particle.unique()):
        tracked = tracked.sort_values(by=['frame','particle'],
                                      ignore_index=True)
    elif not (all(tracked[:nParticles].frame == 0) and
              all(tracked.frame[list(range(0,nFrames*nParticles,nParticles))]
                  == tracked.frame.unique())):
        tracked = tracked.sort_values(by=['frame','particle'],
                                      ignore_index=True)
    if x_max == np.nan:
        x_max = [np.max(tracked[direction][tracked.particle == i])
                 for i in range(nParticles)]
    max_length = np.tile(x_max, nFrames)
    current_length = np.array(tracked[direction])
    
    if Cauchy:
        Cauchy_strain = -(current_length - max_length)/max_length
        tracked['Cauchy_strain'] = Cauchy_strain
    if Hencky:
        Hencky_strain = -np.log(current_length/max_length)
        tracked['Hencky_strain'] = Hencky_strain
    return tracked

def calc_strain_rate(tracked, Cauchy=False, Hencky=False, dt=1, framerate=1, direction='x'): 
    if (Cauchy and 'Cauchy_strain' not in tracked.columns) or (Hencky and 'Hencky_strain' not in tracked.columns):
        print('''Required strain not found. Therefore, the strain is \n calculated first. The calculation will take a bit longer.''')
        tracked = calc_strain(tracked, Cauchy=Cauchy, Hencky=Hencky)
    
    nParticles = len(tracked.particle.unique())
    if Cauchy:
        Cauchy_strain_min_dt = np.array(tracked.Cauchy_strain[:-2*dt*nParticles])
        Cauchy_strain_plus_dt = np.array(tracked.Cauchy_strain[2*dt*nParticles:])
        Cauchy_strain_rate = ((Cauchy_strain_plus_dt - Cauchy_strain_min_dt)*framerate)/(2*dt)
        Cauchy_strain_rate = np.pad(Cauchy_strain_rate, dt*nParticles, mode='constant', constant_values=np.nan)
        tracked['Cauchy_strain_rate'] = Cauchy_strain_rate
        
    if Hencky:
        Hencky_strain_min_dt = np.array(tracked.Hencky_strain[:-2*dt*nParticles])
        Hencky_strain_plus_dt = np.array(tracked.Hencky_strain[2*dt*nParticles:])
        Hencky_strain_rate = ((Hencky_strain_plus_dt - Hencky_strain_min_dt)*framerate)/(2*dt)
        Hencky_strain_rate = np.pad(Hencky_strain_rate, dt*nParticles, mode='constant', constant_values=np.nan)
        tracked['Hencky_strain_rate'] = Hencky_strain_rate
    return tracked

if __name__ == '__main__':
    tracked = pd.read_pickle(r'E:\Lars\Oscillatory Compression\20201103 Creep\Decompression_125Start_26End_25Back\Preprocessed\V1\tracked.pkl')
    nParticles = len(tracked.particle.unique())
    x_max = [np.max(tracked['y'][tracked.particle == i] 
                    )
                 for i in range(nParticles)] 
    tracked = calc_strain(tracked, Cauchy=True, Hencky=True, direction='y', x_max=x_max)
    tracked = calc_strain_rate(tracked, Cauchy=True, Hencky=True, dt=1, framerate=10, direction='y')
    tracked.to_pickle(r'E:\Lars\Oscillatory Compression\20201103 Creep\Decompression_125Start_26End_25Back\Preprocessed\V1\tracked2.pkl')