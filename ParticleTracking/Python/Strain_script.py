# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:08:35 2020

@author: Lars Kool
"""

import os
import numpy as np
import pandas as pd

def calc_strain(tracked, Cauchy=False, Hencky=False, direction='y', no_overwrite=False):
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
    
    x_max = [np.max(tracked[tracked.particle == i][direction])
             for i in range(nParticles)]
    max_length = np.tile(x_max, nFrames)
    current_length = np.array(tracked[direction])
    
    if Cauchy and not ('Cauchy_strain' in tracked.columns and no_overwrite):
        Cauchy_strain = -(current_length - max_length)/max_length
        tracked['Cauchy_strain'] = Cauchy_strain
    if Hencky and not ('Hencky_strain' in tracked.columns and no_overwrite):
        Hencky_strain = -np.log(current_length/max_length)
        tracked['Hencky_strain'] = Hencky_strain
    return tracked

def calc_strain_rate(tracked, Cauchy=False, Hencky=False, direction='y', dt=1, framerate=1, no_overwrite=False): 
    if (Cauchy and 'Cauchy_strain' not in tracked.columns) or (Hencky and 'Hencky_strain' not in tracked.columns):
        print('''Required strain not found. Therefore, the strain is calculated first. The calculation will take a bit longer.''')
        tracked = calc_strain(tracked, Cauchy=Cauchy, Hencky=Hencky, direction='y', no_overwrite=no_overwrite)
    
    nParticles = len(tracked.particle.unique())
    if Cauchy and not ('Cauchy_strain_rate' in tracked.columns and no_overwrite):
        Cauchy_strain_min_dt = np.array(tracked.Cauchy_strain[:-2*dt*nParticles])
        Cauchy_strain_plus_dt = np.array(tracked.Cauchy_strain[2*dt*nParticles:])
        Cauchy_strain_rate = ((Cauchy_strain_plus_dt - Cauchy_strain_min_dt)*framerate)/(2*dt)
        Cauchy_strain_rate = np.pad(Cauchy_strain_rate, dt*nParticles, mode='constant', constant_values=np.nan)
        tracked['Cauchy_strain_rate'] = Cauchy_strain_rate
        
    if Hencky and not ('Hencky_strain_rate' in tracked.columns and no_overwrite):
        Hencky_strain_min_dt = np.array(tracked.Hencky_strain[:-2*dt*nParticles])
        Hencky_strain_plus_dt = np.array(tracked.Hencky_strain[2*dt*nParticles:])
        Hencky_strain_rate = ((Hencky_strain_plus_dt - Hencky_strain_min_dt)*framerate)/(2*dt)
        Hencky_strain_rate = np.pad(Hencky_strain_rate, dt*nParticles, mode='constant', constant_values=np.nan)
        tracked['Hencky_strain_rate'] = Hencky_strain_rate
    return tracked

if __name__ == '__main__':
    path = r'F:\Lars\Oscillatory Compression\20201103 Creep\Compression_26Start_125End_25_Back'
    
    tracked = pd.read_pickle(os.path.join(path, r'Preprocessed\V1\tracked.pkl'))
    nParticles = len(tracked.particle.unique())
    nFrames = len(tracked.frame.unique())
    
    tracked = calc_strain(tracked, Cauchy=True, Hencky=True,direction='y')
    tracked = calc_strain_rate(tracked, Cauchy=True, Hencky=True, dt=2, framerate=0.25, no_overwrite=False)
    tracked.to_pickle(os.path.join(path, r'Preprocessed\V1\tracked.pkl'))