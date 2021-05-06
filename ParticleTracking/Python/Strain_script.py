# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:08:35 2020

@author: Lars Kool
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def calc_strain(tracked, Cauchy=False, Hencky=False, direction='y', no_overwrite=False):
    nParticles = len(tracked.particle.unique())
    nFrames = len(tracked.frame.unique())

    # # Set the expected sorting
    # if all(tracked[:nParticles].particle == tracked.particle.unique()):
    #     tracked = tracked.sort_values(by=['frame', 'particle'],
    #                                   ignore_index=True)
    # elif not (all(tracked[:nParticles].frame == 0) and
    #           all(tracked.frame[list(range(0,nFrames*nParticles,nParticles))]
    #               == tracked.frame.unique())):
    #     tracked = tracked.sort_values(by=['frame', 'particle'],
    #                                   ignore_index=True)

    # x_max = np.array([np.abs(np.max(tracked[tracked.particle == i][direction]))
                      # for i in range(nParticles)])
    f_max = np.argmin([np.mean(tracked[tracked.frame == i]) for i in range(50)])
    x_max = np.array(tracked[tracked.frame == f_max][direction])
    max_length = np.tile(x_max, nFrames)
    current_length = np.abs(np.array(tracked[direction]))

    if Cauchy and not ('Cauchy_strain' in tracked.columns and no_overwrite):
        Cauchy_strain = -(current_length - max_length) / max_length
        tracked.loc[:, 'Cauchy_strain'] = Cauchy_strain
    if Hencky and not ('Hencky_strain' in tracked.columns and no_overwrite):
        Hencky_strain = -np.log(current_length / max_length)
        tracked.loc[:, 'Hencky_strain'] = Hencky_strain
    return tracked


def calc_strain_rate(tracked, Cauchy=False, Hencky=False, direction='y', dt=1, framerate=1, no_overwrite=False):
    if (Cauchy and 'Cauchy_strain' not in tracked.columns) or (Hencky and 'Hencky_strain' not in tracked.columns):
        print('''Required strain not found. Therefore, the strain is calculated first. The calculation will take a bit longer.''')
        tracked = calc_strain(tracked, Cauchy=Cauchy, Hencky=Hencky, direction='y', no_overwrite=no_overwrite)

    nParticles = len(tracked.particle.unique())
    if Cauchy and not ('Cauchy_strain_rate' in tracked.columns and no_overwrite):
        Cauchy_strain_min_dt = np.array(tracked.Cauchy_strain[:-2 * dt * nParticles])
        Cauchy_strain_plus_dt = np.array(tracked.Cauchy_strain[2 * dt * nParticles:])
        Cauchy_strain_rate = ((Cauchy_strain_plus_dt - Cauchy_strain_min_dt) * framerate) / (2 * dt)
        Cauchy_strain_rate = np.pad(Cauchy_strain_rate, dt * nParticles, mode='constant', constant_values=np.nan)
        tracked['Cauchy_strain_rate'] = Cauchy_strain_rate

    if Hencky and not ('Hencky_strain_rate' in tracked.columns and no_overwrite):
        Hencky_strain_min_dt = np.array(tracked.Hencky_strain[:-2 * dt * nParticles])
        Hencky_strain_plus_dt = np.array(tracked.Hencky_strain[2 * dt * nParticles:])
        Hencky_strain_rate = ((Hencky_strain_plus_dt - Hencky_strain_min_dt) * framerate) / (2 * dt)
        Hencky_strain_rate = np.pad(Hencky_strain_rate, dt * nParticles, mode='constant', constant_values=np.nan)
        tracked['Hencky_strain_rate'] = Hencky_strain_rate
    return tracked


if __name__ == '__main__':
    path = r'G:\Lars\Oscillatory Compression\20210421\Avg125_Amp100_Back25_Per600_C40'
    file = r'Preprocessed\V1\tracked.pkl'

    tracked = pd.read_pickle(os.path.join(path, file))
    nParticles = len(tracked.particle.unique())
    nFrames = len(tracked.frame.unique())

    tracked = calc_strain(tracked, Cauchy=True, Hencky=True, direction='y')
    tracked = calc_strain_rate(tracked, Cauchy=True, Hencky=True, dt=2, framerate=0.2, no_overwrite=False)
    tracked.to_pickle(os.path.join(path, file))

    strain_avg = [np.mean(tracked[tracked.frame == i].Cauchy_strain) for i in range(nFrames)]
    
    i = 0
    plt.figure(dpi=500)
    plt.plot(-100*np.cos(2*np.pi*np.arange(1079 - i, nFrames - i)/(600*0.2))+100, strain_avg[1079:],
             linestyle = '-',
             linewidth = 1,
             marker = 'o',
             markersize = 2)
    plt.xlabel('P (mBar)')
    plt.ylabel('$\\leftangle \epsilon_{Hencky} \\rightangle$ (-)')