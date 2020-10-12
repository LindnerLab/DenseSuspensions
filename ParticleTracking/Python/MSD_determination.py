# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:49:39 2020

@author: Lars Kool
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    directory = r'F:\Lars\Oscillatory Compression\20200921 Instrumental Noise\000mbar\0.1fps\Preprocessed\V3'
    Plarge = pd.read_pickle(directory+r'\Complete_tracks_large.pkl')
    Psmall = pd.read_pickle(directory+r'\Complete_tracks_small.pkl')
    Psmall.particle = Psmall.particle+np.max(Plarge.particle)+1
    Pall = pd.concat([Plarge,Psmall],ignore_index=True)
    del Psmall, Plarge
    nParticles = len(Pall.particle.unique())
    nFrames = len(Pall.frame.unique())
    
    SD = np.zeros([nFrames-1,nParticles])
    
    for i, p in enumerate(Pall.particle.unique()):
        x = np.array(Pall[Pall.particle == p].x)
        y = np.array(Pall[Pall.particle == p].y)
        
        dx = x[:-1]-x[1:]
        dy = y[:-1]-y[1:]
    
        SD[:,i] = np.square(dx) + np.square(dy)
    
    MSD = [np.mean(SD,axis=0), np.std(SD,axis=0)]
    
    MSD_tot = [np.sqrt(np.mean(SD)),np.sqrt(np.std(SD))]
        
    H1 = np.histogram(SD, bins=1000)
    x = H1[1][:-1]+H1[1][1]/2
    y = H1[0]/np.max(H1[0])
    plt.figure(dpi=500)
    plt.plot(x, y,
             marker='.',
             linestyle='')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim([0, 0.03])
    plt.ylim([0, 1.1])
    plt.xlabel('Square Displacement (px^2)')
    plt.ylabel('Freq/max(Freq) (-)')

    H2 = np.histogram(np.sqrt(SD), bins=100)
    x = H2[1][:-1]+H2[1][1]/2
    y = H2[0]/np.max(H2[0])
    plt.figure(dpi=500)
    plt.plot(x, y,
             marker='.',
             linestyle='')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim([0, 0.3])
    plt.ylim([0, 1.1])
    plt.xlabel('Displacement (px)')
    plt.ylabel('Freq/max(Freq) (-)')