# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 14:25:01 2020

@author: Lars Kool
"""
import numpy as np
import pandas as pd
import scipy.signal as scisig
import matplotlib as mpl
import matplotlib.pyplot as plt
import Pair_Correlation as pc
import scipy

if __name__ == '__main__':
    tracked_complete = pd.read_pickle('F:/Lars/Oscillatory Compression/20200713 Soft Particle/Avg100_Amp80_Per120_Back25/Preprocessed/Complete_tracks_large.pkl')
    nFiles = 1384
    per = 120
    nCycles = np.floor(nFiles/per).astype('int')
    
    track_i = tracked_complete[tracked_complete.particle == 10]
    
    loc_max = scisig.find_peaks(track_i.x,distance=110)[0]
    loc_min = scisig.find_peaks(-track_i.x, distance=110)[0]
    
    D = np.zeros([nCycles-1,])
    
    # for i in range(nCycles-1):
    #     cycle1 = np.arange(loc_min[0]+per*(i+0),loc_min[0]+per*(i+1))
    #     track1 = np.array([tracked_complete[np.isin(tracked_complete.frame,cycle1) & np.array(tracked_complete.particle == 0)].x,tracked_complete[np.isin(tracked_complete.frame,cycle1) & np.array(tracked_complete.particle == 0)].y])
        
    #     cycle2 = np.arange(loc_min[0]+per*(i+1),loc_min[0]+per*(i+2))
    #     track2 = np.array([tracked_complete[np.isin(tracked_complete.frame,cycle2) & np.array(tracked_complete.particle == 0)].x,tracked_complete[np.isin(tracked_complete.frame,cycle2) & np.array(tracked_complete.particle == 0)].y])
        
    #     D[i] = np.sum(np.sqrt(np.sum((track2-track1)**2)))
    
    plt.figure(dpi=500)    
    for i in range(0,2):
        idx = np.arange(loc_min[0]+per*(i+0),loc_min[0]+per*(i+1))
        x = track_i.x[idx]
        y = track_i.y[idx]
        
        plt.plot(x,y)
    plt.xlabel('x-pos (px)')
    plt.ylabel('y-pos (px)')
        
    t = pd.read_pickle('F:/Lars/Oscillatory Compression/20200713 Soft Particle/Avg100_Amp80_Per120_Back25/Preprocessed/Complete_tracks_large.pkl')
    t = t.append(pd.read_pickle('F:/Lars/Oscillatory Compression/20200713 Soft Particle/Avg100_Amp80_Per120_Back25/Preprocessed/Complete_tracks_small.pkl'))
    
    plt.figure(dpi=500)
    for i in [0,10]:
        x_pc = np.array(t.x[t.frame == loc_max[i]])
        y_pc = np.array(t.y[t.frame == loc_max[i]])
        
        [g_r, r, ref_idx] = pc.pairCorrelationFunction_2D(x_pc,y_pc,2000,200,5)
        plt.plot(r/(20.5*2), g_r)
    plt.xlabel(r'$r/D_{small}$ (-)')
    plt.ylabel('g(r) (-)')
    plt.legend(['nCycles = 1', 'nCycles = 8'])
    
    plt.figure(dpi=500)
    x_pc = np.array(t.x[t.frame == 107])
    y_pc = np.array(t.y[t.frame == 107])
    
    [g_r, r, ref_idx] = pc.pairCorrelationFunction_2D(x_pc,y_pc,2000,200,5)
    plt.plot(r/(20.5*2),g_r)
    
    x_pc = np.array(t.x[t.frame == 161])
    y_pc = np.array(t.y[t.frame == 161])
    
    [g_r, r, ref_idx] = pc.pairCorrelationFunction_2D(x_pc,y_pc,2000,200,5)
    plt.plot(r/(20.5*2), g_r)
    
    plt.xlabel(r'$r/D_{small}$ (-)')
    plt.ylabel('g(r) (-)')
    plt.legend([r'$\phi_{min}$', r'$\phi_{max}$'])
    
    
    
    Img = plt.imread('F:/Lars/Oscillatory Compression/20200713 Soft Particle/Avg100_Amp80_Per120_Back25/RawData/_00107.tif')
    Img = Img[85:85+2047,25:25+2279]
    
    plt.figure(dpi=500)
    plt.imshow(Img,cmap='gray')
    
    for i in np.unique(tracked_complete.particle):
        idx = np.array(tracked_complete.particle == i)
        plt.plot(tracked_complete.x[idx], tracked_complete.y[idx])
    
    
    # for i in np.unique(tracked_complete.particle).astype('int'):
    #     track_i = tracked_complete[tracked_complete.particle == i]
        
    #     loc_max = scisig.find_peaks(track_i.x,distance=110)[0]
    #     loc_min = scisig.find_peaks(-track_i.x, distance=110)[0]
        
    #     C1 = np.array([track_i.x[loc_min[0]+per*(i+0):loc_min[0]+per*(i+1)],track_i.y[loc_min[0]+per*(i+0):loc_min[0]+per*(i+1)]])
    #     C2 = np.array([track_i.x[loc_min[0]+per*(i+1):loc_min[0]+per*(i+2)],track_i.y[loc_min[0]+per*(i+1):loc_min[0]+per*(i+2)]])
        
    #     D = np.sum(np.sqrt(np.sum((C2-C1)**2)))