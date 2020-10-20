# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:03:48 2020

@author: lars Kool
"""
import numpy as np
import pandas as pd
import os
from Track_repair import trackRepair
import trackpy as tp

def DataParser(datapath, dataFormat):
    result = pd.DataFrame(columns=['y', 'x', 'r','frame'])
    nParticles = []
    
    if dataFormat == 'csv':
        files = [file for file in os.listdir(datapath) if file.lower().endswith('.csv')]
        nFiles = len(files)
        
        result = pd.concat([pd.read_csv(datapath + file, names=['x','y', 'r','frame']) for file in files], ignore_index = True)
    
    elif dataFormat == 'pkl':
        files = [file for file in os.listdir(datapath) if file.endswith('.pkl')]
        nFiles = len(files)
        
        for idx, file in enumerate(files):
            nParticles.append(len(pd.read_pkl(datapath + file)))
            result = pd.concat([result, pd.read_pkl(datapath + file)], ignore_index = True)
    else:
        print('No supported file format selected!')
    
    return result, nFiles

if __name__ == "__main__":   
    datapath = r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25\Preprocessed\V1'
    prefix1 = '0'
    dataFormat = 'csv'
    minLength = 3500
    
    # # Import and parse data
    # [imported, nFiles] = DataParser(datapath+r'\CSV\\', dataFormat)
    # print('Particles imported')
    
    # # Separate small from large particles, this makes tracking faster (as there are less permutations) and more reliable (as small particles should never by matched with large ones and viceversa)
    # Psmall = imported[imported.r == 20.5].reset_index(drop=True)
    # Plarge = imported[imported.r == 28].reset_index(drop=True)
    # del imported
    
    # # Track the particles
    # search_range = 10
    # tracked_small = tp.link_df(Psmall, search_range, memory=5).sort_values(by=['particle','frame'], ignore_index=True)
    # tracked_large = tp.link_df(Plarge, search_range, memory=5).sort_values(by=['particle','frame'], ignore_index=True)

    # # Remove small tracks and repair small holes (linear interpolation) and fill large holes with NaN
    # complete_small = trackRepair(tracked_small,10,nFiles, minLength)
    # complete_large = trackRepair(tracked_large,10,nFiles, minLength)
    
    # complete_small.to_pickle(datapath+r'\Complete_tracks_small.pkl')
    # complete_large.to_pickle(datapath+r'\Complete_tracks_large.pkl')
    
    complete_small = pd. read_pickle(os.path.join(datapath+r'\Complete_tracks_small.pkl'))
    complete_large = pd. read_pickle(os.path.join(datapath+r'\Complete_tracks_large.pkl'))
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
    Crop = [0, 2303, 158, 2178]
    Img = plt.imread(r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25\RawData\V1\_00099.tif')[Crop[0]:Crop[1]][Crop[2]:Crop[3]]
    fig,ax = plt.subplots(1)
    ax.imshow(Img)
    ax.add_patch(patches.Circle((494.097, 662.771),28))
    idx = complete_small.frame == 100
    for xx, yy in zip(complete_small[idx].x, complete_small[idx].y):
        circ = patches.Circle((xx,yy), 20)
        ax.add_patch(circ)
        
    idx = complete_large.frame == 100
    for xx, yy in zip(complete_large[idx].x, complete_large[idx].y):
        circ = patches.Circle((xx,yy), 20)
        ax.add_patch(circ)
    
    fig,ax = plt.subplots(1)
    ax.imshow(Img,cmap='gray')
    C1 = plt.Circle((complete_large[complete_large.particle == 227].x[0]-15,
                     complete_large[complete_large.particle == 227].y[0]+3),
                    30,
                    fc = 'none',
                    ec = 'red',
                    lw = 1)
    C2 = plt.Circle((complete_large[complete_large.particle == 725].x[0],
                     complete_large[complete_large.particle == 725].y[0]+25),
                    30,
                    fc = 'none',
                    ec = 'red',
                    lw = 1)
    ax.add_patch(C1)
    ax.add_patch(C2)
    plt.axis('off')
    plt.show()
    
    
    # Slip particle 227
    # Late rearrangement 725
    px_size = 5000/2020 # um/px
    x = complete_large[complete_large.particle == 227].x * px_size
    y = complete_large[complete_large.particle == 227].y * px_size
    plt.figure(dpi=500)
    plt.plot(x[50:171], y[50:171])
    plt.plot(x[170:291], y[170:291])
    plt.xlabel('X-position ($\mu m$)')
    plt.ylabel('Y-position ($\mu m$)')
    plt.xlim([1060, 1240])
    plt.ylim([1635, 1680])
    
    plt.figure(dpi=500)
    plt.plot(x[1750:1871], y[1750:1871])
    plt.plot(x[1870:1991], y[1870:1991])
    plt.xlabel('X-position ($\mu m$)')
    plt.ylabel('Y-position ($\mu m$)')
    plt.xlim([1060, 1240])
    plt.ylim([1635, 1680])
    
    plt.figure(dpi=500)
    plt.plot(x[1630:1751], y[1630:1751])
    plt.plot(x[1750:1871], y[1750:1871])
    plt.xlabel('X-position ($\mu m$)')
    plt.ylabel('Y-position ($\mu m$)')
    plt.xlim([1060, 1240])
    plt.ylim([1635, 1680])
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    