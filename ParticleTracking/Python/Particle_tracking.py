# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:03:48 2020

@author: lars Kool
"""
import numpy as np
import pandas as pd
import os
from Track_repair import trackRepair

def DataParser(datapath, dataFormat):
    result = pd.DataFrame(columns=['y', 'x', 'r','frame'])
    nParticles = []
    
    if dataFormat == 'csv':
        files = [file for file in os.listdir(datapath) if file.endswith('.csv')]
        nFiles = len(files)
        
        for idx, file in enumerate(files):
            particles_temp = pd.read_csv(datapath + file, names=['x','y', 'r','frame'])
            particles_temp.frame = idx
            nParticles.append(len(particles_temp))
            result = pd.concat([result, particles_temp], ignore_index = True)
    
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
    datapath = 'F:/Lars/Oscillatory Compression/20200713 Soft Particle/Avg100_Amp80_Per120_Back25/Preprocessed/'
    prefix = '0'
    dataFormat = 'csv'
    minLength = 1000
    
    # Import and parse data
    [imported, nFiles] = DataParser(datapath+'CSV', dataFormat)
    
    # Separate small from large particles, this makes tracking faster (as there are less permutations) and more reliable (as small particles should never by matched with large ones and viceversa)
    Psmall = imported[imported.r == 20.5].reset_index(drop=True)
    Plarge = imported[imported.r == 28].reset_index(drop=True)
    
    # Track the particles
    search_range = 10
    tracked_small = tp.link_df(Psmall, search_range, memory=5).sort_values(by=['particle','frame'], ignore_index=True)
    tracked_large = tp.link_df(Plarge, search_range, memory=5).sort_values(by=['particle','frame'], ignore_index=True)

    complete_small = trackRepair(tracked_small,5,1384, minLength)
    complete_large = trackRepair(tracked_large,5,1384, minLength)
    
    complete_small.to_pickle(datapath+'Complete_tracks_small.pkl')
    complete_large.to_pickle(datapath+'Complete_tracks_large.pkl')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    