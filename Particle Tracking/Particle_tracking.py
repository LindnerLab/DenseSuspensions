# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:03:48 2020

@author: lars Kool
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
from scipy import ndimage
from skimage import morphology, util, filters
import skimage
import matplotlib.patches as mpatches
from numba import jit
import multiprocessing as mp
import Pair_Correlation as pc
import os

import pims
import trackpy as tp

@pims.pipeline
def crop(img):
    """
    Crop the image to select the region of interest
    """
    x_min = 155
    x_max = 2304
    y_min = 444
    y_max = 2028 
    return img[y_min:y_max,x_min:x_max]

@pims.pipeline
def preprocess_foam(img):
    """
    Apply image processing functions to return a binary image
    """ 
    # Apply thresholds
    gauss = ndimage.gaussian_filter(np.float64(img),100)
    img_filt = np.subtract(img,gauss)
    threshold = 5000
    idx = img_filt < threshold
    idx2 = img_filt >= threshold
    img_filt[idx] = 0
    img_filt[idx2] = 1
    
    # Crop the pictures as for raw images.
    img_filt = crop(img_filt)
    return img_filt

def main():
    pass

if __name__ == "__main__":
    datapath = 'E:/Lars/Failure/20200310/v2/Preprocessed/'
    prefix = 'Data00002_'
    files = [file for file in os.listdir(datapath) if file.endswith('.pkl')]
    nFiles = len(files)
    
    result = pd.DataFrame(columns=['y', 'x', 'r','frame'])
    nParticles = []
    
    for idx, file in enumerate(files):
        nParticles.append(len(pd.read_pickle(datapath + file)))
        result = pd.concat([result, pd.read_pickle(datapath + file)], ignore_index = True)
        
    search_range = 60
    tracked = tp.link_df(result, search_range, memory=0)
    tracked = tracked.sort_index()
    
        # Optimize this, as it is incredibly slow!!!
    tracked_complete = pd.DataFrame(columns = ['y','x','r','frame','particle'])
    for i in range(np.max(np.array(tracked.particle))):
        if np.sum(np.array(tracked.particle) == i) == nFiles:
            tracked_complete = pd.concat([tracked_complete, pd.concat([tracked.iloc[[j]] for j in tracked.index if tracked.particle[j] == i])])
    tracked_complete = tracked_complete.reset_index(drop=True)
    
    
    # img = preprocess_foam(pims.open('E:/Lars/Failure/20200310/v2/Raw_Data/Data00002_00500.tif'))
    # img = img[0]
    # fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(14, 14))
    # tp.plot_traj(t,superimpose=img)
    # plt.savefig('E:/Lars/Failure/20200310/v2/Processed/Track_v1.png',dpi=1000)
            
            

    
    # x = Centers['x'].values
    # y = Centers['y'].values
    # S = 1500
    # rMax = 300
    # dr = 10
    
    # (g_average, radii, interior_indices) = pc.pairCorrelationFunction_2D(x, y, S, rMax, dr)
    # plt.plot(radii/90, g_average)
    # plt.xlabel('d/Dsmall')
    # plt.ylabel('g(r)')
    # plt.savefig('E:/Lars/Failure/20200310/v2/Processed/g_r.png',dpi=1000)