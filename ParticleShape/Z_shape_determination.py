# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 11:42:01 2020

@author: Lars Kool
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import os

def cart2pol(x, y):
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(r, phi)

def pol2cart(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return(x, y)

def radialAvg(r, Img, bins):
    if len(Img.shape) == 3:
        nFiles = np.shape(Img)[2]
    else:
        nFiles = 1
        
    Avg_I = np.zeros([len(bins)-1, nFiles],dtype='float64')
    for i in range(len(bins)-1):
        idx = np.array([r > bins[i], r < bins[i+1]]).all(axis=0)
        Avg_I [i,:] = np.mean(Img[idx], axis=0)
        
    Avg_centered = Avg_I - np.mean(Avg_I,axis=0)
    return Avg_centered

if __name__ == '__main__':
    path = r'F:\Lars\Particle Shape\20201015 Size Determination\20x\RawData\\'
    files = [file for file in os.listdir(path) if file.endswith('.tif')]
    nFiles = len(files)
    pixel_size = 1/(0.3974*8) #um/px
    # offset = [1182.1, 1232.8]
    offset = [[1176.505, 1148.349],[1158.660, 1424.716],[1143.457, 1191.207]]
    # bins = np.arange(175, 250, 2)
    bins = np.arange(150, 300, 0.5)
    # Load example image (to determine size)
    Img = plt.imread(path+files[0])
    Img_size = np.shape(Img)
    # Load the images and determine the radial distribution of the intensities
    # around offset
    Img = np.zeros([2304,2304,nFiles])
    Avg_I = np.zeros([len(bins)-1,nFiles])
    for i in range(0,nFiles):
        Img[:,:,i] = plt.imread(path+files[i])
        # Create meshgrid for cart -> pol conversion, and subtract the centoid of the ROI
        [xv, yv] = np.meshgrid(range(Img_size[0]),range(Img_size[1]))
        xv = xv - offset[i][0]
        yv = yv - offset[i][1]
        # Convert the xy positions of the pixels to polar coordinates
        [r, phi] = cart2pol(xv, yv)
        Avg_I[:,i] = radialAvg(r, Img[:, :, i], bins).reshape(len(bins)-1)
    
    # plt.figure(dpi=500)
    # plt.plot(bins[1:]*pixel_size,Avg_I[:,0])
    # plt.xlabel(r'r $(\mu m)$')
    # plt.ylabel('I (a.u.)')
    

    
    
        
        
        
        

    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
