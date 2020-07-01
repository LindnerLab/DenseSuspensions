# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

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
import os

import pims
import trackpy as tp

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

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

def find_particles(img, num, name):
    features = pd.DataFrame()
    # Label elements on the picture
    label_image, number_of_labels = skimage.measure.label(img, background=0, return_num=True)       
    for region in skimage.measure.regionprops(label_image,img):
        # Everywhere, skip small and large areas
        if region.area < 1000 or region.area > 100000:
            continue
        # Only black areas
        if region.major_axis_length/region.minor_axis_length > 1.5:
            continue
        # On the top, skip small area with a second threshold
        if region.centroid[0] < 260 and region.area < 80:
            continue
        # Store features which survived to the criterions
        r = region.bbox
        r1 = r[2]-r[0]
        r2 = r[3]-r[1]
        features = features.append([{'y': region.centroid[0],
                                     'x': region.centroid[1],
                                     'r': np.ceil(np.mean([r1,r2])),
                                     'frame': num,
                                     },])
    features.to_pickle(name + str(num).zfill(5) + '.pkl')

def find_particles_parallel(num):
    import numpy as np
    import pandas as pd
    import skimage
    import pims
    
    filename = 'Data00002_'
    inp_loc = 'E:/Lars/Failure/20200310/v2/Raw_Data/'
    exp_loc = 'E:/Lars/Failure/20200310/v2/Preprocessed/'
    
    if filename + str(num+1).zfill(5) + '.pkl' in os.listdir(exp_loc):
        pass
    else:
        features = pd.DataFrame()
        
        img = preprocess_foam(pims.open(inp_loc + filename + str(num+1).zfill(5) + '.tif'))
        
        label_image, number_of_labels = skimage.measure.label(img[0], background=0, return_num=True)       
        for region in skimage.measure.regionprops(label_image,img[0]):
            # Everywhere, skip small and large areas
            if region.area < 1000 or region.area > 100000:
                continue
            # Only black areas
            if region.major_axis_length/region.minor_axis_length > 1.5:
                continue
            # On the top, skip small area with a second threshold
            if region.centroid[0] < 260 and region.area < 80:
                continue
            # Store features which survived to the criterions
            r = region.bbox
            features = features.append([{'y': region.centroid[0],
                                         'x': region.centroid[1],
                                         'r': np.ceil(np.mean([r[2]-r[0],r[3]-r[1]])),
                                         'frame': num +1
                                         },])
        features.to_pickle(exp_loc + filename + str(num+1).zfill(5) + '.pkl')
    return
    
def save_parallel(features, name, num):
    features.to_pickle(name + str(num).zfill(5) + '.pkl')
    return

if __name__ == "__main__":
    pool = mp.Pool(processes = 16)
    pool.map_async(find_particles_parallel, range(0,1378))
    
#    datapath = 'E:/Lars/Failure/20200310/v2/Raw_Data/'
#    prefix = 'Data00002_0000'
#    frames = preprocess_foam(pims.open(datapath + prefix + '*.tif'))
#    name = datapath + prefix
    
#    pool = mp.Pool(processes = 5)
#    with mp.Pool(processes=5) as pool:
#        pool.map(lambda num: find_particles_parallel(num), range(11))
        
#    client = ipp.Client()
#    view = client.load_balanced_view()
    
#    curried_locate = lambda image, num: find_particles_parallel(image, num, name)
#    view.map(curried_locate, frames[:-1:5], range(1,len(frames),5))
#    amr = view.map_async(curried_locate, frames[:10], range(1,len(frames),5))
#    amr.wait_interactive()
#    results = amr.get()
#    
#    features = pd.concat(results, ignore_index=True)
#    features.head()
#
#    img = frames[0]
#    search_range = 75
#    t = tp.link_df(features, search_range, memory=2)
#    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(14, 14))
#    tp.plot_traj(t,superimpose=img)
#    plt.savefig('C:/Users/Charly/Pictures/test.png',dpi=1000)
    
    
#    # Pair correlation function
#    x = t.x[0:365]
#    y = t.y[0:365]
#    S = 1500
#    rMax = 300
#    dr = 10
#    
#    (g_average, radii, interior_indices) = pc.pairCorrelationFunction_2D(x, y, S, rMax, dr)
#    plt.plot(radii/90, g_average)
#    plt.xlabel('d/Dsmall')
#    plt.ylabel('g(r)')
#    plt.savefig('E:/Lars/Failure/20200310/v2/g_r.png',dpi=1000)
    
    
    
    
    
    
    
    
    
    
    
    
    