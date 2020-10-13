# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 15:15:25 2020

@author: Lars Kool
"""
import os
import pandas as pd
import numpy as np
import scipy.ndimage as scipyimage
import matplotlib as mpl
import matplotlib.pyplot as plt
import skimage.measure as measure
import ray

#%%
def check_file_structure(path):
    '''
    Insert description here
    '''
    if not os.path.isdir(os.path.join(path[0], path[2])):
        os.mkdir(os.path.join(path[0], path[2]))
        os.mkdir(os.path.join(path[0], path[2:4]))
        os.mkdir(os.path.join(path[0], path[2:4], 'PKL'))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], path[2:4])):
        os.mkdir(os.path.join(path[0], path[2:4]))
        os.mkdir(os.path.join(path[0], path[2:4], 'PKL'))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], path[2:4],
                                        'PKL')):
        os.mkdir(os.path.join(path[0], path[2:4], 'PKL'))
        print('Output File structure created')
    else:
        print('Output File structure already present')
        
#%%
def image_pretreatment(path, fname, settings):
    '''
    Insert description here
    '''
    # Only pretreat image if the image isn't already analyzed
    if not os.path.isfile(os.path.join(path[0],
                                       path[2],
                                       'PKL',
                                       'P1'+fname[0:-4]+'.pkl')):
        crop = np.array(settings['crop'])
        img = plt.imread(os.path.join(path[0:2],
                                      path[3],
                                      fname))
        img = img[crop[0]:crop[1]][crop[2]:crop[3]]
        # Divide image by gaussian blur of image to correct for lumination
        # differences across the sample
        if settings['Div_gauss']:
            img = img/scipyimage.gaussian_filter(img, sigma=50)
            img = img/np.max(img)
        # Invert the image
        if settings['Inv_img']:
            img = 1-img
            
    return img

#%%
def find_serial(img, path, fname, settings):
    '''
    Insert description here
    '''
    # Determine some parameters for further use
    radii = np.array(settings['R'])
    nTypes = np.shape(radii)[0]
    nFill = int(np.ceil(np.log10(nTypes)))
    selection_criteria = pd.DataFrame(settings['selection_criteria'])
    # Set all pixels below a threshold to 0 (increase contrast), and pad image
    img[img < settings['thresh']] = 0
    img = np.pad(img, (np.max(radii), np.max(radii)))
    # For each particle type, perform particle convolution to find the centers
    locs_output = []
    for i in range(nTypes):
        mask = annulus_mask(radii[i], np.shape(img))
        locs = find_particles_convolution(img,
                                          np.max(radii[i]),
                                          selection_criteria,
                                          settings['thresh_conv'],
                                          mask)
        nParticles = len(locs)
        # Save results as .pkl, if desired
        if settings['save_files']:
            particles = pd.DataFrame(np.array([locs[:, 0],
                                               locs[:, 1],
                                               np.ones((nParticles,))*np.mean(radii[i])]
                                              ).transpose(),
                                     columns=['x', 'y', 'r'])
            out_file_name = ''.join(['P',
                                     str(i+1).zfill(nFill),
                                     '_',
                                     fname[0:-4],
                                     '.pkl'])
            particles.to_pickle(os.path.join(path[0],
                                             path[2],
                                             'PKL',
                                             out_file_name))
        # Also, output the found centers as array for easy debugging
        if settings['verbose']:
            locs_output.append([locs])
    return locs_output

#%%
@ray.remote
def find_parallel(path, fname, settings):
    '''
    Insert description here
    '''
    # Set verbose to False (parallel doesn' work with verbose), pretreat image
    # and then apply the 'normal' find_serial to the pretreated image
    settings['verbose'] = False
    img = image_pretreatment(path, fname, settings)
    find_serial(img, path, fname, settings)
    
#%%
def annulus_mask(radius, img_size):
    '''
    Insert description here
    '''
    # Find center of object, (r,r) for a circle meaning that the peaks in the
    # convoluted image will be at (x+r,y+r), instead of (x,y)! Also, note if
    # the firste radius is the inner or outer circle of the annulus
    x_center = np.max(radius)
    y_center = np.max(radius)
    idx = np.where(radius == x_center)[0][0]
    # Determine which pixels are inside the circles
    mask_temp = np.zeros([img_size[0], img_size[1], 2])
    nPx = img_size[0]*img_size[1]
    for i in range(2):
        x_circle = radius[i]*np.cos(np.linspace(0, 2*np.pi, num=16))+ x_center
        y_circle = radius[i]*np.sin(np.linspace(0, 2*np.pi, num=16)) + y_center
        xv, yv = np.meshgrid(range(img_size[1]), range(img_size[0]))
        
        path_circle = mpl.path.Path(np.transpose([x_circle, y_circle]))
        mask_temp[:, :, i] = path_circle.contains_points(np.transpose([xv.reshape(nPx),
                                                                       yv.reshape(nPx)])
                                                         ).reshape((img_size[0],
                                                                    img_size[1]))       # Check for each x and y in the mesh if it is inside or outside the 'object' on the mask.
    # Build the mask such that I(outside) = 0, I(annulus) = 1, and
    # I(inside) = -0.5
    if idx == 0:
        mask = mask_temp[:, :, 0] - 1.5*mask_temp[:, :, 1]
    else:
        mask = mask_temp[:, :, 1] - 1.5*mask_temp[:, :, 0]
    return mask

#%%
def find_particles_convolution(img, radius, selection_criteria, threshold, mask):
    '''
    Insert description here
    '''
    # Perform 2D Fast Fourier Transform on image and mask
    img_fft = np.fft.fft2(img)
    mask_fft = np.fft.fft2(mask)
    # Convolute and normalize image, and remove padding
    img_conv_fft = img_fft*mask_fft
    img_conv = np.fft.ifft2(img_conv_fft).real
    img_conv_norm = img_conv/np.max(img_conv)
    img_conv_norm = img_conv_norm[radius:-radius, radius:-radius] 
    # Binarize the image, such that it can be used in regionprops
    img_conv_bin = img_conv_norm
    img_conv_bin[img_conv_bin < threshold] = 0
    img_conv_bin[img_conv_bin >= threshold] = 1
    # Determine locations of peaks in convoluted image
    particles = particle_selection(img_conv_bin,
                                   img_conv_norm,
                                   selection_criteria)
    return particles

#%%
def particle_selection(img_conv_bin, img_conv_norm, selection_criteria):
    '''
    Insert description here
    '''
    # Determine the number of criteria given
    nCriteria = len(selection_criteria)    
    # Determine regionprops of all labelled regions
    img_label = measure.label(img_conv_bin)
    props = measure.regionprops(img_label, intensity_image=img_conv_norm)
    nRegions = np.shape(props)[0]
    # Extract centroids from regionprops
    centroids = np.array([prop.weighted_centroid for prop in props])

    if nCriteria == 0:
        print('No selection criteria were given,so all centroids are returned.\n')
        particles = centroids
    elif nRegions == 0:
        print('No particles were found.\n')
        particles = []
    else:
        # Determine if each region complies with all selection criteria
        selection = np.zeros((nCriteria, nRegions)).astype(bool)
        for i in range(nCriteria):
            property_ = selection_criteria.property_[i]
            value = selection_criteria.value[i]
            criteria = selection_criteria.criteria[i]
            region_props = np.array([prop[property_] for prop in props])
            # Adds a vector with logical values to Selection, where 1 means the
            # regionproperty of that region complies with the given criteria
            # and 0 means that the regionproperty doesnt comply with the
            # criteria.
            if criteria.lower() == 'greater':
                selection[i, :] = np.array([region_props > value])
            elif criteria.lower() == 'smaller':
                selection[i, :] = np.array([region_props < value])
        # AND gate with nCriteria inputs, output == 1 if all inputs are 1, else output == 0
        particles = centroids[np.prod(selection, axis=0).astype(bool)][:]
        if np.shape(particles)[1] == 0:
            print('No particles met the requirements')
            particles = []
            
    return particles
    
#%% Debugging
if __name__ == '__main__':
    pass
