# -*- coding: utf-8 -*-
""" Particle finding functions
The functions in this file are all part of the particle finding algorithm
developed by the contributors mentioned below.

Contributors : Lars Kool
Affiliations : Laboratoire Physiqe et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France
               
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162               
"""

import os
import sys
import pandas as pd
import numpy as np
import scipy.ndimage as scipyimage
import matplotlib as mpl
import matplotlib.pyplot as plt
import skimage.measure as measure
import ray

#%%
def check_file_structure(path):
    """ This function checks if all the input and output folders are present.
    And if not, it will create them.
    
    Parameters
    ----------
    path : List of strings
        List of foldernames:
            Path[0]: General folder of the data to be analyzed \n
            Path[1]: Name of the folder containing the input images \n
            Path[2]: Name of the output folder \n
            Path[3]: Name of the version of the dataset (allows to distinguish
                     multiple datasets taken on the same day with identical
                     settings)

    Returns
    -------
    None.

    """
    if not os.path.isdir(os.path.join(path[0], path[2])):
        os.mkdir(os.path.join(path[0], path[2]))
        os.mkdir(os.path.join(path[0], ''.join(path[2:4])))
        os.mkdir(os.path.join(path[0], ''.join(path[2:4]), 'PKL'))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], ''.join(path[2:4]))):
        os.mkdir(os.path.join(path[0], ''.join(path[2:4])))
        os.mkdir(os.path.join(path[0], ''.join(path[2:4]), 'PKL'))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], ''.join(path[2:4]),
                                        'PKL')):
        os.mkdir(os.path.join(path[0], ''.join(path[2:4]), 'PKL'))
        print('Output File structure created')
    else:
        print('Output File structure already present')
        
    if not os.path.isdir(os.path.join(''.join(path[0:2]), path[3])):
        print('''Input files not present in the expected folder:
              DIRECTORY\\INPUT\\VERSION)''')
        sys.exit()
        
#%%
def image_pretreatment(path, fname, settings):
    """ Pretreatment of an image.
    Prepare a raw image for the particle finding algorithm. Current treatment
    options are:\n
    - Divide by gaussian blur
    - Invert image              

    Parameters
    ----------
    path : List of strings
        List of foldernames:
            Path[0]: General folder of the data to be analyzed \n
            Path[1]: Name of the folder containing the input images \n
            Path[2]: Name of the output folder \n
            Path[3]: Name of the version of the dataset (allows to distinguish
                     multiple datasets taken on the same day with identical
                     settings)
    fname : String
        Name of the file (including extension)to be preprocessed that can be
        found in folder: path[0]\\path[1]\\path[3]
    settings : dict
        Dict with settings to be used in the data analysis. See main particle
        finding script for a detailed explanation of all keys.

    Returns
    -------
    img : 2D numpy array with identical shape to input image.
        All pixel values are normalized, so that all intensity values are [0-1]
    """
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
    """
    

    Parameters
    ----------
    img : TYPE
        DESCRIPTION.
    path : List of strings
        List of foldernames:
            Path[0]: General folder of the data to be analyzed \n
            Path[1]: Name of the folder containing the input images \n
            Path[2]: Name of the output folder \n
            Path[3]: Name of the version of the dataset (allows to distinguish
                     multiple datasets taken on the same day with identical
                     settings)
    fname : String
        Name of the file (including extension)to be preprocessed that can be
        found in folder: path[0]\\path[1]\\path[3]
    settings : dict
        Dict with settings to be used in the data analysis. See main particle
        finding script for a detailed explanation of all keys.

    Returns
    -------
    locs_output : list of 2D numpy arrays
        List of arrays, where every item in list is a two-column numpy array.
        The first column contains x-coordinates of the particles found, and
        the second column contains the y-coordinates. The index of the list
        corresponds to the particle size belonging to settings['R'][index]
        (0,0) is in the top left corner of the image, and locations are in
        pixels.
    """
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
        mask = create_annulus_mask(radii[i], np.shape(img))
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
    """ Parallized version of find_serial. It is parallelized using the ray
    library https://github.com/ray-project/ray
    
    Parameters
    ----------
    path : List of strings
        List of foldernames:
            Path[0]: General folder of the data to be analyzed \n
            Path[1]: Name of the folder containing the input images \n
            Path[2]: Name of the output folder \n
            Path[3]: Name of the version of the dataset (allows to distinguish
                     multiple datasets taken on the same day with identical
                     settings)
    fname : String
        Name of the file (including extension) to be preprocessed that can be
        found in folder: path[0]\\path[1]\\path[3]
    settings : dict
        Dict with settings to be used in the data analysis. See main particle
        finding script for a detailed explanation of all keys.

    Returns
    -------
    None.

    """
    # Set verbose to False (parallel doesn' work with verbose), pretreat image
    # and then apply the 'normal' find_serial to the pretreated image
    settings['verbose'] = False
    img = image_pretreatment(path, fname, settings)
    find_serial(img, path, fname, settings)
    
#%%
def create_annulus_mask(radii, img_size):
    """ Creation of an annulus mask.
    Entire mask size is img_size, and outer radius of annulus is max(radii)
    and inner radius is min(radii)

    Parameters
    ----------
    radii : List of length=2 of numerics (float, int)
        List of radii to describe the annulus (max radius is the outer ring
        and min radius is the inner ring)
    img_size : List of length=2
        List containing (width, height) of the image to which this mask will
        be applied. Size is given in pixels.

    Returns
    -------
    mask : 2D numpy array
        2D numpy array of width img_size[0] and height img_size[1], where
        the background pixels are 0, the annulus pixels are 1 and the pixels
        inside the annulus are -0.5 (this prevents masks of bigger particles
        to easily find smaller particles)
    """
    # Find center of object, (r,r) for a circle meaning that the peaks in the
    # convoluted image will be at (x+r,y+r), instead of (x,y)! Also, note if
    # the first radius is the inner or outer circle of the annulus
    x_center = np.max(radii)
    y_center = np.max(radii)
    idx = np.where(radii == x_center)[0][0]
    # Determine which pixels are inside the circles
    mask_temp = np.zeros([img_size[0], img_size[1], 2])
    nPx = img_size[0]*img_size[1]
    for i in range(2):
        x_circle = radii[i]*np.cos(np.linspace(0, 2*np.pi, num=16))+ x_center
        y_circle = radii[i]*np.sin(np.linspace(0, 2*np.pi, num=16)) + y_center
        xv, yv = np.meshgrid(range(img_size[1]), range(img_size[0]))
        path_circle = mpl.path.Path(np.transpose([x_circle, y_circle]))
        mask_temp[:, :, i] = path_circle.contains_points(np.transpose([xv.reshape(nPx),
                                                                       yv.reshape(nPx)])
                                                         ).reshape((img_size[0],
                                                                    img_size[1]))
    # Build the mask such that I(outside) = 0, I(annulus) = 1, and
    # I(inside) = -0.5
    if idx == 0:
        mask = mask_temp[:, :, 0] - 1.5*mask_temp[:, :, 1]
    else:
        mask = mask_temp[:, :, 1] - 1.5*mask_temp[:, :, 0]
    return mask

#%%
def find_particles_convolution(img, radius, selection_criteria, threshold, mask):
    """ Convolutional based method to find particles
    Algorith that uses FFT convolution to find rotationally symmetric
    particles using a mask.    

    Parameters
    ----------
    img : 2D numpy array of float
        Pre-treated image (0 <= pixel values <= 1), in which particles of
        radius radius are present
    radius : int
        Radius of the particles to be found (here only used to remove padding)
    selection_criteria : dict
        Dict with criteria that will be used to distinguish particles from
        false positives.
        Contains three keys (property, value, criteria), all of which are
        associated with a list of length nCriteria.\n
        properties : List of regionprops properties to be used to distinguish
            particles and false positives.\n
        value : List of threshold values to which to compare the regionprops
            properties.\n
        criteria : List of 'greater' or 'smaller' strings. Used to indicate if
            the regionprops properties should be greater or smaller than value.\n
    threshold : float
        Value to convert the convoluted image into a binary image, which can
        then be used in regionprops. As the convoluted image is normalized,
        0 <= threshold <= 1.
    mask : 2D numpy array of numerics (float, int)
        Mask of identical size as img which is used as a template to find
        particles. Features to be found should correspond to positive
        intensities, background to 0, and features to exclude should
        correspond to negative values.

    Returns
    -------
    particles : 2D numpy array of floats
        2-column numpy array of floats, where the first column corresponds to
        the x-coordinates of the particles found, and the second column to the
        y-coordinates of the pixels found. (0,0) is the top-left corner of the
        image and all distances are in pixels
    """
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
    """ Convert a binary image into a list of particle coordinates that
    comply with the criteria in 'selection_criteria' using regionprops.   

    Parameters
    ----------
    img_conv_bin : 2D numpy array (binary values only)
        Binarized version of the image after convolution with a mask. Used
        with scipy.ndimage.measure regionprops functions
    img_conv_norm : 2D numpy array of float
        Normalized image after convolution (0 <= I <= 1) used to determine the
        weighted centroid of particles.
    selection_criteria : dict
        Dict with criteria that will be used to distinguish particles from
        false positives.
        Contains three keys (property, value, criteria), all of which are
        associated with a list of length nCriteria.\n
        properties : List of regionprops properties to be used to distinguish
            particles and false positives.\n
        value : List of threshold values to which to compare the regionprops
            properties.\n
        criteria : List of 'greater' or 'smaller' strings. Used to indicate if
            the regionprops properties should be greater or smaller than value.\n
    threshold : float
        Value to convert the convoluted image into a binary image, which can
        then be used in regionprops. As the convoluted image is normalized,
        0 <= threshold <= 1.
    mask : 2D numpy array of numerics (float, int)
        Mask of identical size as img which is used as a template to find
        particles. Features to be found should correspond to positive
        intensities, background to 0, and features to exclude should
        correspond to negative values.

    Returns
    -------
    particles : 2D numpy array of floats
        2-column numpy array of floats, where the first column corresponds to
        the x-coordinates of the particles found, and the second column to the
        y-coordinates of the pixels found. (0,0) is the top-left corner of the
        image and all distances are in pixels./n
        Or, if no regions met the conditions/no regions were found an empty
        array will be returned and an appropriate message will be printed
    """
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
        particles = np.array([])
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
            particles = np.array([])
            
    return particles
    
#%% Debugging
if __name__ == '__main__':
    pass
