# -*- coding: utf-8 -*-
"""
Convert a binary image into a list of particle coordinates that
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
    
Author information
------------------
Contributors : Lars Kool
Affiliations : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France
               
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162
More info on this programme: https://caliper-itn.org/
"""

import numpy as np
import skimage.measure as measure

def particle_selection(img_conv_bin, img_conv_norm, selection_criteria, radius):
    # Determine the number of criteria given
    nCriteria = len(selection_criteria)
    # Determine regionprops of all labelled regions
    img_label = measure.label(img_conv_bin)
    props = np.array(measure.regionprops(img_label, intensity_image=img_conv_norm))
    centroids = np.array([prop.weighted_centroid for prop in props])
    _temp = np.all(centroids > radius*1.5, axis=1)
    props = props[_temp]
    nRegions = np.shape(props)[0]
    # Extract centroids from regionprops
    centroids = centroids[_temp]
    prop_intensity = np.array([[prop.max_intensity,
                                prop.mean_intensity] for prop in props])

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
        fit_vals = prop_intensity[np.prod(selection, axis=0).astype(bool)][:]
        if np.shape(particles)[1] == 0:
            print('No particles met the requirements')
            particles = np.array([])
            fit_vals = np.array([])
    return particles, fit_vals