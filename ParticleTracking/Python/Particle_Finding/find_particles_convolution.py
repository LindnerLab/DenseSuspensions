# -*- coding: utf-8 -*-
"""Convolutional based method to find particles
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
from particle_selection import particle_selection

def find_particles_convolution(img, radius, selection_criteria, threshold, mask):
    """ 
    """
    # Perform 2D Fast Fourier Transform on image and mask
    img_fft = np.fft.fft2(img)
    mask_fft = np.fft.fft2(mask)
    # Convolute and normalize image, and remove padding
    img_conv_fft = img_fft*mask_fft
    img_conv = np.fft.ifft2(img_conv_fft).real
    img_conv_norm = img_conv/np.max(img_conv)
    img_conv_norm = img_conv_norm[5*radius:-5*radius, 5*radius:-5*radius] 
    # Binarize the image, such that it can be used in regionprops
    img_conv_bin = np.copy(img_conv_norm)
    img_conv_bin[img_conv_bin < threshold] = 0
    img_conv_bin[img_conv_bin >= threshold] = 1
    # Determine locations of peaks in convoluted image
    particles, fit_vals  = particle_selection(img_conv_bin,
                                              img_conv_norm,
                                              selection_criteria,
                                              radius)
    # Correct for offset of the mask (r,r)
    particles = particles - radius
    return particles, fit_vals