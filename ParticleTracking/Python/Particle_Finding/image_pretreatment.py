# -*- coding: utf-8 -*-
"""Pretreatment of an image.
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
    found in folder: path[0]/path[1]/path[3]
settings : dict
    Dict with settings to be used in the data analysis. See main particle
    finding script for a detailed explanation of all keys.

Returns
-------
img : 2D numpy array with identical shape to input image.
    All pixel values are normalized, so that all intensity values are [0-1]
        
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

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as scipyimage

def image_pretreatment(path, fname, settings):
    # Only pretreat image if the image isn't already analyzed
    if not os.path.isfile(os.path.join(path[0],
                                       path[2],
                                       'PKL',
                                       'P1'+fname[0:-4]+'.pkl')):
        crop = np.array(settings['crop'])
        img = plt.imread(os.path.join(path[0],
                                      path[1],
                                      path[3],
                                      fname))
        img = img[crop[0]:crop[1]][crop[2]:crop[3]]
        # Divide image by gaussian blur of image to correct for lumination
        # differences across the sample
        if settings['div_gauss']:
            img = img/scipyimage.gaussian_filter(img, sigma=50)
            img = img/np.max(img)
        # Invert the image
        if settings['inv_img']:
            img = 1-img
            
    return img