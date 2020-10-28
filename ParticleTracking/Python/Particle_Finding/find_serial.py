# -*- coding: utf-8 -*-
""" Main particle finding function
Returns and/or saves centroids of particles found in img, using settings.
This algorithm uses a convolution of the image and a mask fabricated from
settings provided in settings to determine the particle locations.

Parameters
----------
img : 2D numpy array of floats
    All pixel values are expected to be normalized, such that all intensity
    values are [0-1]
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
locs_output : list of 2D numpy arrays
    List of arrays, where every item in list is a two-column numpy array.
    The first column contains x-coordinates of the particles found, and
    the second column contains the y-coordinates. The index of the list
    corresponds to the particle size belonging to settings['R'][index]
    (0,0) is in the top left corner of the image, and locations are in
    pixels.
    
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
import pandas as pd
from create_annulus_mask import create_annulus_mask
from find_particles_convolution import find_particles_convolution
from remove_overlapping_particles import remove_overlapping_particles

def find_serial(img, path, fname, settings):
    # Determine some parameters for further use
    radii = np.array(settings['R'])
    nTypes = np.shape(radii)[0]
    nFill = int(np.ceil(np.log10(nTypes)))
    selection_criteria = pd.DataFrame(settings['selection_criteria'])
    # Set all pixels below a threshold to 0 (increase contrast), and pad image
    img[img < settings['thresh_img']] = 0
    # For each particle type, perform particle convolution to find the centers
    locs_output = []
    fit_vals_output = []
    for i in range(nTypes):
        img_pad = np.pad(img, (np.max(radii[i]*5), np.max(radii[i]*5)))
        mask = create_annulus_mask(radii[i], np.shape(img_pad))
        locs, fit_vals = find_particles_convolution(img_pad,
                                                    np.max(radii[i]),
                                                    selection_criteria,
                                                    settings['thresh_conv'],
                                                    mask)
        locs_output.append(locs)
        fit_vals_output.append(fit_vals)
    
    locs_output = remove_overlapping_particles(locs_output, fit_vals_output, radii[:,1])
    for i in range(nTypes):
        nParticles = len(locs_output[i])
        # Save results as .pkl, if desired
        if settings['save_files']:
            particles = pd.DataFrame(np.array([locs_output[i][:, 0],
                                                locs_output[i][:, 1],
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
                                              path[3],
                                              'PKL',
                                              'P'+str(i+1),
                                              out_file_name))
    # Also, output the found centers as array for easy debugging
    if not settings['verbose']:
        locs_output = []
    return locs_output