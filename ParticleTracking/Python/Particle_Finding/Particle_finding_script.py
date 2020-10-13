# -*- coding: utf-8 -*-
"""
This is the main script belonging to the particle finding algorith. This
script 

@author: Lars Kool
"""

import os
import sys
import json
import ray
import particle_finding_functions as Pff

#%% Indicate data locations
DIRECTORY = r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25'
VERSION = r'V1'
INPUT = r'RawData'
OUTPUT = r'Preprocessed'
PATH = [DIRECTORY, INPUT, OUTPUT, VERSION]
files = [f for f in os.listdir(os.path.join(PATH[0:2])) if f[-4:] == '.tif']
nFiles = len(files)
load_settings = False
# Set file of interest to analyse verbosely (Set settings['verbose'] to True!)
file_of_interest = 0

#%% Indicate data processing parameters
if not load_settings:
    # Create a settings dict and fill in the values below
    settings = {}
    # list of lists of radii of the inner and outer circle used to create an
    # annulus mask, every list is a particle type to detect
    settings['R'] = [[19, 22], [26, 30]]
    # x,y coordinates of the top left and bottom right pixel to crop the image
    # (X-TL,Y-TL,X-BR,Y-BR)
    settings['crop'] = [0, 2303, 158, 2178]
    # Value to threshold the image before convolution (main way to adjust
    # sensitivity of the particle finding algorithm, after the mask has been
    # optimized)
    settings['thresh_img'] = 0.536
    # Threshold used for binarization of the convoluted image to prep it for
    # regionprops. Typically, the SN-ratio after convolution is quite high,
    # so tinkering with this value isn't needed very often. Generally,
    # lower = better, as this will give larger regions to calculate the centroid
    settings['thresh_conv'] = 0.2
    # Invert the image before particle finding
    settings['inv_img'] = True
    # Divide image by a coarse gaussian blur of the image (to correct for
    # background illumination inhomogeneities)
    settings['div_gauss'] = True
    settings['save_files'] = True
    settings['parallel_processing'] = False
    settings['nCores'] = 10
    # Dict with criteria to use to select particles from possible particle
    # locations. Every key contains a list of length Ncriteria, where identical
    # index means identical criterium.
    selection_criteria = {}
    selection_criteria['property'] = ['Area', 'Area']
    selection_criteria['value'] = [35, 200]
    selection_criteria['criteria'] = ['greater', 'smaller']
    settings['selection_criteria'] = selection_criteria
else:
    if os.path.exist:
        settings = json.load(open(os.path.join(PATH[0],
                                               PATH[2],
                                               r'\settings.json')
                                  , 'x'))
    else:
        print('No settings file present')
        sys.exit()
        
#%% Find Particles, below this line no input is required!
# First, create the folder structure, if not already present
Pff.check_file_structure(PATH)

if settings['parallel_processing']:
    # Initialize the ray parallel workflow with nCores
    ray.init(num_cpus=settings['nCores'])
    # Initialize functions to be executed in parallel, the FindParallel
    # function is just a ray.remote wrapper around the FindSerial function
    results = [Pff.FindParallel.remote(PATH,
                                       files[i],
                                       settings)
               for i in range(nFiles)]
    # Execute the parallel functions in parallel the function itself has no
    # output. Data is saved to file directly to avoid memory overflow if
    # Save_files, otherwise data is discarded
    ray.get(results)
    ray.shutdown()
elif settings['verbose']:
    img = Pff.image_pretreatment(PATH,
                                 files[file_of_interest],
                                 settings)
    particles = Pff.find_serial(img,
                                PATH,
                                files[file_of_interest],
                                settings)
else:
    for i in range(nFiles):
        img = Pff.image_pretreatment(PATH,
                                     files[i],
                                     settings)
        Pff.find_serial(img,
                        PATH,
                        files[i],
                        settings)
# Save settings dict to json file
if settings['save_files']:
    settings_json = open(os.path.join(PATH[0],
                                      PATH[2],
                                      r'settings.json'),
                         'r+')
    json.dump(settings, settings_json, indent=4)
    settings_json.close()
    