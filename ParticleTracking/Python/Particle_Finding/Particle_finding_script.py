# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 13:30:42 2020

@author: Lars Kool
"""

import os
import json
import ray
import Particle_finding_functions as Pff

#%% Indicate data locations
directory = r'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25'
version = r'V1'
files = [f for f in os.listdir(os.path.join(directory,'RawData',version)) if f[-4:] == '.tif']
output = 'Preprocessed\\'+ version
nFiles = len(files)

#%% Indicate data processing parameters
load_settings = False
verbose = True
Dev_gauss = True
Save_files = False
Parallel_processing = False
global nCores #Only used when Parallel_processing
nCores = 10

if not load_settings:
    # Create a settings dict and fill in the values below
    settings = {}
    
    settings['R'] = [[19, 22],[26, 30]] # list of lists of radii of the inner and outer circle used to create an annulus mask, every list is a particle type to detect
    settings['crop'] = [0, 2303, 158, 2178] # x,y coordinates of the top left and bottom right pixel to crop the image (X-TL,Y-TL,X-BR,Y-BR)
    settings['thresh_img'] = 0.53618 # Value to threshold the convoluted image (main way to adjust sensitivity of the particle finding algorithm, after the mask has been optimized)
    settings['thresh_conv'] = 0.2
    settings['inv_img'] = True # Invert the image before particle finding
    settings['div_gauss'] = True # Divide image by a coarse gaussian blur of the image (to correct for background lighting)
    settings['save_files'] = True
    settings['parallel_processing'] = False
    settings['nCores'] = 10
    
    selection_criteria = {} # Dict with criteria to use to select particles from possible particle locations. Every key contains a list of length Ncriteria, where identical index means same criteria.
    selection_criteria['Property'] = ['Area','Area']
    selection_criteria['Value'] = [35, 200]
    selection_criteria['criteria'] = ['Greater', 'Smaller']
    
    settings['selection_criteria'] = selection_criteria
else:
    if os.path.exist:
        settings = json.load(open(os.path.join(directory,output,r'\settings.json'),'x'))
    else:
        print('No settings file present')
        quit

# Save settings dict to json file
if Save_files:
    settings_json = open(os.path.join(directory,output,r'settings.json'),'r+')
    json.dump(settings,settings_json, indent=4)
    settings_json.close()

#%% Find Particles, below this line no input is required!

# First, create the folder structure, if not already present
Pff.check_file_structure(directory, version)

if settings['parallel_processing']:
    # Initialize the ray parallel workflow with nCores
    ray.init(num_cpus=nCores)
    # Initialize functions to be executed in parallel, the FindParallel
    # function is just a ray.remote wrapper around the FindSerial function
    results = [Pff.FindParallel.remote(directory, version, output, files[i], settings) for i in range(302)]
    # Execute the parallel functions in parallel the function itself has no
    # output. Data is saved to file directly to avoid memory overflow if
    # Save_files, otherwise data is discarded
    ray.get(results)
    # Shut down the parallel workflow
    ray.shutdown()
elif settings['verbose']:
    i = 0
    Img = Pff.image_pretreatment(directory, version, output, files[i], settings)
    Particles = Pff.find_serial(Img, directory, version, output, files[i], settings)
else:
    i = 0
    Img = Pff.image_pretreatment(directory, version, output, files[i], settings)
    Pff.find_serial(Img, directory, version, output, files[i], settings)