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
    Settings = {}
    
    Settings['R'] = [[19, 22],[26, 30]] # list of lists of radii of the inner and outer circle used to create an annulus mask, every list is a particle type to detect
    Settings['Crop'] = [0, 2303, 158, 2178] # x,y coordinates of the top left and bottom right pixel to crop the image (X-TL,Y-TL,X-BR,Y-BR)
    Settings['thresh_img'] = 0.53618 # Value to threshold the convoluted image (main way to adjust sensitivity of the particle finding algorithm, after the mask has been optimized)
    Settings['Thresh_conv'] = 0.2
    Settings['Inv_img'] = True # Invert the image before particle finding
    Settings['Div_gauss'] = True # Divide image by a coarse gaussian blur of the image (to correct for background lighting)

    SelectionCriteria = {}
    SelectionCriteria['Property'] = ['Area','Area']
    SelectionCriteria['Value'] = [35, 200]
    SelectionCriteria['Criteria'] = ['Greater', 'Smaller']
    
    Settings['SelectionCriteria'] = SelectionCriteria
else:
    if os.path.exist:
        Settings = json.load(open(os.path.join(directory,output,r'\Settings.json'),'x'))

if Save_files:
    Settings_json = open(os.path.join(directory,output,r'Settings.json'),'r+')
    json.dump(Settings,Settings_json, indent=4)
    Settings_json.close()

#%% Find Particles, below this line no input is required!

# First, create the folder structure, if not already present
Pff.CheckFileStructure(directory, version)

if Parallel_processing:
    # Initialize the ray parallel workflow with nCores
    ray.init(num_cpus=nCores)
    # Initialize functions to be executed in parallel, the FindParallel
    # function is just a ray.remote wrapper around the FindSerial function
    results = [Pff.FindParallel.remote(i, directory, version, output, files[i], Settings) for i in range(302)]
    # Execute the parallel functions in parallel the function itself has no
    # output. Data is saved to file directly to avoid memory overflow if
    # Save_files, otherwise data is discarded
    ray.get(results)
    # Shut down the parallel workflow
    ray.shutdown()
elif verbose:
    i = 0
    Img = Pff.ImagePretreatment(directory, version, output, files[i], Settings)
    Particles = Pff.FindSerial(Img, directory, version, output, files[i], Settings, Save_files, verbose)
else:
    i = 0
    Img = Pff.ImagePretreatment(directory, version, output, files[i], Settings)
    Pff.FindSerial(Img, directory, version, output, files[i], Settings, Save_files, verbose)
    
    # if Settings['Div_gauss']:
    #     Img = Img/scipyimage.gaussian_filter(Img,sigma=50)
    #     Img = np.uint16(Img*(39050/1.0245))
        
    # if Settings['Inv_img']:
    #     Img_inv = 2**16 - Img - 1
    #     Img = Img_inv