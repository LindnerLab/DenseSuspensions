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
def CheckFileStructure(directory, version):
    if not os.path.isdir(os.path.join(directory,'Preprocessed')):
        os.mkdir(os.path.join(directory,'Preprocessed'))
        os.mkdir(os.path.join(directory,'Preprocessed',version))
        os.mkdir(os.path.join(directory,'Preprocessed',version,'PKL'))
        print('Output File structure created')
        
    elif not os.path.isdir(os.path.join(directory,'Preprocessed',version)):
        os.mkdir(os.path.join(directory,'Preprocessed',version))
        os.mkdir(os.path.join(directory,'Preprocessed',version,'PKL'))
        print('Output File structure created')
    
    elif not os.path.isdir(os.path.join(directory,'Preprocessed',version,'PKL')):
        os.mkdir(os.path.join(directory,'Preprocessed',version,'PKL'))
        print('Output File structure created')
    else:
        print('Output File structure already present')



#%%
def ImagePretreatment(directory, version, output, fname, Settings):
    if not os.path.isfile(os.path.join(directory,output,'PKL','P1'+fname[0:-4]+'.pkl')):
        Crop = np.array(Settings['Crop'])
        Img = plt.imread(os.path.join(directory, 'RawData', version, fname))[Crop[0]:Crop[1]][Crop[2]:Crop[3]]
        
        if Settings['Div_gauss']:
            Img = Img/scipyimage.gaussian_filter(Img,sigma=50)
            Img = Img/np.max(Img)
        
        if Settings['Inv_img']:
            Img = 1-Img
            
    return Img



#%%
def FindSerial(Img, directory, version, output, fname, Settings, Save_files, verbose):
    R = np.array(Settings['R'])
    nTypes = np.shape(R)[0]
    SelectionCriteria = pd.DataFrame(Settings['SelectionCriteria'])
    
    Img[Img < Settings['thresh']] = 0
    Img = np.pad(Img,(np.max(R),np.max(R)))
    Img_size = np.shape(Img)
    
    if verbose:
        Locs_output = []
    
    for i in range(nTypes):
        mask = AnnulusMask(R[i], Img_size)
        Locs = FindParticlesConvolution(Img,np.max(R[i]),SelectionCriteria,Settings['thresh_conv'],mask)
        Nparticles = len(Locs)
    
        if Save_files:
            P = pd.DataFrame(np.array([Locs[:,0], Locs[:,1], np.ones((Nparticles,))*np.mean(R[i])]).transpose(),columns=['x', 'y', 'r'])
            P.to_pickle(os.path.join(directory,
                                     output,
                                     'PKL',
                                     'P' + str(i+1).zfill(int(np.ceil(np.log10(nTypes)))) + '_' + fname[0:-4]+'.pkl'))
        
        if verbose:
            Locs_output.append([Locs])
    if verbose:
        return Locs_output



#%%
@ray.remote
def FindParallel(directory, version, output, fname, Settings, Save_files):
    verbose = False
    Img = ImagePretreatment(directory, version, fname, Settings)
    FindSerial(Img, directory, version, output, fname, Settings, Save_files, verbose)
    
    
    
#%%
def AnnulusMask(r, Img_size):
    xCenter = np.max(r)                                                         # X-position of the center of the object to look for. It is centered at (r,r), meaning that the peaks in the convoluted image will be at (x+r,y+r), instead of (x,y)!
    yCenter = np.max(r)                                                         # Y-position of the center of the object to look for.
    idx = np.where(r == xCenter)[0][0]
    
    mask_temp = np.zeros([Img_size[0],Img_size[1],2])
    nPx = Img_size[0]*Img_size[1]
    
    for i in range(2):
        xCircle = r[i]*np.cos(np.linspace(0,2*np.pi,num=16))+ xCenter           # X-coordinates of a circle centered at (xCenter, yCenter) with radius r.
        yCircle = r[i]*np.sin(np.linspace(0,2*np.pi,num=16)) + yCenter          # Y-coordinates of a circle centered at (xCenter, yCenter) with radius r.
        xv, yv = np.meshgrid(range(Img_size[1]),range(Img_size[0]))
        
        p = mpl.path.Path(np.transpose([xCircle,yCircle]))
        mask_temp[:,:,i] = p.contains_points(np.transpose([xv.reshape(nPx),yv.reshape(nPx)])).reshape((Img_size[0],Img_size[1]))       # Check for each x and y in the mesh if it is inside or outside the 'object' on the mask.
    
    if idx == 0:
        mask = mask_temp[:,:,0] - 1.5*mask_temp[:,:,1];
    else:
        mask = mask_temp[:,:,1] - 1.5*mask_temp[:,:,0];
    return mask



#%%
def FindParticlesConvolution(Img,r,SelectionCriteria, threshold, mask):
    Img_fft = np.fft.fft2(Img)                                     # Perform a 2D-fast fourier transform on the image
    mask_fft = np.fft.fft2(mask)                                   # Perform a 2D-fast fourier transform on the mask
    
    Img_conv_fft = Img_fft*mask_fft;                                # Calculate the FT convolution of the image by multiplying the FT image and FT mask
    Img_conv = np.fft.ifft2(Img_conv_fft).real                     # Obtain the convoluted image by performing an 2D-inverse FFT
    Img_conv_norm = Img_conv/np.max(Img_conv)                      # Normalize the image such that the maximum value will always be 1 (background removal in pretreatment makes sure that the lowest value is always 0)
    Img_conv_norm = Img_conv_norm[r:-r,r:-r] 
    Img_conv_bin = Img_conv_norm                                   # Binarize the image, such that it can be used in regionprops
    Img_conv_bin[Img_conv_bin < threshold] = 0
    Img_conv_bin[Img_conv_bin >= threshold] = 1
    
    Pparticles = ParticleSelection(Img_conv_bin, Img_conv_norm, r, SelectionCriteria) # Determine if the found patches are particles according to the SelectionCriteria
    return Pparticles



#%%
def ParticleSelection(Img_conv_bin, Img_conv_norm, r, SelectionCriteria):
    Ncriteria = len(SelectionCriteria)                                     # Determine the number of criteria given
    Img_label = measure.label(Img_conv_bin)
    Props = measure.regionprops(Img_label, intensity_image=Img_conv_norm)                                      # Determine the centroids of all potential particles found
    Nregions = np.shape(Props)[0]                                               # Determine the number of potential particles found
    Centroids = np.array([Prop.weighted_centroid for Prop in Props])

    if Ncriteria == 0:
        print('No selection criteria were given, so all centroids are returned.\n')
        Pparticles = Centroids
        return Pparticles
    elif Nregions == 0:
        print('No particles were found.\n')
        Pparticles = []
        return

    Selection = np.zeros((Ncriteria,Nregions)).astype(bool)
    for i in range(Ncriteria):
        Property = SelectionCriteria.Property[i]                          # Region property to check
        Value = SelectionCriteria.Value[i]                             # Value to check the property against
        Criteria = SelectionCriteria.Criteria[i]                          # Criteria to be a valid particle (either 'greater' or 'smaller')
        P_properties = np.array([Prop[Property] for Prop in Props])                              # Determine the values of the regionproperties to compare with Value
        # Adds a vector with logical values to Selection, where 1 means the
        # regionproperty of that region complies with the given criteria
        # and 0 means that the regionproperty doesnt comply with the
        # criteria.
        if Criteria.lower() == 'greater':
            Selection[i,:] = np.array([P_properties > Value])
        elif Criteria.lower() == 'smaller':
            Selection[i,:] = np.array([P_properties < Value])

    # AND gate with Ncriteria inputs, output == 1 if all inputs are 1, else output == 0
    Pparticles = Centroids[np.prod(Selection,axis=0).astype(bool)][:]

    if np.shape(Pparticles)[1] == 0:
        print('No particles met the requirements');
        return
    else:
        return Pparticles
    
    
    
    
#%% Debugging
if __name__ == '__main__':
    pass