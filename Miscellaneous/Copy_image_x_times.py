# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:17:40 2020

@author: Charly
"""

import PIL as PIL

path = 'F:/Lars/Oscillatory Compression/20200901 Fitting Noise/RawData/'
path_im = path+'original.tif'
# Load example image (to determine size)
Img = PIL.Image.open(path_im,mode='r')

for i in range(1000):
    Img.save(path+'{:04d}'.format(i)+'.tif')


