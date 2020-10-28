# -*- coding: utf-8 -*-
"""Creation of an annulus mask.
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
import matplotlib as mpl

def create_annulus_mask(radii, img_size):
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