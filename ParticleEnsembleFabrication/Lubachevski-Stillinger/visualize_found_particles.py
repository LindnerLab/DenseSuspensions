# -*- coding: utf-8 -*-
"""
Function to plot the found particles over the original image. This makes
it easy to see if the particles are identified correctly and none of the
particles are overlapping

Parameters
----------
particles : list of nParticles-by-2 numpy float array
    Locations of the centers of the particles. Where each array [i] in the
    list are particles of the same type (where r=radii[i]), and each row
    are the [x,y] coordinates of that particle..
img : y-by-x numpy float array
    Original image containing the particles to be found.
r : TYPE
    DESCRIPTION.

Returns
-------
None.

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
import matplotlib.pyplot as plt

def visualize_found_particles(particles, r, dpi=500, img=None):
    patches = []
    fig, ax = plt.subplots(dpi=dpi)
    if not img == None:
        ax.imshow(img, cmap='gray')
    color = ['red', 'blue']
    for i in range(len(particles)):
        for x, y in zip(particles[i][:, 0], particles[i][:, 1]):
            circle = mpl.patches.Circle((y, x), r[i],
                                        facecolor=color[i],
                                        alpha=0.4)
            patches.append(circle)
    p = mpl.collections.PatchCollection(patches, match_original=True)
    ax.add_collection(p)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.show()