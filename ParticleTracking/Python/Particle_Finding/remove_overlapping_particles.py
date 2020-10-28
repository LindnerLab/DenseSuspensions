# -*- coding: utf-8 -*-
""" Remove particles that are overlapping. 

Parameters
----------
locs : list of nParticles-by-2 numpy float array
    Locations of the centers of the particles. Where each array [i] in the
    list are particles of the same type (where r=radii[i]), and each row
    are the [x,y] coordinates of that particle.
fit_vals : list of nParticles-by-2 numpy float array
    Normalized peak height of the particle center after convolution, where
    higher values are a better fit. Each array [i] in the list are
    particles of the same type (where r=radii[i]), and each row are the
    [regionprop.max_intensity, regionprops.mean_intensity] of that peak.
    When comparing particles, the particle with the highest fit_val value
    is taken as the 'correct' particle.
radii : list of floats
    List of particle radii, where each item in the list is a different
    particle type, such that len(locs)=len(fit_vals)=len(radii).

Returns
-------
locs : list of nParticles-by-2 numpy float array
    Locations of the centers of the particles. Where each array [i] in the
    list are particles of the same type (where r=radii[i]), and each row
    are the [x,y] coordinates of that particle. But now, none of the
    particles are overlapping. 
    
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
from neighboring_particles import neighboring_particles
from distance_to_neighbors import distance_to_neighbors

def remove_overlapping_particles(locs, fit_vals, radii):
    nTypes = len(locs)
    temp = [[],[]]
    temp1 = [[],[]]
    particles_removed = [[] for i in range(nTypes)]
    for i in range(nTypes):
        nParticles = len(locs[i])
        inter_type_overlap = neighboring_particles(locs[i], 1.5*radii[i])
        for j, overlaps in enumerate(inter_type_overlap):
            if overlaps:
                for overlap in overlaps:
                    if fit_vals[i][overlap][0] < fit_vals[i][j][0] and (overlap not in particles_removed[i]):
                        particles_removed[i].append(overlap)

        non_overlap = np.full((nParticles,), True)
        non_overlap[particles_removed[i]] = False
        locs[i] = locs[i][non_overlap]
        fit_vals[i] = fit_vals[i][non_overlap]
    
    particles_removed = [[] for i in range(nTypes)]
    nParticles = len(locs[0])
    for i in range(nParticles):
        overlaps = distance_to_neighbors(locs[0][i,:], locs[1]) < 0.7*np.sum(radii)
        if any(overlaps):
            fits = fit_vals[1][overlaps,0]
            idx = np.arange(len(locs[1]))[overlaps]
            for fit, idx in zip(fits, idx):
                if fit > fit_vals[0][i,0]:
                    particles_removed[0].append(i)
                else:
                    particles_removed[1].append(idx)
    
    for i in range(nTypes):
        nParticles = len(locs[i])
        non_overlap = np.full((nParticles,), True)
        non_overlap[particles_removed[i]] = False
        locs[i] = locs[i][non_overlap]
        fit_vals[i] = fit_vals[i][non_overlap]
                
    return locs