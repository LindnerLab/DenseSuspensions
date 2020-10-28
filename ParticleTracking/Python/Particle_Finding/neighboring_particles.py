# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:49:44 2020

@author: Charly
"""
import numpy as np
from distance_to_neighbors import distance_to_neighbors

def neighboring_particles(xy_particles, cutoff):
    """
    Determines which particles are within the cutoff distance of each particle.

    Parameters
    ----------
    xy_particles : n-by-2 numpy array of floats
        Coordinatees of n particles, where the column indicates x/y
        (col=0 -> x, col=1 -> y) and the row indicates the particle index.
    cutoff : float
        Maximum distance for which two particles can be considered neighbors.

    Returns
    -------
    neighbors : list of lists of ints
        Lists of lists of indices of neighoring particles, where
        len(neighbors[i]) will give the number of neighbors of particle i, and
        xy_particles[neighbors[i],:] will give a list of xy-coordinates of the
        particles that neighbor particle i.

    """
    nParticles = xy_particles.shape[0]
    neighbors = [[] for i in range(nParticles)]
    for i in range(nParticles):
        distances = distance_to_neighbors(xy_particles[i, :],
                                          xy_particles[i+1:, :])
        neighbors_particle = distances < cutoff
        idx = [j+i+1 for j, x in enumerate(neighbors_particle) if x == True]
        neighbors[i] = neighbors[i]+idx
        for j in idx:
            neighbors[j] = neighbors[j]+[i]
    return neighbors