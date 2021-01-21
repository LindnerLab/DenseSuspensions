# -*- coding: utf-8 -*-
"""
Contributors : Lars Kool
Affiliations : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France
               
This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162
More info on this programme: https://caliper-itn.org/
"""
import numpy as np

def distance_to_neighbors(xy_particle, xy_neighbors):
    """
    Euclidian distance between a reference point (xy_particle [x,y]) and a
    list of n other points (xy_neighbors [[x1,y1],[x2,y2],...,[x_n, y_n]])

    Parameters
    ----------
    xy_particle : 1-by-2 array of floats
        Coordinates of the reference point [x,y].
    xy_neighbors : n-by-2 array of floats
        Coordinatees of n other points, where the column indicates x/y
        (col=0 -> x, col=1 -> y) and the row indicates the point.

    Returns
    -------
    distance : n-by-1 numpy array of floats
        Euclidian distance between the reference point and the other points,.

    """
    distance = np.sqrt(np.square(xy_neighbors[:, 0] - xy_particle[0])
                       + np.square(xy_neighbors[:, 1] - xy_particle[1]))
    return distance