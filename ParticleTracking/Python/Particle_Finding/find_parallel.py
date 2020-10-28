# -*- coding: utf-8 -*-
"""
Parallized version of find_serial. It is parallelized using the ray
library https://github.com/ray-project/ray

Parameters
----------
path : List of strings
    List of foldernames:
        Path[0]: General folder of the data to be analyzed \n
        Path[1]: Name of the folder containing the input images \n
        Path[2]: Name of the output folder \n
        Path[3]: Name of the version of the dataset (allows to distinguish
                 multiple datasets taken on the same day with identical
                 settings)
fname : String
    Name of the file (including extension) to be preprocessed that can be
    found in folder: path[0]/path[1]/path[3]
settings : dict
    Dict with settings to be used in the data analysis. See main particle
    finding script for a detailed explanation of all keys.

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
import ray
from image_pretreatment import image_pretreatment
from find_serial import find_serial

@ray.remote
def find_parallel(path, fname, settings):
    # Set verbose to False (parallel doesn' work with verbose), pretreat image
    # and then apply the 'normal' find_serial to the pretreated image
    settings['verbose'] = False
    img = image_pretreatment(path, fname, settings)
    find_serial(img, path, fname, settings)