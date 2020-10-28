# -*- coding: utf-8 -*-
"""
This function checks if all the input and output folders are present.
And if not, it will create them.

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

import os
import sys

def check_file_structure(path, nTypes):
    if not os.path.isdir(os.path.join(path[0], path[2])):
        os.mkdir(os.path.join(path[0], path[2]))
        os.mkdir(os.path.join(path[0], path[2], path[3]))
        os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL'))
        for i in range(nTypes):
            os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL', 'P'+str(i+1)))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], path[2], path[3])):
        os.mkdir(os.path.join(path[0], path[2], path[3]))
        os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL'))
        for i in range(nTypes):
            os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL', 'P'+str(i+1)))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], path[2], path[3], 'PKL')):
        os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL'))
        for i in range(nTypes):
            os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL', 'P'+str(i+1)))
        print('Output File structure created')
    elif not os.path.isdir(os.path.join(path[0], path[2], path[3], 'PKL','P'+str(nTypes))):
        for i in range(nTypes):
            os.mkdir(os.path.join(path[0], path[2], path[3], 'PKL', 'P'+str(i+1)))       
    else:
        print('Output File structure already present')
        
    if not os.path.isdir(os.path.join(path[0], path[1], path[3])):
        print('''Input files not present in the expected folder:
              DIRECTORY\\INPUT\\VERSION)''')
        sys.exit()