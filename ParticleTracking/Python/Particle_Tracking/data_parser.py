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
import os
import pandas as pd


def data_parser(datapath, dataFormat):
    result = pd.DataFrame(columns=['y', 'x', 'r', 'frame'])

    if dataFormat == 'csv':
        files = [file for file in os.listdir(os.path.join(datapath)) if file.lower().endswith('.csv')]
        nFiles = len(files)

        result = pd.concat([pd.read_csv(datapath + file, names=['x', 'y', 'r', 'frame']) for file in files], ignore_index=True)

    elif dataFormat == 'pkl':
        files = [file for file in os.listdir(os.path.join(datapath)) if file.endswith('.pkl')]
        _temp = [[] for file in files]
        for idx, file in enumerate(files):
            new_particles = pd.read_pickle(os.path.join(datapath, file))
            new_particles['frame'] = idx
            _temp[idx] = new_particles
        result = pd.concat(_temp, ignore_index=True)
    else:
        print('No supported file format selected!')
    return result
