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

from track_particles import track_particles

datapath = r'G:\Lars\Oscillatory Compression\20210222\Avg125_Amp100_Back25_Per600_C50\Preprocessed\V1'
dataFormat = 'pkl'
minLength = 1100
search_range = 35

track_particles(datapath, dataFormat, minLength, search_range)