# -*- coding: utf-8 -*-
"""
Author Information
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
import numpy as np
import pandas as pd
import trackpy as tp
from data_parser import data_parser
from Track_repair import track_repair
from renumber_particles import renumber_particles


def track_particles(datapath, dataFormat, minLength, search_range):
    types = os.listdir(os.path.join(datapath, dataFormat))
    # Import and parse data
    untracked = [data_parser(os.path.join(datapath, dataFormat, type_), dataFormat) for i, type_ in enumerate(types)]
    print('Particles imported')
    # Track the particles
    tracked_particle = [tp.link_df(particles, search_range, memory=5).sort_values(by=['particle', 'frame'], ignore_index=True)
                        for particles in untracked]
    nFiles = len(tracked_particle[0].frame.unique())
    del untracked
    # Remove small tracks and repair small holes (linear interpolation) and fill large holes with NaN
    tracked_repaired = [track_repair(track, 10, nFiles, minLength)
                        for track in tracked_particle]
    del tracked_particle
    # Save sorted, renumbered, and complete tracks
    tracked_complete = [[] for track in tracked_repaired]
    nParticles = [len(track.particle.unique()) for track in tracked_repaired]
    for i, track in enumerate(tracked_repaired):
        tracked_complete[i] = renumber_particles(track, 'particle', int(np.sum(nParticles[:i])))
    tracked_complete = pd.concat(tracked_repaired, ignore_index=True)
    tracked_complete = tracked_complete.sort_values(by=['frame', 'particle'], ignore_index=True)
    tracked_complete.to_pickle(os.path.join(datapath, r'tracked.pkl'))
