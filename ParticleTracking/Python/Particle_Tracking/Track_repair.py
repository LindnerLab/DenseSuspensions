# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 10:51:21 2020
This script repairs broken tracks (tracks where in some frames no particle was
found) by linearly interpolating the position of the particles using the last
known position and the next known position. It also removes all tracks that
which contain holes that are too big to repair (Gap > maxMissing)

Dependencies:
    - Pandas
    - Numpy
    - neighboringElements

Contributors : Lars Kool
Affiliations : Laboratoire Physique et Méchanique des Milieux Hétérogènes
(PMMH), ESPCI, Paris, France

This project has received funding from the European Union’s Horizon 2020
research and innovation programme under the Marie Skłodowska-Curie grant
agreement No. 813162
More info on this programme: https://caliper-itn.org/
"""
import pandas as pd
import numpy as np
from neighboring_elements import neighboring_elements


def track_repair(tracked, maxMissing, nFiles, minLength):
    # Determine the length of each found track
    track_lengths = tracked.particle.value_counts()
    track_lengths = track_lengths.sort_index()
    nParticles = len(track_lengths)

    # Remove all short tracks
    tracked_temp = pd.concat([tracked[tracked.particle == i] for i in np.unique(tracked.particle) if track_lengths[i] >= minLength])
    tracked_temp.reset_index(drop=True)
    nParticles = len(tracked.particle.value_counts())
    print('Short tracks removed')

    # Determine the length of each long track found
    track_lengths = tracked_temp.particle.value_counts()
    track_lengths = track_lengths.sort_index()
    nParticles = len(track_lengths)

    # Determine the unique particle ID's
    P_ids = tracked_temp.particle.unique()
    ints = pd.DataFrame(columns=['ints'])
    ints.ints = list(range(nFiles))
    track_complete = pd.DataFrame({'y': [],
                                   'x': [],
                                   'r': [],
                                   'frame': [],
                                   'particle': []
                                   })

    for i in range(nParticles):
        # Slice all particles with Particle ID = P_ids[i]
        idx = P_ids[i]
        track_incomplete = tracked_temp[tracked_temp.particle == idx]
        track_incomplete = track_incomplete.reset_index(drop=True)
        track_incomplete.particle = i
        r = track_incomplete.r[0]

        if len(track_incomplete) < nFiles:  # Check if all particles are present, if not continue
            if list(track_incomplete.frame[-1:]) != [nFiles - 1]:  # If last particle is missing, copy the last particle found and add it at the end (last particles will be omitted from analysis anyway)
                new_particles = pd.DataFrame({'y': list(track_incomplete.y[-1:]),
                                              'x': list(track_incomplete.x[-1:]),
                                              'r': r,
                                              'frame': nFiles - 1,
                                              'particle': i
                                              })
                track_incomplete = track_incomplete.append(new_particles, ignore_index=True)

            if track_incomplete.frame[0] != 0:  # If first particle is missing, copy the first particle found and add it. It will be omitted in analysis as well).
                new_particles = pd.DataFrame({'y': [track_incomplete.y[0]],
                                              'x': [track_incomplete.x[0]],
                                              'r': r,
                                              'frame': 0,
                                              'particle': i
                                              })
                track_incomplete = track_incomplete.append(new_particles, ignore_index=True)
            track_incomplete = track_incomplete.sort_values(by=['frame'], ignore_index=True)

            # Determine in which frames the particle is not found, group neighboring frames and determine the length of the sections missing).
            missing = ~ints.ints.isin(track_incomplete.frame)
            missing = missing[missing == True].index

            if len(missing) > 0:
                missing_grouped = neighboring_elements(missing)

                LengthMissing = [len(missing_grouped[j]) for j in range(len(missing_grouped))]

                for j in range(len(missing_grouped)):
                    nMissing = LengthMissing[j]
                    if nMissing < maxMissing + 1:  # Check if any of the missing sections is too big to repair, if not go and repair!
                        idx_start = track_incomplete.index[track_incomplete.frame == missing_grouped[j][0] - 1]
                        idx_end = track_incomplete.index[track_incomplete.frame == missing_grouped[j][-1] + 1]

                        x_start = track_incomplete.x[idx_start]
                        x_end = track_incomplete.x[idx_end]
                        y_start = track_incomplete.y[idx_start]
                        y_end = track_incomplete.y[idx_end]

                        new_x = np.linspace(x_start, x_end, nMissing + 2)
                        new_x = np.transpose(new_x)
                        new_x = new_x[0][1:-1]
                        new_y = np.linspace(y_start, y_end, nMissing + 2)
                        new_y = np.transpose(new_y)
                        new_y = new_y[0][1:-1]
                    else:
                        new_x = np.full([nMissing], np.nan)
                        new_y = np.full([nMissing], np.nan)

                    new_particles = pd.DataFrame({'y': new_y,
                                                  'x': new_x,
                                                  'r': np.full([nMissing], r),
                                                  'frame': np.array(missing_grouped[j]),
                                                  'particle': np.full([nMissing], i)
                                                  })
                    track_incomplete = track_incomplete.append(new_particles, ignore_index=True)
                    track_incomplete = track_incomplete.sort_values(by=['frame'], ignore_index=True)

        # Once the incomplete track is repaired or already complete, add it to the completed tracks
        track_complete = track_complete.append(track_incomplete)

        if i % 100 == 0:
            doneness = (i / nParticles) * 100
            print(f'{doneness:{3}.{3}} of the tracks repaired')
    return track_complete
