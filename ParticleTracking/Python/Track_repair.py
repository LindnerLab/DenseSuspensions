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

@author: Lars Kool
"""
def trackRepair(tracked, maxMissing, nFiles, minLength):
    import pandas as pd
    import numpy as np
    from Neighboring_elements import neighboringElements
    
    # Determine the length of each found track
    track_lengths = tracked.particle.value_counts()
    track_lengths = track_lengths.sort_index()
    nParticles = len(track_lengths)
    
    # Remove all short tracks
    tracked_temp = pd.DataFrame(columns = ['y','x','r','frame','particle'])
    for i in range(nParticles):
        if track_lengths[i] >= minLength:
            tracked_temp = pd.concat([tracked_temp, tracked[tracked.particle == i]])
        doneness = (i/nParticles)*100
        if i % 50 == 0:
            print(f'{doneness:{3}.{3}} of the particles checked for completeness' )
    tracked_temp.reset_index(drop=True)
    nParticles = len(tracked.particle.value_counts())
    
    # Determine the length of each long track found
    track_lengths = tracked_temp.particle.value_counts()
    track_lengths = track_lengths.sort_index()
    nParticles = len(track_lengths)
    
    # Determine the unique particle ID's
    P_ids = tracked_temp.particle.unique()
    ints = pd.DataFrame(columns=['ints'])
    ints.ints = list(range(nFiles))
    track_complete = pd.DataFrame({'y':[],
                                   'x':[],
                                   'r':[],
                                   'frame':[],
                                   'particle':[]})
    
    for i in range(nParticles):
        # Slice all particles with Particle ID = P_ids[i]
        idx = P_ids[i]
        track_incomplete = tracked_temp[tracked_temp.particle == idx]
        track_incomplete = track_incomplete.reset_index(drop=True)
        track_incomplete.particle = i
        r = track_incomplete.r[0]
        
        if len(track_incomplete) < nFiles: # Check if all particles are present, if not con
            if list(track_incomplete.frame[-1:]) != [nFiles-1]: # If last particle is missing, copy the last particle found and add it at the end (last particles will be omitted from analysis anyway)
                new_particles = pd.DataFrame({'y':list(track_incomplete.y[-1:]),
                                              'x':list(track_incomplete.x[-1:]),
                                              'r':r,
                                              'frame':nFiles-1,
                                              'particle':i})
                track_incomplete = track_incomplete.append(new_particles,ignore_index = True)
            elif track_incomplete.frame[0] != 0: # If first particle is missing, copy the first particle found and add it. It will be omitted in analysis as well).
                new_particles = pd.DataFrame({'y':[track_incomplete.y[0]],
                                              'x':[track_incomplete.x[0]],
                                              'r':r,
                                              'frame':0,
                                              'particle':i})
                track_incomplete = track_incomplete.append(new_particles,ignore_index = True)
                track_incomplete = track_incomplete.sort_values(by=['frame'], ignore_index=True)
            
            # Determine in which frames the particle is not found, group neighboring frames and determine the length of the sections missing).
            missing = ~ints.ints.isin(track_incomplete.frame)
            missing = missing[missing == True].index
            missing_grouped = neighboringElements(missing)
            nMissing = []
            
            for j in range(len(missing_grouped)):
                nMissing.append(len(missing_grouped[j]))
                
            if max(nMissing) < maxMissing + 1: # Check if any of the missing sections is too big to repair, if not go and repair!
                # Fill in a missing section using linear interpolation
                for j in range(len(missing_grouped)):
                    nMissing = len(missing_grouped[j])
                    idx_start = track_incomplete.index[track_incomplete.frame == missing_grouped[j][0] -1]
                    idx_end = track_incomplete.index[track_incomplete.frame == missing_grouped[j][-1] +1]
                                        
                    x_start = track_incomplete.x[idx_start]
                    x_end = track_incomplete.x[idx_end]
                    y_start = track_incomplete.y[idx_start]
                    y_end = track_incomplete.y[idx_end]
                                   
                    new_x = np.linspace(x_start,x_end,nMissing+2)
                    new_x = np.transpose(new_x)
                    new_x = new_x[0][1:-1]
                    new_y = np.linspace(y_start,y_end,nMissing+2)
                    new_y = np.transpose(new_y)
                    new_y = new_y[0][1:-1]
                    
                    new_particles = pd.DataFrame({'y':new_y,
                                                 'x':new_x,
                                                 'r':[r]*nMissing,
                                                 'frame':missing_grouped[j],
                                                 'particle':[i]*nMissing})
                    track_incomplete = track_incomplete.append(new_particles,ignore_index = True)
                    track_incomplete = track_incomplete.sort_values(by=['frame'], ignore_index=True)
                # Once the incomplete track is repaired, add it to the completed tracks    
                track_complete = track_complete.append(track_incomplete)
            
        else: # If the track was already complete, add it to the completed tracks as well
            track_complete = track_complete.append(track_incomplete)
            
        if i % 100 == 0:
            print(f'{doneness:{3}.{3}} of the tracks repaired' )  
    return track_complete