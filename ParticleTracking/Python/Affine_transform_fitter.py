# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:40:47 2020

@author: Lars Kool
"""

import numpy as np
import pandas as pd
import math as m
import matplotlib.pyplot as plt

# Euclidian distance between points [x1[i],y1[i]] and [x2[i],y2[i]]
def distance(x1,x2,y1,y2):
    distance = np.sqrt(np.square(x1 - x2) + np.square(y1 - y2))
    return distance

def d2min(p0, p1, MaxD):
    result = np.zeros((len(p0),1)) # Initiate the output matrix
    test = np.empty([len(p0.index)])
    
    for i in p0.index: # So, do this loop for every input particle
        # Determine the neighbors for this particle
        distances = distance(p0.x[:],p0.x[i],p0.y[:],p0.y[i])
        neighbors = np.array(distances < MaxD)
        neighbors[i] = False # As Rii = 0, the particle itself is also counted as a neighbor, and should be removed from the list
        nNeighbors = sum(neighbors)
            
        if nNeighbors > 0:
            ref = np.array([p0.x[i], p0.y[i], 1.]) # P_i(t)
            ref_next = np.array([p1.x[i], p1.y[i], 1.]) # P_i(t+dt)
            
            before = np.array([np.array(p0.x[neighbors]), np.array(p0.y[neighbors]), np.ones([nNeighbors,])]) # P_ij(t)
            next_frame = np.array([np.array(p1.x[neighbors]), np.array(p1.y[neighbors]), np.ones([nNeighbors,])]) # P_ij(t+dt)
        
            M =  np.linalg.lstsq(np.array(np.transpose(before)),  np.transpose(np.array(next_frame)), rcond=None) # Best fit affine transformation such that Pij_(t+dt)-M*P_ij(t) is minimalized
            
            after = np.dot(np.transpose(M[0]), before) # M*P_ij(t)
            ref_after = np.dot(np.transpose(M[0]), ref) # M*P_i(t)
            
            Rij_next = np.linalg.norm(np.transpose(next_frame) - ref_next, ord=2, axis=1) # R_ij(t+dt)
            Rij_after = np.linalg.norm(np.transpose(after) - ref_after, ord=2, axis=1) # M*R_ij(t)
            
            # Calculation of D2min  
            result[i] = np.sum(np.square(Rij_next - Rij_after)) / (nNeighbors)
        else:
            result[i] = 0
    return result


if __name__ == '__main__':
    tracked_complete = pd.read_pickle('F:/Lars/Oscillatory Compression/20200713 Soft Particle/Avg100_Amp80_Per120_Back25/Preprocessed/Complete_tracks_large.pkl')
    MaxD = 200
    result = []
    p0 = tracked_complete[tracked_complete.frame == 106].reset_index(drop=True)
    p1 = tracked_complete[tracked_complete.frame == 106+120].reset_index(drop=True)
    
    for i in range(2):
        p0 = tracked_complete[tracked_complete.frame == 106+120*i].reset_index(drop=True)
        p1 = tracked_complete[tracked_complete.frame == 106+120*(i+1)].reset_index(drop=True)
    
        result.append(d2min(p0,p1, MaxD))
    
    # test = np.empty([9,2])
    # test[:,0] = np.mean(result,axis=1).reshape(9,)
    # test[:,1] = np.std(result,axis=1).reshape(9,)

    px_size = 5000/1997 #(um/particle)
    p = 10
    cycle = 8
    start = 46
    period = 120
    
    particle = tracked_complete[tracked_complete.particle == p].reset_index(drop=True)
    
    plt.figure(dpi=500)
    plt.plot(particle.x[start+cycle*period:start+(cycle+1)*period]*px_size,particle.y[start+cycle*period:start+(cycle+1)*period]*px_size)
    plt.plot(particle.x[start+(cycle+1)*period:start+(cycle+2)*period]*px_size,particle.y[start+(cycle+1)*period:start+(cycle+2)*period]*px_size)
    plt.xlim((3000,3450))
    plt.ylim((248,260))
    plt.xlabel('x-position $(\mu m)$')
    plt.ylabel('y-position $(\mu m)$')
    
    plt.figure(dpi=500)
    plt.plot(particle.x*px_size)
    plt.xlabel('Time (sec)')
    plt.ylabel('x-position $(\mu m)$')





































        
        