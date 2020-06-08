# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:40:47 2020

@author: Lars Kool
"""

import numpy as np
import pandas as pd
import math as m
from scipy.optimize import newton_krylov

def affinetrans(particle_id,frame1,frame2,matrix):
    pass

def distance(x1,x2,y1,y2):
    distance = np.sqrt(np.square(x1 - x2) + np.square(y1 - y2))
    return distance

def affine_transformation(x,y,M):
    result = np.zeros([len(x), 3])
    for i in range(len(x)):
        result[i] = np.matrix([x[i], y[i], 1])*M
    
    return result[:,0], result[:,1]

def d2min(p0, p1):
    result = np.zeros((len(p0),1))
    for i in p0.index:
        distances = distance(p0.x[:],p0.x[i],p0.y[:],p0.y[i])
        neighbors = np.array(distances < 200)
        
        ones = np.array([1 for i in range(sum(neighbors))])
        ref = np.array([p0.x[i], p0.y[i], 1.])
        ref_next = np.array([p1.x[i], p1.y[i], 1.])
        before = np.array([p0.x[neighbors], p0.y[neighbors], ones])
        next_frame = np.array([p1.x[neighbors], p1.y[neighbors], ones])
    
        M = np.linalg.lstsq(np.array(np.transpose(before)),  np.transpose(np.array(next_frame)), rcond=None)
        
        after = np.dot(np.transpose(M[0]), before)
        ref_after = np.dot(np.transpose(M[0]), ref)
        
        Rij_next = np.linalg.norm(np.transpose(next_frame) - ref_next, ord=2, axis=1)
        Rij_after = np.linalg.norm(np.transpose(after) - ref_after, ord=2, axis=1)
        
        result[i] = np.sum(np.square(Rij_next - Rij_after)) / (len(Rij_after) - 1)
    return result


if __name__ == '__main__':
    tracked_complete = pd.read_pickle('E:/Lars/Failure/20200310/v2/Processed/tracked_complete.pkl')
   
    #Change this to plain indexing, as number of particles is known
    p0 = pd.concat([tracked_complete.iloc[[i]] for i in tracked_complete.index if tracked_complete.frame[i] == 1])
    p0 = p0.reset_index(drop=True)
    p1 = pd.concat([tracked_complete.iloc[[i]] for i in tracked_complete.index if tracked_complete.frame[i] == 537])
    p1 = p1.sort_values('particle')
    p1 = p1.reset_index(drop=True)
    
    result = d2min(p0,p1)
    
    
    
    
    
    # for i in range(len(p0)):
    # i = 100
    # distances = distance(p0.x[:],p0.x[i],p0.y[:],p0.y[i])
    # neighbors = np.array(distances < 150)
    
    # ones = np.array([1 for i in range(sum(neighbors))])
    # ref = np.array([p0.x[i], p0.y[i], 1.])
    # ref_next = np.array([p1.x[i], p1.y[i], 1.])
    # before = np.array([p0.x[neighbors], p0.y[neighbors], ones])
    # next_frame = np.array([p1.x[neighbors], p1.y[neighbors], ones])

    # result = np.linalg.lstsq(np.array(np.transpose(before)),  np.transpose(np.array(next_frame)), rcond=None)
    
    # after = np.dot(np.transpose(result[0]), before)
    # ref_after = np.dot(np.transpose(result[0]), ref)
    
    # Rij_next = np.linalg.norm(np.transpose(next_frame) - ref_next, ord=2, axis=1)
    # Rij_after = np.linalg.norm(np.transpose(after) - ref_after, ord=2, axis=1)
    
    # test = np.sum(np.square(Rij_next - Rij_after)) / (len(Rij_after) - 1)
    
    
    
    
    
        
        
        
        
        
        
        