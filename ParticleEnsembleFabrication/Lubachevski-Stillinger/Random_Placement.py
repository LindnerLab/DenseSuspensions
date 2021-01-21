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
import pandas as pd
from distance_to_neighbors import distance_to_neighbors
from neighboring_particles import neighboring_particles

class particles:
    def __init__(self, xmax, ymax, R, N):
        self.pos = np.array([np.random.rand(np.sum(N))*xmax,
                             np.random.rand(np.sum(N))*ymax]).transpose()
        self.r = np.concatenate([np.ones([Nparticles[i], ])*Rparticles[i] for i in range(len(Nparticles))], axis=0)
        self.v = np.random.rand(np.sum(N),2)*2-1
        self.m = np.ones([np.sum(N),])
        self.neighbors = []
        self.forces = [[] for i in range(np.sum(N))]
        self.Fsum = np.zeros([np.sum(N),2])
        self.a = np.zeros([np.sum(N),2])
        self.a_dt = np.zeros([np.sum(N),2])
        self.Ntot = np.sum(N)
        
    # Generate a new list of neighbors for each particle. This list does not
    # have to be updated every cycle (frequency depends on the speed and
    # average interparticle distance)
    def update_neighbors(self, cutoff):
        self.neighbors = neighboring_particles(self.pos, cutoff)

    def update_position(self, dt):
        self.pos = np.array([self.pos[:,0] + self.v[:,0]*dt + 0.5*self.a[:,0]*(dt**2),
                             self.pos[:,1] + self.v[:,1]*dt + 0.5*self.a[:,1]*(dt**2)]).transpose()
        
    def update_force(self):
        for i in range(self.Ntot):
            forces = []
            xy_particle = self.pos[i,:]
            xy_neighbors = self.pos[self.neighbors[i]]
            d_neighbors = distance_to_neighbors(xy_particle, xy_neighbors)
            
            for j in range(len(self.neighbors[i])):
                if d_neighbors > (self.r[i] + self.r[self.neighbor[i][j]]):
                    forces.append()

    def update_acceleration(self):
        self.a = self.a_dt
        self.a_dt = self.force/self.m
    
    def update_velocity(self):
        self.v = self.v + 0.5*(self.a + self.a_dt)*dt
    


if __name__ == '__main__':
    Nparticles = [1500, 1000]
    Rparticles = [0.03, 0.045]
    length = 25
    width = 5
    cutoff = 0.3
    dt = 0.001
    
    packing = particles(length, width, Rparticles, Nparticles)
    packing.update_neighbors(cutoff)
    packing.update_force()