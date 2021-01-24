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
import matplotlib as mpl
import matplotlib.pyplot as plt
from distance_to_neighbors import distance_to_neighbors
from neighboring_particles import neighboring_particles


class particles:
    def __init__(self, xmax, ymax, R, N, k_wall, k_particle, T):
        self.xmax = xmax
        self.ymax = ymax
        self.pos = np.array([np.random.rand(np.sum(N)) * xmax,
                             np.random.rand(np.sum(N)) * ymax]).transpose()
        self.r = np.concatenate([np.ones([Nparticles[i], ]) * Rparticles[i] for i in range(len(Nparticles))], axis=0)
        # self.v = np.random.rand(np.sum(N), 2) * 2 - 1
        self.v = initial_velocity(N)
        self.m = np.ones([np.sum(N), ])
        self.neighbors = []
        self.forces = [[] for i in range(np.sum(N))]
        self.Fsum = np.zeros([np.sum(N), 2])
        self.a = np.zeros([np.sum(N), 2])
        self.a_dt = np.zeros([np.sum(N), 2])
        self.Ntot = np.sum(N)
        self.k_wall = k_wall
        self.k_particle = k_particle
        
    def initial_velocity(self, N):
        phi = np.random.rand(np.sum(N)) * 2 * np.pi
        v = np.array([np.cos(phi), np.sin(phi)])
        pass

    # Generate a new list of neighbors for each particle. This list does not
    # have to be updated every cycle (frequency depends on the speed and
    # average interparticle distance)
    def update_neighbors(self, cutoff):
        self.neighbors = neighboring_particles(self.pos, cutoff)

    def update_position(self, dt):
        self.pos = np.array([self.pos[:, 0] + self.v[:, 0] * dt + 0.5 * self.a[:, 0] * (dt ** 2),
                             self.pos[:, 1] + self.v[:, 1] * dt + 0.5 * self.a[:, 1] * (dt ** 2)]).transpose()

    def update_force(self):
        for i in range(self.Ntot):
            forces = np.array([0, 0]).reshape([1, 2])
            xy_particle = self.pos[i, :]
            xy_neighbors = self.pos[self.neighbors[i]]
            d_neighbors = distance_to_neighbors(xy_particle, xy_neighbors)

            # Interparticle forces
            r1 = self.r[i]
            r2 = self.r[self.neighbors[i]]
            for j in range(len(self.neighbors[i])):
                if d_neighbors[j] < (self.r[i] + self.r[self.neighbors[i][j]]):
                    dpos = self.pos[i] - self.pos[self.neighbors[i][j]]
                    Fmag = (d_neighbors[j] - r1 - r2[j]) * self.k_particle
                    F = -np.array([(dpos[0] / d_neighbors[j]) * Fmag,
                                  (dpos[1] / d_neighbors[j]) * Fmag]).reshape([1, 2])
                    forces = np.concatenate([forces, F], axis=0)

            # Particle-wall forces
            if xy_particle[0] < 0:
                F = -np.array([xy_particle[0] * self.k_wall, 0]).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)
            elif xy_particle[0] > self.xmax:
                F = np.array([(self.xmax - xy_particle[0]) * self.k_wall, 0]).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)

            if xy_particle[1] < 0:
                F = -np.array([0, xy_particle[1] * self.k_wall]).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)
            elif xy_particle[1] > self.ymax:
                F = np.array([0, (self.ymax - xy_particle[1]) * self.k_wall]).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)

            self.forces[i] = forces
            self.Fsum[i, :] = np.sum(forces, axis=0)

    def update_acceleration(self):
        self.a = self.a_dt
        self.a_dt = self.Fsum / self.m[:, None]

    def update_velocity(self, dt):
        self.v = self.v + 0.5 * (self.a + self.a_dt) * dt


if __name__ == '__main__':
    Nparticles = [10]
    Rparticles = [0.03]
    Ncycles = 100000
    length = 1
    width = 1
    cutoff = 0.3
    dt = 0.0001
    k_particle = 100000
    k_wall = 10000000
    T = 10

    packing = particles(length, width, Rparticles, Nparticles, k_wall, k_particle, T)
    packing.update_neighbors(cutoff)

    trajectory = np.zeros((Ncycles, 10, 2))
    Ekin = np.array([])
    Epot = np.array([])
    Etot = np.array([])
    i = 0

    while i < Ncycles:
        Ekin = np.concatenate([Ekin, [0.5 * np.sum(packing.v ** 2)]])
        packing.update_position(dt)
        packing.update_force()
        packing.update_acceleration()
        packing.update_velocity(dt)
        trajectory[i, :, :] = packing.pos
        F = packing.forces
        Fsum = packing.Fsum
        v = packing.v
        pos = packing.pos
        i = i + 1

    plt.figure(dpi=500)
    ax = plt.subplot(111)
    for i in range(10):
        ax.plot(trajectory[:, i, 0], trajectory[:, i, 1])
    plt.xlim([0, 1])
    plt.ylim([0, 1])
