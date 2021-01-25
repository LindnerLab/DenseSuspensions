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
from pair_correlation_function import pairCorrelationFunction_2D
from visualize_found_particles import visualize_found_particles


class particles:
    def __init__(self, xmax, ymax, R, N, k_wall, k_particle):
        self.xmax = xmax
        self.ymax = ymax
        self.pos = np.array([np.random.rand(np.sum(N)) * xmax,
                             np.random.rand(np.sum(N)) * ymax]).transpose()
        self.r = np.concatenate([np.ones([Nparticles[i], ]) * Rparticles[i] for i in range(len(Nparticles))], axis=0)
        # self.v = np.random.rand(np.sum(N), 2) * 2 - 1
        self.v = np.zeros([np.sum(N), 2])
        self.m = np.ones([np.sum(N), ])
        self.neighbors = []
        self.forces = [[] for i in range(np.sum(N))]
        self.Fsum = np.zeros([np.sum(N), 2])
        self.a = np.zeros([np.sum(N), 2])
        self.a_dt = np.zeros([np.sum(N), 2])
        self.Ntot = np.sum(N)
        self.k_wall = k_wall
        self.k_particle = k_particle
        self.Epot = 0
        self.Ekin = 0
        self.Etot = 0

    def initial_velocity(self):
        phi = np.random.rand(self.Ntot) * 2 * np.pi
        self.v = np.array([np.cos(phi), np.sin(phi)]).transpose()

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
            forces = np.array([0, 0, 1]).reshape([1, 3])
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
                                   (dpos[1] / d_neighbors[j]) * Fmag,
                                   self.k_particle]).reshape([1, 3])
                    forces = np.concatenate([forces, F], axis=0)

            # Particle-wall forces
            if xy_particle[0] < 0:
                F = -np.array([xy_particle[0] * self.k_wall, 0,
                               self.k_wall]).reshape([1, 3])
                forces = np.concatenate([forces, F], axis=0)
            elif xy_particle[0] > self.xmax:
                F = np.array([(self.xmax - xy_particle[0]) * self.k_wall, 0,
                              self.k_wall]).reshape([1, 3])
                forces = np.concatenate([forces, F], axis=0)

            if xy_particle[1] < 0:
                F = -np.array([0, xy_particle[1] * self.k_wall,
                               self.k_wall]).reshape([1, 3])
                forces = np.concatenate([forces, F], axis=0)
            elif xy_particle[1] > self.ymax:
                F = np.array([0, (self.ymax - xy_particle[1]) * self.k_wall,
                              self.k_wall]).reshape([1, 3])
                forces = np.concatenate([forces, F], axis=0)

            self.forces[i] = forces
            self.Fsum[i, :] = np.sum(forces[:, 0:2], axis=0)

    def update_acceleration(self):
        self.a = self.a_dt
        self.a_dt = self.Fsum / self.m[:, None]

    def update_velocity(self, dt):
        self.v = self.v + 0.5 * (self.a + self.a_dt) * dt

    def update_energies(self):
        self.Ekin = 0.5 * np.sum(self.m[:, None] * self.v ** 2)

        Epot = np.zeros([self.Ntot, ])
        for i, force in enumerate(self.forces):
            dx = np.divide(force[:, 0:2], force[:, 2, None])
            Epot[i] = 0.25 * np.sum(abs(force[:, 2, None]) * (dx[:, 0:2] ** 2))
        self.Epot = np.sum(Epot)
        self.Etot = self.Ekin + self.Epot


if __name__ == '__main__':
    Nparticles = [20]
    Rparticles = [0.03]
    length = 1
    width = 1
    cutoff = 0.15
    dr = 0.01
    dt = 0.0002
    k_particle = 1000
    k_wall = 1000000
    T = 1
    
    """
    No user input required below this point
    """

    packing = particles(length, width, Rparticles, Nparticles, k_wall, k_particle)
    packing.initial_velocity()
    packing.update_neighbors(cutoff)

    Ekin = np.array([0])
    Epot = np.array([0])
    Etot = np.array([0])
    phi = np.array([0])
    packings = []
    r = []
    i = 0

    while packing.Epot < 1000:
        l = length + 2 * np.mean(packing.r)
        w = width + 2 * np.mean(packing.r)
        phi = np.concatenate([phi, [np.sum(packing.r**2 * np.pi) / (l * w)]])
        packing.update_position(dt)
        packing.update_force()
        packing.update_acceleration()
        packing.update_velocity(dt)
        packing.update_energies()
        Ekin = np.concatenate([Ekin, [packing.Ekin]])
        Epot = np.concatenate([Epot, [packing.Epot]])
        Etot = np.concatenate([Etot, [packing.Etot]])

        if i % 100 == 0:
            packing.update_neighbors(2 * np.max(packing.r))
            print(f'Ekin={Ekin[-1]:.2f}, Epot={Epot[-1]:.2f}, phi={phi[-1]:.2f}')
            packing.v = packing.v / np.sqrt(Ekin[-1] / T)

        if i % 20000 == 0:
            packing.r = packing.r + dr
            packings.append(packing.pos)
            r.append(packing.r)

        i += 1

    plt.figure(dpi=500)
    plt.plot(Ekin)

    i = 9
    patches = []
    fig, ax = plt.subplots(dpi=500)
    color = ['red']
    for x, y, radius in zip(packings[i][:, 0], packings[i][:, 1], r[i]):
        circle = mpl.patches.Circle((x, y), radius,
                                    facecolor='r',
                                    alpha=0.4)
        patches.append(circle)
    p = mpl.collections.PatchCollection(patches, match_original=True)
    ax.add_collection(p)
    ax.set_aspect('equal', 'box')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    plt.xlim([-np.mean(r[i]), 1 + np.mean(r[i])])
    plt.ylim([-np.mean(r[i]), 1 + np.mean(r[i])])
    plt.show()
