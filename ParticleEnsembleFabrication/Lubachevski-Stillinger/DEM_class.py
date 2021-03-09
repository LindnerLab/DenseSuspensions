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
from distance_to_neighbors import distance_to_neighbors
from neighboring_particles import neighboring_particles
import matplotlib as mpl
import matplotlib.pyplot as plt


class particles:
    def __init__(self, xmax, ymax, R, N, k_wall, k_particle, T, damp):
        self.xmax = xmax
        self.ymax = ymax
        self.r = np.concatenate([np.ones([N[i], ]) * R[i] for i in range(len(N))], axis=0)
        self.pos = np.array([self.r + np.random.rand(np.sum(N)) * (xmax - 2 * self.r),
                             self.r + np.random.rand(np.sum(N)) * (ymax - 2 * self.r)]).transpose()
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
        self.T = T
        self.Fbox = np.array([0, 0]).astype(float)
        self.damp = damp * 2 * np.sqrt(0.5 * k_particle)

    def initial_velocity(self):
        v_mag = np.sqrt((2 * self.T) / self.Ntot)
        v_phi = np.random.rand(self.Ntot) * 2 * np.pi
        self.v = np.array([np.cos(v_phi), np.sin(v_phi)]).transpose() * v_mag

    # Generate a new list of neighbors for each particle. This list does not
    # have to be updated every cycle (frequency depends on the speed and
    # average interparticle distance)
    def update_neighbors(self, cutoff):
        self.neighbors = neighboring_particles(self.pos, cutoff)

    # Update the position of the particles using velocity-Verlet integration
    def update_position(self, dt):
        self.pos = np.array([self.pos[:, 0] + self.v[:, 0] * dt + 0.5 * self.a[:, 0] * (dt ** 2),
                             self.pos[:, 1] + self.v[:, 1] * dt + 0.5 * self.a[:, 1] * (dt ** 2)]).transpose()

    def update_force(self):
        self.forces = [[] for i in range(np.sum(self.Ntot))]
        self.Fbox = np.array([0, 0]).reshape([1, 2]).astype(float)
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
                    dpos = self.pos[self.neighbors[i][j]] - self.pos[i]
                    dv = self.v[self.neighbors[i][j]] - self.v[i]
                    theta = np.arctan2(dpos[1], dpos[0])
                    Fmag = (d_neighbors[j] - r1 - r2[j]) * self.k_particle
                    F = np.array([np.cos(theta) * Fmag + self.damp * dv[0],
                                  np.sin(theta) * Fmag + self.damp * dv[1]]
                                 ).reshape([1, 2])
                    forces = np.concatenate([forces, F], axis=0)

            # Particle-wall forces
            if xy_particle[0] < r1:
                F = -np.array([(xy_particle[0] - r1) * self.k_wall + self.damp * self.v[i, 0],
                               0]).reshape([1, 2])

                forces = np.concatenate([forces, F], axis=0)
                self.Fbox += F
            elif xy_particle[0] > self.xmax - r1:
                F = np.array([(self.xmax - r1 - xy_particle[0]) * self.k_wall - self.damp * self.v[i, 0],
                              0]).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)
                self.Fbox -= F

            if xy_particle[1] < r1:
                F = -np.array([0,
                               (xy_particle[1] - r1) * self.k_wall + self.damp * self.v[i, 1]]
                              ).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)
                self.Fbox += F
            elif xy_particle[1] > self.ymax - r1:
                F = np.array([0,
                              (self.ymax - r1 - xy_particle[1]) * self.k_wall] - self.damp * self.v[i, 1]
                             ).reshape([1, 2])
                forces = np.concatenate([forces, F], axis=0)
                self.Fbox -= F

            if not np.any(self.forces[i]):
                self.forces[i] = forces
            else:
                self.forces[i] = np.concatenate(self.forces[i], forces)

    def update_acceleration(self):
        self.Fsum = [np.sum(self.forces[i], axis=0) for i in range(len(self.forces))]
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
    Nparticles = [100]
    Rparticles = [0.05]
    nCores = 12
    length = 2
    width = 2
    cutoff = 0.15
    k_particle = 10000
    k_wall = 1000000
    dt = np.pi * np.sqrt(1 / k_wall) / 50
    T = 1
    damp = 0.7
    I_num = 1E-3

    """
    No user input required below this point
    """
    packing = particles(length, width, Rparticles, Nparticles, k_wall, k_particle, T, damp)
    packing.initial_velocity()
    cutoff = 5 * np.mean(packing.r)
    packing.update_neighbors(cutoff)
    packing.update_force()
    packing.update_acceleration()

    phi = np.array([])
    packings = []
    r = []

    i = 0
    print('Start equilibration')
    while np.sum(np.abs(packing.Fsum)) > 10:
        packing.update_position(dt)
        packing.update_force()
        packing.update_acceleration()
        packing.update_velocity(dt)

        if i % 50 == 0:
            packing.update_neighbors(5 * np.max(packing.r))

        if i % 1000 == 0:
            Epot = np.sum(np.abs(packing.Fsum))
            print(f'Epot = {Epot:.1f}, equilibrated when Epot < 10')

        i += 1

    P = 0
    i = 0
    packing.initial_velocity()
    while P < 1E3 and i < 1E6:
        packing.update_position(dt)
        packing.update_force()
        packing.update_acceleration()
        packing.update_velocity(dt)

        P = np.sum(packing.Fbox) / (2 * (packing.xmax + packing.ymax))

        # if P <= 1:
        #     dxmax = 1E-8 * packing.xmax
        # else:
        #     dxmax = I_num * dt * packing.xmax / np.sqrt(1 / P)
        dxmax = 1E-5

        packing.xmax -= dxmax

        if i % 100 == 0:
            packing.update_neighbors(5 * np.max(packing.r))

        if i % 1000 == 0:
            packings.append(packing.pos)
            r.append(packing.r)
            phi = np.sum(np.pi * packing.r ** 2) / (packing.xmax * packing.ymax)
            print(f'P={P:.1f}, phi={phi:.2f}')

            patches = []
            fig, ax = plt.subplots(dpi=500)
            color = ['red']
            for x, y, radius in zip(packing.pos[:, 0], packing.pos[:, 1], packing.r):
                r_temp = np.mean(radius)
                circle = mpl.patches.Circle((x, y), radius,
                                            facecolor='r',
                                            alpha=0.4)
                patches.append(circle)
            p = mpl.collections.PatchCollection(patches, match_original=True)
            ax.add_collection(p)
            ax.set_aspect('equal', 'box')
            plt.xlim([0, length])
            plt.ylim([0, width])
            # plt.text(0.79*length, 0.95*width, '$\phi$={:.3f}'.format(phi[2000*i+1]))
            plt.show()

        i += 1
