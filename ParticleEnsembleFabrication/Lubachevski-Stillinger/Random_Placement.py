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
import matplotlib as mpl
import matplotlib.pyplot as plt
from DEM_class import particles


if __name__ == '__main__':
    Nparticles = [1500, 1000]
    Rparticles = [0.06, 0.09]
    length = 30
    width = 4.98
    # Nparticles = [100]
    # Rparticles = [0.06]
    nCores = 12
    # length = 2
    # width = 2
    cutoff = 0.15
    dt = 0.0002
    k_particle = 10000
    k_wall = 10000000
    T = 10
    
    """
    No user input required below this point
    """   
    packing = particles(length, width, Rparticles, Nparticles, k_wall, k_particle, T)
    packing.initial_velocity()
    cutoff = 2.5*np.max([packing.v[:,0]**2 + packing.v[:,1]**2])*dt*100
    packing.update_neighbors(cutoff)
    packing.update_force()
    packing.update_energies()

    Ekin = np.array([])
    Epot = np.array([])
    Etot = np.array([])
    phi = np.array([])
    packings = []
    r = []

    i = 0
    print('Start equilibration')
    while packing.Epot > 0.1 or i > 1E5:
        packing.update_position(dt)
        packing.update_force()
        packing.update_acceleration()
        packing.update_velocity(dt)
        packing.update_energies()
        
        if i % 10 == 0:
            packing.update_neighbors(2 * np.max(packing.r))
            # packing.v = packing.v / np.sqrt(packing.Ekin / T)
            packing.initial_velocity()
            
        if i % 100 == 0:
            print(f'Epot = {packing.Epot:.3f}, equilibrated when Epot < 0.1')
        i += 1

    packing.initial_velocity()
    i = 0
    cont = True
    while cont:
        phi = np.concatenate([phi, [np.sum(packing.r**2 * np.pi) / (length * width)]])
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
            # packing.v = packing.v / np.sqrt(Ekin[-1] / T)
            packing.initial_velocity()

        if i % 2000 == 0:
            packings.append(packing.pos)
            r.append(packing.r)
            
        if phi[i] > 0.83:
            cont = False

        packing.r = packing.r + r[0]/30000
        i += 1

    for i in range(len(r)):
        patches = []
        fig, ax = plt.subplots(dpi=500)
        color = ['red']
        for x, y, radius in zip(packings[i][:, 0], packings[i][:, 1], r[i]):
            r_temp = np.mean(radius)
            circle = mpl.patches.Circle((x, y), radius,
                                        facecolor='r',
                                        alpha=0.4)
            patches.append(circle)
        # rectangle = mpl.patches.Rectangle([-r_temp, -r_temp], 1 + 2*r_temp, 1 + 2*r_temp, angle=0,
        #                                   linestyle='--',
        #                                   edgecolor='grey',
        #                                   facecolor='none')
        # patches.append(rectangle)
        p = mpl.collections.PatchCollection(patches, match_original=True)
        ax.add_collection(p)
        ax.set_aspect('equal', 'box')
        plt.xlim([0, length])
        plt.ylim([0, width])
        plt.text(0.79*length, 0.95*width, '$\phi$={:.3f}'.format(phi[2000*i+1]))
        plt.show()
