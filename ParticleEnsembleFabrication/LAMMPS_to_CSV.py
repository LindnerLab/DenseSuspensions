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

def data_parser(path, columns, start):
    with open(path, 'r') as lines:
        raw_input = [line for line in lines]
    start_idx = [i+1 for i, line in enumerate(raw_input) if start in line][0]
    
    raw_input = [line.replace('\n', '').split(' ') for line in raw_input[start_idx:]]
    raw_input = pd.DataFrame([list(map(float, row)) for row in raw_input],
                             columns=columns)
    types = raw_input.diameter.unique()
    parsed = [[] for _type in types]
    for i, _type in enumerate(types):
        parsed[i] = raw_input[raw_input.diameter == _type]
    return parsed

def Shortest_Path(data):
    nodes = np.array([data.x, data.y]).transpose()
    Nnodes = len(nodes)
    nodes_sorted = [[] for node in nodes]
    p = np.array([0, 0])
    
    for i in range(Nnodes):
        d = np.sqrt(np.sum((nodes - p) ** 2, 1))
        idx = np.where(d == np.min(d))
        p = nodes[idx][0]
        nodes_sorted[i] = p
        nodes = np.delete(nodes, idx, 0)
    return np.array(nodes_sorted)

def LAMMPS_to_CSV(path, columns, start):
    data_parsed = data_parser(path, columns, start)
    data_sorted = [[] for _type in data_parsed]
    for i, _type in enumerate(data_parsed):
        data_sorted[i] = Shortest_Path(_type)
        np.savetxt(path[:-4] + '_P' + str(i+1) + '.csv', data_sorted[i], delimiter=",")
    return data_parsed, data_sorted

def visualize_found_particles(particles, r, dpi=500, img=[]):
    patches = []
    fig, ax = plt.subplots(dpi=500)
    if img:
        ax.imshow(img, cmap='gray')
    color = ['red', 'blue']
    for i in range(len(particles)):
        for x, y in zip(particles[i][:, 0], particles[i][:, 1]):
            circle = mpl.patches.Circle((y, x), r[i],
                                        facecolor=color[i],
                                        alpha=0.4)
            patches.append(circle)
    p = mpl.collections.PatchCollection(patches, match_original=True)
    ax.add_collection(p)
    plt.axis('equal')
    plt.xlim([0, 100])
    plt.ylim([0, 50])
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    plt.show()

if __name__ == '__main__':
    path = r'E:\Lars\LAMMPS\20210518\V2\dump.25000.txt'
    columns = ['atom', 'group', 'x', 'y', 'z', 'diameter']
    start = r'ITEM: ATOMS'
    data_parsed, data_sorted = LAMMPS_to_CSV(path, columns, start)
    plt.figure(dpi=500)
    plt.plot(data_parsed[0].y, data_parsed[0].x,
              marker='o',
              markersize=2,
              linewidth=0.5)
    plt.axis('equal')
    
    plt.figure(dpi=500)
    plt.plot(data_sorted[0][:,1], data_sorted[0][:,0],
              marker='o',
              markersize=2,
              linewidth=0.5)
    plt.axis('equal')
    plt.plot(data_sorted[1][:,1], data_sorted[1][:,0],
              marker='o',
              markersize=2,
              linewidth=0.5)
    plt.axis('equal')
    
    visualize_found_particles(data_sorted, [0.6, 0.8])