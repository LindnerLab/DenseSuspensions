# function that takes an N-by-2 numpy matrix, and sorts them to yield a short path to pass through all nodes
import numpy as np

def short_path(nodes):
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