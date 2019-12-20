#! /usr/bin/python


import numpy as np
import networkx as nx


def get_p0 (G, selection):

    """
    Returns list of network nodes with starting probabilities considering
    nodes in selection (each node in selection with same probability)
    selection: list of start nodes
    """

    nodes = list(G.nodes())

    p0 = np.empty(len(nodes))
    for i, item in enumerate(nodes):
        if item in selection:
            p0[i] = 1./len(selection)
        else:
            p0[i] = 0
    return (p0)


def get_transition_matrix(G):

    adj = nx.adjacency_matrix(G)
    adj = adj.todense()
    ksum = adj.sum(axis=1)
    D = np.identity(adj.shape[0])
    np.fill_diagonal(D, ksum)
    T = np.dot(np.linalg.inv(D),adj)

    return (T)



def run_walk (p0, T, r = 0.2, conv_threshold = 0.000001):

    """
    Run random walk

    p0: starting probabilities
    T: transition matrix
    r: restart probability
    """

    p_t = np.copy(p0)
    diff_norm = 1
    while diff_norm > conv_threshold:
        p_t_1 = ((1-r) * (p_t * T)) + r * p0
        p_t_1 = np.asarray(p_t_1)
        diff_norm = np.linalg.norm(np.subtract(p_t_1,p_t),1)
        p_t = p_t_1
        p_t = np.squeeze(np.asarray(p_t))
    return (p_t)




if __name__ == '__main__':

    seeds = [query]
    p0 = get_p0(G,seeds)
    T = get_transition_matrix(G)
    rw = run_walk(p0, T, r = 0.1, n_steps = None, conv_threshold = 0.000001)
