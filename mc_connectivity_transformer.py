# coding=utf-8

import numpy as np
from itertools import product


def find_K(Cold, Npre, Npost):
    """
    compute indegree * Npost
    """
    h1 = np.log(1. - Cold)
    h2 = (Npre * Npost - 1.) / (Npre * Npost)
    h3 = h1 / np.log(h2) / Npost
    return h3


def find_new_C(K, Npre, Npost):
    """
    find new connectivity
    """
    h1 = (Npre * Npost - 1.) / (Npre * Npost)
    h2 = np.power(h1, K * Npost)
    h3 = 1. - h2
    return h3


def co(index):
    """
    for given index of old matrix compute new one
    """
    index2 = np.floor(index[ 0 ] / 3) * 2 + (index[ 0 ] % 3)
    index3 = np.floor(index[ 1 ] / 3) * 2 + (index[ 1 ] % 3)
    return np.array([ index2, index3 ])


def find_connectivity_from_connmat(connmat, popvec):
    """
    :param connmat: the OLD connectivity matrix
    :param popvec: the NEW population vector
    :return:
    """

    numpops = popvec.shape[ 0 ]
    newconnmat = np.zeros((numpops, numpops))
    popindex = np.arange(0, numpops)

    for i in product(popindex, popindex):
        if popvec[ i[ 0 ] ] == 0 or popvec[ i[ 1 ] == 0 ]:
            newconnmat[ i ] = 0.0
        else:
            K = find_K(connmat[ co(i) ], popvec[ i[ 0 ] ], popvec[ i[ 1 ] ])
            C = find_new_C(K, popvec[ i[ 0 ] ], popvec[ i[ 1 ] ])
            newconnmat[ i ] = C
    return newconnmat

if __name__ == "__main__":
    #             2/3e      2/3i    4e      4i      5e      5i      6e      6i
    conn_probs = [ [ 0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0., 0.0076, 0. ],
                   [ 0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0., 0.0042, 0. ],
                   [ 0.0077, 0.0059, 0.0497, 0.135, 0.0067, 0.0003, 0.0453, 0. ],
                   [ 0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0., 0.1057, 0. ],
                   [ 0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0. ],
                   [ 0.0548, 0.0269, 0.0257, 0.0022, 0.06, 0.3158, 0.0086, 0. ],
                   [ 0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252 ],
                   [ 0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443 ] ]

    N_full = {
        'L23': {'E': 20683, 'I': 5834},
        'L4':  {'E': 21915, 'I': 5479},
        'L5':  {'E': 4850, 'I': 1065},
        'L6':  {'E': 14395, 'I': 2948}
    }
    N_full = [  20683,
                5834,
                0,
                21915,
                5479,
                0,
                2425,
                1065,
                2425,
                14395,
                2948 ]

    newCa = find_connectivity_from_connmat(conn_probs, N_full)
    print(newCa)