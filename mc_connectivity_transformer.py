# coding=utf-8

import numpy as np
from itertools import product

def find_K(Cold, Npre, Npost):
    h1 = np.log(1. - Cold)
    h2 = (Npre * Npost - 1.) / (Npre * Npost)
    h3 = h1 / np.log(h2) / Npost
    return h3

def find_new_C(K, Npre, Npost):
    h1 = (Npre * Npost - 1.) / (Npre * Npost)
    h2 = np.power(h1, K * Npost)
    h3 = 1. - h2
    return h3

def co(index):
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



    # if newpopsize == 0.0:
    #     newconnmat[ :, 2 ] = 0.0
    #     newconnmat[ 2, 0 ] = 0.0
    #     newconnmat[ :2, :2 ] = connmat
    #     return newconnmat
    # else:
    #     factor = popvec[ 0, 0 ] / newpopsize
    #
    #     newconnmat[ 1, 1 ] = connmat[ 1, 1 ]  # I -> I
    #     K = find_K(connman[ 0, 0 ], popvec[ 0 ], popvec[ 0 ])
    #     newconnmat[ 0, 0 ] = find_new_C(K, )

