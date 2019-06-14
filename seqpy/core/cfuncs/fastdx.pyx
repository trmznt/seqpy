
from seqpy import cerr
import numpy as np
import cython

def pwdist(haplotypes):
    """ calculate pairwise differences between haplotypes
    """

    cdef int n = len(haplotypes)
    cdef int m = len(haplotypes[0])
    cdef int i, j
    cdef float d, c
    cdef char *x_i
    cdef char *x_j
    distm = np.zeros( (n,n) )

    for i in range(n):
        cerr('I: pairwising %d' % i)
        for j in range(i+1, n):
            x_i = haplotypes[i]
            x_j = haplotypes[j]
            d = c = 0
            for idx in range(m):
                if x_i[idx] == b'-' or x_j[idx] == b'-':
                    continue
                if x_i[idx] != x_j[idx]:
                    d += 1
                c+= 1
            distm[i,j] = distm[j,i] = d/c

    return distm

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def pwdistm(haplotype_m):
    """ calculate pairwise differences between haplotypes in numpy matrix format
    """

    cdef int n = len(haplotype_m)
    cdef int m = len(haplotype_m[0])
    cdef int i, j
    cdef float d, c
    cdef short x_i
    cdef short x_j
    cdef signed char[:,:] M = haplotype_m
    distm = np.zeros( (n,n) )

    for i in range(n):
        cerr('I: pairwising %d' % i)
        for j in range(i+1, n):
            d = c = 0
            for idx in range(m):
                if M[i, idx] != M[j, idx]:
                    d += 1
                c+= 1
            distm[i,j] = distm[j,i] = d/c

    return distm


def dxy(distm, groups, group_keys):
    """ calculate average pairwise differences inter and intra populations
    """

    cdef int n = len(group_keys)
    cdef float d = 0
    cdef int c = 0
    dxym = np.zeros( (n,n) )

    for i in range(n):
        p_i = groups[group_keys[i]]
        for j in range(i, j):
            p_j = groups[group_keys[j]]
            d = 0.0
            c = 0
            for idx_i in p_i:
                for idx_j in p_j:
                    d += distm[idx_i, idx_j]
                    c += 1
            dxym[i,j] = dxym[j,i] = d/c

    return dxym


def da(dxym):
    """ calculate net differences between populations
    """

    cdef int n = len(dxym)
    da = np.zeros( (n,n) )

    for i in range(n):
        for j in range(i+1, n):
            da[i,j] = da[i,j] = dxym[i,j] - 0.5 * (dxym[i,i] + dxym[j,j])

    return da
