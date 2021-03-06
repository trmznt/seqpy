#cython: language_level=3

from seqpy import cerr, cexit
import numpy as np
import cython
from numpy cimport int8_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def ralt(genotypes, mindepth=1):

    cdef int shape = len(genotypes)
    cdef int gt_tot = 0
    data = np.empty(shape = shape, dtype=np.double)
    read = np.empty(shape = shape, dtype=np.short)
    depth = np.empty(shape = shape, dtype=np.ushort)
    cdef short[:, :] genotype_view = genotypes
    cdef double[:] data_view = data
    cdef short[:] read_view = read
    cdef unsigned short [:] depth_view = depth
    cdef short[:] gt_view
    cdef short gt_0
    cdef short gt_1
    cdef short min_depth = mindepth

    for i in range(shape):
        gt_view = genotype_view[i]
        #gt = genotypes[i]
        #tot = gt[0] + gt[1]
        gt_0 = gt_view[0]
        gt_1 = gt_view[1]
        gt_tot = gt_0 + gt_1
        if gt_tot < min_depth:
            data_view[i] = -1
        else:
            data_view[i] = <double> gt_0/gt_tot
            #dis = (data_view[i] - gt[0]/tot) ** 2
            #if  dis > 1e-5:
            #    cerr('Data: %d %d %f| %d %d %f'% (gt_view[0], gt_tot, <double> gt_view[0]/gt_tot, gt[0], tot, gt[0]/tot))

#                cerr('Disreparancy: %f' % dis)
        read_view[i] = gt_0 if gt_0 < gt_1 else gt_1
        depth_view[i] = gt_tot

    return data, read, depth

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def ralt_to_nalt(r_alt, n_mdp, threshold=-1, minmdp=3):

    cdef float i_threshold = threshold
    cdef double[:,:] ralt_v = r_alt
    cdef int[:,:] nmdp_v = n_mdp
    cdef float r
    cdef int d
    cdef int N = r_alt.shape[0]
    cdef int L = r_alt.shape[1]
    cdef int n
    cdef int i
    nalt = np.empty(shape = r_alt.shape, dtype=np.int8)
    cdef int8_t[:,:] nalt_v = nalt
    cdef short min_mdp = minmdp


    if i_threshold < 0:
        for n in range(N):
            for i in range(L):
                r = ralt_v[n,i]
                if r < 0:
                    nalt_v[n,i] = -1
                elif r > 0.5:
                    nalt_v[n,i] = 2
                else:
                    nalt_v[n,i] = 0

    elif min_mdp <= 0:
        for n in range(N):
            for i in range(L):
                r = ralt_v[n,i]
                if r < 0:
                    nalt_v[n,i] = -1
                elif i_threshold < r < (1.0 - i_threshold):
                    nalt_v[n,i] = 1
                elif r < 0.5:
                    nalt_v[n,i] = 0
                else:
                    nalt_v[n,i] = 2

    else:
        for n in range(N):
            for i in range(L):
                r = ralt_v[n,i]
                d = nmdp_v[n,i]
                if r < 0:
                    nalt_v[n,i] = -1
                elif i_threshold < r < (1.0 - i_threshold) and d > min_mdp:
                    nalt_v[n,i] = 1
                elif r < 0.5:
                    nalt_v[n,i] = 0
                else:
                    nalt_v[n,i] = 2

    return nalt
