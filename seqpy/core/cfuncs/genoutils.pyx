
from seqpy import cerr, cexit
import numpy as np
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def ralt(genotypes):

    cdef int shape = len(genotypes)
    cdef int gt_tot = 0
    data = np.empty(shape = shape, dtype=np.double)
    cdef short[:, :] genotype_view = genotypes
    cdef double[:] data_view = data
    cdef short[:] gt_view

    for i in range(shape):
        gt_view = genotype_view[i]
        #gt = genotypes[i]
        #tot = gt[0] + gt[1]
        gt_tot = gt_view[0] + gt_view[1]
        if gt_tot == 0:
            data_view[i] = -1
        else:
            data_view[i] = <double> gt_view[0]/gt_tot                        
            #dis = (data_view[i] - gt[0]/tot) ** 2
            #if  dis > 1e-5: 
            #    cerr('Data: %d %d %f| %d %d %f'% (gt_view[0], gt_tot, <double> gt_view[0]/gt_tot, gt[0], tot, gt[0]/tot))

#                cerr('Disreparancy: %f' % dis)

    return data

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def ralt_to_nalt(array, threshold=-1):

    cdef int shape = len(array)
    cdef float i_threshold = threshold
    cdef double[:] array_v = array
    cdef float r
    nalt = np.empty(shape = shape, dtype=np.short)
    cdef short[:] nalt_v = nalt

    if i_threshold < 0:
        for i in range(shape):
            r = array_v[i]
            if r < 0:
                nalt_v[i] = -1
            elif r > 0.5:
                nalt_v[i] = 2
            else:
                nalt_v[i] = 0

    else:
        for i in range(shape):
            r = array_v[i]
            if r < 0:
                nalt_v[i] = -1
            elif r < i_threshold:
                nalt_v[i] = 0
            elif r > (1.0 - i_threshold):
                nalt_v[i] = 2
            else:
                nalt_v[i] = 1
        print(1-i_threshold)
        print(array)
        print(nalt)

    return nalt
