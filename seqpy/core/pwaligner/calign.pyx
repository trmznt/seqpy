#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as np
from libc.string cimport strlen
import os.path as op
import sys

cimport cython



cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject *PyByteArray_FromStringAndSize(char *, size_t)
    int PyByteArray_Resize(PyObject *, size_t)
    char * PyByteArray_AS_STRING(PyObject *)
    int PyByteArray_Check(PyObject *)
    int PyBytes_Check(PyObject *)
    char * PyBytes_AS_STRING(PyObject *)
#    PyObject *PyString_FromStringAndSize(char *, size_t)
#    int _PyString_Resize(PyObject **, size_t)
#    char * PyString_AS_STRING(PyObject *)

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.float32_t DTYPE_FLOAT

cdef inline DTYPE_FLOAT max3(DTYPE_FLOAT a, DTYPE_FLOAT b, DTYPE_FLOAT c):
    if c > b:
        return c if c > a else a
    return b if b > a else a

cdef inline DTYPE_FLOAT max2(DTYPE_FLOAT a, DTYPE_FLOAT b):
    return b if b > a else a

cdef object read_matrix(path, object cache={}):
    """
    so here, we read a matrix in the NCBI format and put
    it into a numpy array. so the score for a 'C' changing
    to an 'A' is stored in the matrix as:
        mat[ord('C'), ord('A')] = score
    as such, it's a direct array lookup from each pair in the alignment
    to a score. this makes it very fast. the cost is in terms of space.
    though it's usually less than 100*100.
    """
    if path in cache: return cache[path]


    cdef np.ndarray[DTYPE_INT, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v, mat_size
    if not op.exists(path):
        if "/" in path: raise Exception("path for matrix %s doest not exist" \
                                        % path)
        mat_path = op.abspath(op.join(op.dirname(__file__), "data"))
        fh = open(op.join(mat_path, path))
    else:
        fh = open(path)


    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [ord(x) for x in line.split(' ') if x]
    mat_size = max(headers) + 1

    a = np.zeros((mat_size, mat_size), dtype=np.int)

    line = fh.readline()
    while line:
        line_vals = [int(x) for x in line[:-1].split(' ')[1:] if x]
        for ohidx, val in zip(headers, line_vals):
            a[headers[ai], ohidx] = val
        ai += 1
        line = fh.readline()

    cache[path] = a
    return a

cdef char* as_string(object s, txt='expecting bytes or bytearray'):
    if PyByteArray_Check( <PyObject *> s):
        return PyByteArray_AS_STRING( <PyObject *> s)
    elif PyBytes_Check( <PyObject *> s):
        return PyBytes_AS_STRING( <PyObject *> s)
    raise RuntimeError( txt )




def max_index(array):
    """
    """
    max_value = array.argmax()
    idx = np.unravel_index(max_value, array.shape)
    return idx

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def aligner(object _seqj, object _seqi, \
            DTYPE_FLOAT gap_open=-7, DTYPE_FLOAT gap_extend=-7, DTYPE_FLOAT gap_double=-7,\
            method="global", matrix="BLOSUM62", fmt='alignment', debug=False):
    """
    Calculates the alignment of two sequences. The supported "methods" are
    "global" for a global Needleman-Wunsh algorithm, "local" for a local
    Smith-Waterman alignment, "global_cfe" for a global alignment with cost-free
    ends and "glocal" for an alignment which is "global" only with respect to
    the shorter sequence, this is also known as a "semi-global" alignment."
    Returns the aligned (sub)sequences as character arrays.

    Gotoh, O. (1982). J. Mol. Biol. 162, 705-708.
    Needleman, S. & Wunsch, C. (1970). J. Mol. Biol. 48(3), 443-53.
    Smith, T.F. & Waterman M.S. (1981). J. Mol. Biol. 147, 195-197.

    Arguments:

        - seqj (``sequence``) First aligned iterable object of symbols.
        - seqi (``sequence``) Second aligned iterable object of symbols.
        - method (``str``) Type of alignment: "global", "global_cfe", "local",
          "glocal".
        - gap_open (``float``) The gap-opening cost.
        - gap_extend (``float``) The cost of extending an open gap.
        - gap_double (``float``) The gap-opening cost if a gap is already open
          in the other sequence.
        - matrix (``dict``) A score matrix dictionary.
    """
    cdef int NONE = 0,  LEFT = 1, UP = 2,  DIAG = 3, p, last_p, counter
    cdef bint flip = 0

    cdef char* seqj = as_string( _seqj, 'expecting bytes or bytearray for first sequence')
    cdef char* seqi = as_string( _seqi, 'expecting bytes or bytearray for second sequence')

    cdef size_t align_counter = 0

    cdef int imethod = {"global": 0, "local": 1, "glocal": 2, "global_cfe": 3}[method]

    cdef size_t max_j = strlen(seqj)
    cdef size_t max_i = strlen(seqi)
    if max_i == max_j == 0:
        return "", "", 0

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i

    cdef char *align_j, *align_i
    cdef unsigned int i, j
    cdef char ci, cj
    cdef PyObject *ai, *aj
    cdef float max_score = 0

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    if debug: print("allocating memory...", file=sys.stderr)

    cdef np.ndarray[DTYPE_FLOAT, ndim=2] agap_i = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    cdef np.ndarray[DTYPE_FLOAT, ndim=2] agap_j = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    agap_i.fill(-np.inf)
    agap_j.fill(-np.inf)

    cdef np.ndarray[DTYPE_FLOAT, ndim=2] score = np.zeros((max_i + 1, max_j + 1), dtype=np.float32)

    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.zeros((max_i + 1, max_j + 1), dtype=np.uint)
    cdef np.ndarray[DTYPE_INT, ndim=2] amatrix = read_matrix(matrix)


    if debug: print("initializing memory...", file=sys.stderr)
    # START HERE:
    #cdef int imethod = {"global": 0, "local": 1, "glocal": 2, "global_cfe": 3}[method]
    if imethod == 0:
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
        score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
        score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)
    elif imethod == 3:
        pointer[0, 1:] = LEFT
        pointer[1:, 0] = UP
    elif imethod == 2:
        pointer[0, 1:] = LEFT
        score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)

    if debug: print("filling DP matrices...", file=sys.stderr)

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]
            # agap_i
            agap_i[i,j] = max3(
                         score[i, j - 1] + gap_open,
                         agap_i[i, j - 1] + gap_extend,
                         agap_j[i, j - 1] + gap_double)
            # agap_j
            agap_j[i,j] = max3(
                         score[i - 1, j] + gap_open,
                         agap_j[i - 1, j] + gap_extend,
                         agap_i[i - 1, j] + gap_double)
            # score
            diag_score = score[i - 1, j - 1] + amatrix[ci, cj]
            left_score = agap_i[i, j]
            up_score   = agap_j[i, j]
            max_score = max3(diag_score, up_score, left_score)

            score[i, j] = max2(0, max_score) if method == 'local' else max_score

            #cdef int imethod = {"global": 0, "local": 1, "glocal": 2, "global_cfe": 3}[method]
            if imethod == 1:
                if score[i,j] == 0:
                    pass # point[i,j] = NONE
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT
            elif imethod == 2:
                # In a semi-global alignment we want to consume as much as
                # possible of the longer sequence.
                if max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == diag_score:
                    pointer[i,j] = DIAG
                elif max_score == left_score:
                    pointer[i,j] = LEFT
            else:
                # global
                if max_score == up_score:
                    pointer[i,j] = UP
                elif max_score == left_score:
                    pointer[i,j] = LEFT
                else:
                    pointer[i,j] = DIAG

    if debug: print("finding optimal alignment...", file=sys.stderr)
    if method == 'local':
        # max anywhere
        i, j = max_index(score)
    elif method == 'glocal':
        # max in last col
        i, j = (score[:,-1].argmax(), max_j)
    elif method == 'global_cfe':
        # from i,j to max(max(last row), max(last col)) for free
        row_max, col_idx = score[-1].max(), score[-1].argmax()
        col_max, row_idx = score[:, -1].max(), score[:, -1].argmax()
        if row_max > col_max:
            pointer[-1,col_idx+1:] = LEFT
        else:
            pointer[row_idx+1:,-1] = UP

    if fmt == 'code':
        end_i = i
        end_j = j
        code = []
        last_p = NONE
        counter = 0

        p = pointer[i,j]
        while p != NONE:
            if p == DIAG:
                i -= 1
                j -= 1
            elif p == LEFT:
                j -= 1
            elif p == UP:
                i -= 1
            else:
                raise Exception('wtf!:pointer: %i', p)

            if last_p == NONE:
                last_p = p
                counter += 1
            elif last_p == p:
                counter += 1
            else:
                code.append( (last_p, counter) )
                last_p = p
                counter = 1
            p = pointer[i, j]
        code.append( (last_p, counter) )

        return (i, j, end_i, end_j, max_score, reversed(code))

    if debug: print("preparing alignment...", file=sys.stderr)
    seqlen = max_i + max_j

    ai = PyByteArray_FromStringAndSize(NULL, seqlen)
    aj = PyByteArray_FromStringAndSize(NULL, seqlen)

    # use this and PyObject instead of assigning directly...
    align_j = PyByteArray_AS_STRING(aj)
    align_i = PyByteArray_AS_STRING(ai)

    p = pointer[i, j]
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = seqi[i]
        elif p == LEFT:
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = c"-"
        elif p == UP:
            i -= 1
            align_j[align_counter] = c"-"
            align_i[align_counter] = seqi[i]
        else:
            raise Exception('wtf!:pointer: %i', p)
        #print(align_j[align_counter], align_i[align_counter])
        align_counter += 1
        p = pointer[i, j]

    PyByteArray_Resize(aj, align_counter)
    PyByteArray_Resize(ai, align_counter)

    if flip:
        return (<object>ai)[::-1], (<object>aj)[::-1], max_score
    else:
        return (<object>aj)[::-1], (<object>ai)[::-1], max_score

def score_alignment(object a, object b, int gap_open=-7, int gap_extend=-7, matrix='BLOSUM62'):
    cdef char *al = as_string(a, "expecting bytes or bytearray for first sequence" )
    cdef char *bl = as_string(b, "expecting bytes of bytearray for second sequence" )

    cdef size_t l = strlen(al), i
    cdef int score = 0, this_score
    assert strlen(bl) == l, "alignment lengths must be the same"
    cdef np.ndarray[DTYPE_INT, ndim=2] mat
    mat = read_matrix(matrix)

    cdef bint gap_started = 0

    for i in range(l):
        if al[i] == c"-" or bl[i] == c"-":
            score += gap_extend if gap_started else gap_open
            gap_started = 1
        else:
            this_score = mat[al[i], bl[i]]
            score += this_score
            gap_started = 0
    return score



if __name__ == '__main__':
    # global
    a, b = aligner('WW','WEW', method= 'global')
    assert a == 'W-W'
    assert b == 'WEW'
