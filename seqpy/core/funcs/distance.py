
from seqpy import cerr

from math import isqrt
import numpy as np
import sklearn.metrics


# this distance module provides functions to handle proportional genetic distance
# alleles can be coded to up 9999 types

MISSING_CODE = 255


def xxx_calculate_diff_N(X, Y):
    # return pair of (diff, N) as Szudzik's pairing function

    L = X.shape[0]
    X_miss = (X == MISSING_CODE)
    Y_miss = (Y == MISSING_CODE)
    any_missing = (X_miss | Y_miss).sum()
    one_missing = (X_miss ^ Y_miss).sum()
    total_diff = (X != Y).sum()

    diff = total_diff - one_missing
    N = L - any_missing
    #NaN = np.count_nonzero(diff_array > MISSING_THRESHOLD)
    #diff = np.count_nonzero(diff_array != 0) - NaN
    #total = N - NaN

    # return as Szudzik's pairing function
    return (N**2 + diff)


def calculate_diff_N(X, Y):
    # return pair of (diff, N) as Szudzik's pairing function

    L = X.shape[0]
    X_miss = (X == MISSING_CODE)
    Y_miss = (Y == MISSING_CODE)
    any_missing = np.count_nonzero(X_miss | Y_miss)
    one_missing = np.count_nonzero(X_miss ^ Y_miss)
    total_diff = np.count_nonzero(X != Y)

    diff = total_diff - one_missing
    N = L - any_missing

    # return as Szudzik's pairing function
    return (N**2 + diff)


def calculate_diff_N_with_0(X, Y):
    # return pair of (diff, N) as Szudzik's pairing function
    # this only works if MISSING_CODE == 0
    # note: this does not work yet!!

    # X = [1, 2, 1, 4, 0, 1, 0]
    # Y = [1, 2, 4, 1, 1, 0, 0]

    # X | Y = [1, 2, 5, 5, 1, 1, 0]
    # X & Y = [1, 2, 0, 0, 0, 0, 0]
    # X ^ Y = [0, 0, 5, 5, 1, 1, 0]

    # np.count_nonzero(X ^ Y) == np.count_nonzero(X != Y)

    L = X.shape[0]
    X_miss = (X == MISSING_CODE)
    Y_miss = (Y == MISSING_CODE)
    any_missing = L - np.count_nonzero(X | Y)
    one_missing = L - np.count_nonzero(X ^ Y)
    total_diff = np.count_nonzero(X != Y)

    diff = total_diff - one_missing
    N = L - any_missing
    #NaN = np.count_nonzero(diff_array > MISSING_THRESHOLD)
    #diff = np.count_nonzero(diff_array != 0) - NaN
    #total = N - NaN

    # return as Szudzik's pairing function
    return (N**2 + diff)


def decode_distmatrix(distm):
    import numpy as np

    isqrt_vec = np.vectorize(isqrt)
    N = isqrt_vec(distm)
    d = distm - N**2
    return d, N


def pairwise_distances(alleles):
    return sklearn.metrics.pairwise_distances(alleles, metric=calculate_diff_N,
                                              n_jobs=-1,
                                              force_all_finite=True).astype(int)


def AD_to_alleles(allele_depths, min_depth=5, allele_number=-1):
    """ return allele array from xarray allele_depths """

    if allele_number > 0:
        allele_depths = allele_depths[:, :, :allele_number]

    allele_depths[allele_depths < 0] = 0
    alleles = allele_depths.argmax(axis=2)

    # calculate total depths and get missing mask (True == missing)
    variant_depths = allele_depths.sum(axis=2)
    missing_mask = variant_depths < min_depth
    alleles[missing_mask] = MISSING_CODE

    return alleles.T


def any_to_alleles(genotypes, missing=np.nan):
    """ return allele array from any object (char/string/numeric) genotypes """

    translation_dict = {}
    cerr('[Finding unique allele values...]')
    items = np.unique(genotypes)
    cerr('[Building allele translation table...]')
    for idx, item in enumerate(items, 1):
        # we pass missing because the new allele matrix will originally be filled with missing codes
        if item is missing:
            continue
        translation_dict[item] = idx
    alleles = np.full(genotypes.shape, MISSING_CODE, np.uint8)
    for k, i in translation_dict.items():
        cerr(f'[Converting allele {k} to {i}...]')
        alleles[genotypes == k] = i

    return alleles


# EOF
