
from typing import Hashable

import numpy as np
import xarray as xr
import sgkit as sg

from numba import guvectorize

from sgkit.utils import (
    conditional_merge_datasets,
    create_dataset,
    define_variable_if_absent,
)

# get smallest floating point constant to prevent divided-by-zero
SMALLEST_INCR = np.nextafter(0, 1)


def variant_allele_lengths(
    ds: xr.Dataset,
    *,
    merge: bool = True,
) -> xr.Dataset:

    nplen = np.frompyfunc(len, 1, 1)

    new_ds = xr.Dataset(
        {"variant_allele_lengths":
            xr.apply_ufunc(nplen,
                           ds['variant_allele'],
                           output_dtypes=[int],
                           dask='parallelized')
         }
    )
    return conditional_merge_datasets(ds, new_ds, merge)


@guvectorize(
    [
        "(uint64[:], uint64[:])"
    ],
    "(n) -> ()"
)
def nptally_gufunc(x, out):
    out[0] = (x > 0).sum()


def variant_allele_tally(
    ds: xr.Dataset,
    *,
    merge: bool = True,
) -> xr.Dataset:

    new_ds = xr.Dataset(
        {"variant_allele_tally":
            xr.apply_ufunc(nptally_gufunc,
                           ds['variant_allele_count'],
                           input_core_dims=[['alleles']],
                           output_dtypes=[int],
                           dask='parallelized')
         }
    )
    return conditional_merge_datasets(ds, new_ds, merge)


def variant_biallelic(
    ds: xr.Dataset,
    *,
    biallelic_ratio: Hashable = "variant_biallelic_ratio",
    biallelic_ratio_threshold: float = 1.0,
    merge: bool = True,
) -> xr.Dataset:

    ds = define_variable_if_absent(
        ds, "variant_biallelic_ratio", biallelic_ratio, variant_biallelic_ratio
    )

    new_ds = xr.Dataset(
        {"variant_biallelic": (ds['variant_biallele_ratio'] >= biallelic_ratio_threshold)}
    )
    return conditional_merge_datasets(ds, new_ds, merge)


@guvectorize(
    [
        "(uint64[:], float64[:])"
    ],
    "(n) -> ()"
)
def biallelic_ratio_gufunc(x, out):
    out[0] = np.sort(x)[-2:].sum() / (x.sum() + SMALLEST_INCR)


def variant_biallelic_ratio(
    ds: xr.Dataset,
    *,
    merge: bool = True,
) -> xr.Dataset:

    """ Compute ratio of two highest allele count over total allele count,
        suitable for inferring biallelic variants.
        Strict biallelic variants is when biallelic ratio is 1
    """

    new_ds = xr.Dataset(
        {"variant_biallelic_ratio":
            xr.apply_ufunc(biallelic_ratio_gufunc,
                           ds['variant_allele_count'],
                           input_core_dims=[['alleles']],
                           output_dtypes=[float],
                           dask='parallelized')
         }
    )
    return conditional_merge_datasets(ds, new_ds, merge)


def variant_indel_ratio(
    ds: xr.Dataset,
    *,
    allele_lengths: Hashable = "variant_allele_lengths",
    merge: bool = True,
) -> xr.Dataset:

    """ Compute ratio of indel allele count over total allele count,
        suitable for inferring between SNV and indels.
        Strict SNV will have ratio of 0.0
    """

    ds = define_variable_if_absent(
        ds, "variant_allele_lengths", allele_lengths, variant_allele_lengths
    )

    var_allele_count = ds['variant_allele_count']
    var_allele_lengths = ds['variant_allele_lengths']
    indel_count = var_allele_count.where(var_allele_lengths > 1).sum(dim='alleles')

    new_ds = xr.Dataset(
        {"variant_indel_ratio": indel_count / (ds['variant_allele_total'] + SMALLEST_INCR)}
    )
    return conditional_merge_datasets(ds, new_ds, merge)


def variant_stats(
    ds: xr.Dataset
) -> xr.Dataset:

    """ Compute all variant-based statistics """

    _ds = sg.variant_stats(ds)
    _ds = variant_allele_lengths(_ds)
    _ds = variant_allele_tally(_ds)
    _ds = variant_biallelic_ratio(_ds)
    _ds = variant_indel_ratio(_ds)

    return _ds
# EOF
