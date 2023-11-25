
import numpy as np
import xarray as xr
import pandas as pd

from numba import jit, set_num_threads, config, guvectorize
from dask.diagnostics import ProgressBar
from seqpy import cerr


def _allele_for_barcode(dataset, allele_idxes, nonhet_masks=None, nonfailed_masks=None):

    cerr('[Converting allele index to allele genotype...]')
    alleles = xr.apply_ufunc(
        set_alleles,
        allele_idxes,
        dataset.variant_allele,
        input_core_dims=[['samples'], ['alleles']],
        output_core_dims=[['samples']],
        output_dtypes=np.dtype,
        vectorize=True,
        dask='parallelized').compute()

    if nonhet_masks is not None:
        cerr('[Masking het alleles with N...]')
        alleles = alleles.where(nonhet_masks, 'N')

    if nonfailed_masks is not None:
        cerr('[Masking missing alleles with X...]')
        alleles = alleles.where(nonfailed_masks, 'X')

    cerr('[Computing alleles...]')
    return alleles.compute()


def _allele_for_realmccoil(datastore, allele_idxes, nonhet_masks=None, nonfailed_masks=None):

    # for realmccoil, allele_idxes is the alleles
    # make sure we don't have idx > 1, otherwise convert to 1
    cerr('[Converting allele index to TheREALMCCOIL notation...]')
    alleles = allele_idxes.astype(float)
    alleles = alleles.where(alleles <= 1.0, 1.0)

    if nonhet_masks is not None:
        cerr('[Masking het alleles with 0.5...]')
        alleles = alleles.where(nonhet_masks, 0.5)

    if nonfailed_masks is not None:
        cerr('[Masking missing alleles with -1...]')
        alleles = alleles.where(nonfailed_masks, -1)

    return alleles


def _allele_for_hmmibd(datastore, allele_idxes, nonhet_masks=None, nonfailed_masks=None):

    # for hmmibd, allele_idxes is the alleles
    alleles = allele_idxes

    if nonhet_masks is not None:
        # hmmibd does not deal with het_allele, so raise with information
        raise ValueError('hmmIBD can not process het alleles, please use majority option')

    if nonfailed_masks is not None:
        cerr('[Masking missing alleles with -1...]')
        alleles = alleles.where(nonfailed_masks, -1)

    return alleles


def set_alleles(idx, allele):
    return allele[idx]


def calculate_highest_depth(allele_depth, allele_idx):
    # set heterozygosity based on depth ratio
    return np.take_along_axis(
        allele_depth, np.expand_dims(allele_idx, axis=1), axis=1
    ).squeeze(axis=1)


def get_alleles(func, dataset, hetratio=0.67, mindepth=5, minaltdepth=2, useGT=False, threads=-1):
    """ return a list of alleles per variant

        hetratio is the cut-off ratio to call heterozygote, valid values are [0.5, 1.0],
        (major_depths/total_depths) < hetratio will be called hets, values < 0.5 would not
        make sense

        mindepth is the cut-off value for call missing data, depth < mindepth will be mark
        as missing

        minaltdepth is the threshold to call hets, minor_depths >= minaltdepth will be called hets

        hetratio = 0.999 and minaltdepth=2 will indicate that any variants with minor_depth >= 2
        will be called hets, regardless of the hetratio

        to get major allele, use hetratio = -1 and minaltdepth = -1

    """

    if 0 < threads < config.NUMBA_NUM_THREADS:
        set_num_threads(threads)
    else:
        set_num_threads(config.NUMBA_NUM_THREADS)

    ProgressBar().register()

    cerr('[Calculating allele index per variant...]')
    allele_idxes = xr.apply_ufunc(
        np.argmax,
        dataset.call_AD,
        input_core_dims=[['alleles']],
        kwargs={'axis': -1},
        dask='parallelized').compute()

    cerr('[Preparing masks for missing alleles]')
    allele_depths = dataset.call_AD.values
    total_depths = allele_depths.sum(axis=2, where=(allele_depths > 0)) + 1e-3
    nonfailed_masks = (total_depths >= mindepth)

    #import IPython; IPython.embed()
    #total_depths = dataset.call_DP + 1e-3    # this is to prevent "dividing by zero" error
    #nonfailed_masks = (total_depths >= mindepth).compute()

    # import IPython; IPython.embed()

    nonhet_masks = None
    if hetratio > 0 or minaltdepth > 0:
        cerr('[Calculating highest depth per variant...]')
        highest_depths = xr.apply_ufunc(
            calculate_highest_depth,
            dataset.call_AD,
            allele_idxes,
            input_core_dims=[['samples', 'alleles'], ['samples']],
            output_core_dims=[['samples']],
            output_dtypes=np.dtype,
            vectorize=True,
            dask='parallelized').compute()

        if hetratio > 0:
            cerr('[Calculating het ratio...]')
            het_ratios = highest_depths / total_depths
            # if het_ratios >= hetratio cutoff, then we assign as non-het
            nonhet_masks = het_ratios >= hetratio

        if minaltdepth > 0:
            cerr('[Calculating min alt depth...]')
            nonhetminalt_masks = ((total_depths - highest_depths) < minaltdepth).compute()
            if nonhet_masks is None:
                nonhet_masks = nonhetminalt_masks
            else:
                nonhet_masks = (nonhet_masks & nonhetminalt_masks).compute()

    alleles = func(dataset, allele_idxes, nonhet_masks, nonfailed_masks)
    return alleles.compute()


def get_position_tuples(ds):
    cerr('[Generating position labels...]')
    return zip(np.array(ds.contig_id)[ds.variant_contig.values], ds.variant_position.values)


def get_position_ids(ds):
    return [f'{c}:{p}' for c, p in get_position_tuples(ds)]


def has_duplicate_positions(ds):
    positions = list(zip(ds.variant_contig.values, ds.variant_position.values))
    if len(positions) != len(set(positions)):
        return True
    return False


def select_samples(ds, *, samples=None, samplefile=None):

    if samplefile:
        sample_df = pd.read_table(samplefile, header=None)
        samples = sample_df.iloc[:, 0].to_list()

    orig_N = ds.dims['samples']
    ds = ds.sel(samples=ds.sample_id.isin(samples).values)
    curr_N = ds.dims['samples']
    if curr_N != len(samples):
        curr_samples = set(ds.sample_id.values)
        sel_samples = set(samples)
        diff_samples = sel_samples - curr_samples
        raise ValueError(f'Samples not found: {diff_samples}')
    cerr(f'[Subsetting the samples from {orig_N} to {curr_N}]')

    return ds


# EOF
