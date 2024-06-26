
import os
from seqpy import cerr, cexit
import pandas as pd
import numpy as np

try:
    import zarr.storage
    import sgkit
    import sgkit.io.vcf
    import xarray as xr
except ModuleNotFoundError:
    cerr('ERR: require properly installed sgkit with VCF capability')


def _get_zarr_storage(path, mode='r'):
    match os.path.splitext(path)[1].lower():
        case '.zip':
            return zarr.storage.ZipStore(path, mode=mode)
        case '.lmdb':
            return zarr.storage.LMDBStore(path)
        case _:
            return path


def vcf_to_zarr(vcffile,
                outpath,
                fields=['INFO/*', 'FORMAT/*'],
                max_alt_alleles=3,
                ploidy=2):

    cerr(f'INFO: converting from {vcffile} to {outpath}')
    store = _get_zarr_storage(outpath)
    sgkit.io.vcf.vcf_to_zarr(vcffile, store, fields=fields,
                             max_alt_alleles=max_alt_alleles,
                             ploidy=ploidy,
                             regions=None,
                             target_part_size=None,
                             )
    if hasattr(store, 'close'):
        store.close()


def load_dataset(inpath,
                 max_alt_alleles=8,
                 ploidy=2,
                 fields=['INFO/*', 'FORMAT/*']):
    """ load zarr dataset and reset the encoding in all variables
        so that dataset can be saved using save_dataset
    """

    # in case the inpath is still VCF, just use read_vcf
    if inpath.lower().endswith('.vcf') or inpath.lower().endswith('.vcf.gz'):
        return read_vcf(inpath, fields=fields,
                        max_alt_alleles=max_alt_alleles,
                        ploidy=ploidy)

    cerr(f'INFO: loading dataset from {inpath}')
    store = _get_zarr_storage(inpath, 'r')
    ds = sgkit.load_dataset(store)
    for v in ds.variables:
        ds[v].encoding.clear()
    return ds


def save_dataset(datastore, outpath, auto_rechunk=None):

    cerr(f'INFO: saving dataset to {outpath}')
    store = _get_zarr_storage(outpath, 'w')
    sgkit.save_dataset(datastore, store, auto_rechunk=auto_rechunk)
    if hasattr(store, 'close'):
        store.close()


def read_vcf(path,
             fields=['INFO/*', 'FORMAT/*'],
             max_alt_alleles=3,
             ploidy=2,
             outpath=None):
    cerr(f'INFO: reading VCF from {path}')
    if outpath is None:
        memstore = zarr.storage.MemoryStore()
    else:
        memstore = zarr.storage.DirectoryStore(outpath)
    # use sequential readers, because the vcf_to_zarr (the parallel version)
    # has bugs in reading my (Anto's) VCF files (created by bcftools)
    sgkit.io.vcf.vcf_to_zarr(path, memstore, fields=fields,
                             max_alt_alleles=max_alt_alleles,
                             ploidy=ploidy,
                             regions=None,
                             target_part_size=None,
                             )
    ds = sgkit.load_dataset(memstore)
    for v in ds.variables:
        ds[v].encoding.clear()
    return ds


# -- utilities --


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


def prepare_dataset(ds_or_path, *,
                    posfile=None,
                    samplefile=None,
                    metafile=None,
                    max_alt_alleles=8,
                    ploidy=2,
                    duplicate_position_allowed=False):

    from seqpy.core.bioio import posutils, tabutils

    posdf = None
    if posfile:
        posdf = posutils.read_posfile(posfile)

    if not isinstance(ds_or_path, xr.core.dataset.Dataset):
        cerr(f'Loading dataset from {ds_or_path}')
        ds = load_dataset(ds_or_path, fields=['INFO/*', 'FORMAT/GT', 'FORMAT/AD', 'FORMAT/DP'],
                          max_alt_alleles=max_alt_alleles,
                          ploidy=ploidy)

    if not duplicate_position_allowed and has_duplicate_positions(ds):
        cexit('ERROR: duplicate positions found! '
              'Please normalize or filter the VCF/Zarr data first.')

    if posdf is not None:
        ds = posdf.pos.sel_dataset(ds)

    if samplefile:
        ds = select_samples(ds, samplefile=samplefile)

    if metafile:
        samples = ds.sample_id.values
        sample_df, errs = tabutils.join_metafile(samples, metafile, percenttag=True)
    else:
        sample_df = None

    return ds, sample_df


# EOF
