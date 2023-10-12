
import os
from seqpy import cerr


try:
    import zarr.storage
    import sgkit
    import sgkit.io.vcf
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


def save_dataset(datastore, outpath):

    cerr(f'INFO: saving dataset to {outpath}')
    store = _get_zarr_storage(outpath, 'w')
    sgkit.save_dataset(datastore, store)
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

# EOF
