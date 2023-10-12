#!/usr/bin/env spcli

import argparse
import os

import pandas as pd

from seqpy import cerr, cexit
from seqpy.core.bioio import posutils


def init_argparser(p=None):
    p = p if p else argparse.ArgumentParser()
    p = posutils.init_argparser(p)

    p.add_argument('-t', '--threads', type=int, default=-1,
                   help='Number of threads to use')
    p.add_argument('-o', '--outfile')
    p.add_argument('--outtarget', default='')
    p.add_argument('--samplefile', default='',
                   help='A headerless text file containing sample code per single line')
    p.add_argument('-m', '--metafile', default='',
                   help='Metafile containing columns to be added to output file, use eg. '
                   'metafile.tsv:SAMPLE,COUNTRY,REGION to add COUNTRY and REGION to output '
                   'using SAMPLE as index')
    p.add_argument('--useGT', default=False, action='store_true')
    p.add_argument('--mindepth', default=5, type=int,
                   help='Cut-off depth to be called missing variant, eg. mindepth = 5 '
                   'indicates variants with total depth < 5 will be mark as missing. '
                   'Default = 5')
    p.add_argument('--hetratio', default=0.67, type=float,
                   help='The ratio of allele depth over total depth to call hets. '
                   'Value of 0.67 means if depth_of_major_allele/total_depth is < 0.67, '
                   'the genotype will be N. Valid values are between 0.5 to 1.0. '
                   'Set to -1 for obtaining major allele. Default = 0.67')
    p.add_argument('--major', default=False, action='store_true',
                   help='Get major allele only, essentially shortcut for "--hetratio -1 --minaltdepth -1"')
    p.add_argument('--minaltdepth', default=2, type=int,
                   help='Threshold value for minor depth of a variant to be called heterozygote, '
                   'eg. minaltdepth = 2 indicates that variants with alternate reads >= 2 will be '
                   'marked as heterozygote, depending on the hetratio. Use hetratio = 0.9999 if '
                   'hetratio is to be ignored. Default = 2')
    p.add_argument('--ploidy', type=int, default=2,
                   help='Ploidy of samples (in VCF file), default = 2')
    p.add_argument('--max_alt_alleles', type=int, default=8,
                   help='Maximum number of alternate alleles (in VCF file), default = 8')
    p.add_argument('infile')

    return p


def pass_for_duplicate_positions(ds):
    positions = list(zip(ds.variant_contig.values, ds.variant_position.values))
    if len(positions) != len(set(positions)):
        return False
    return True


def prepare_dataset(args):

    from seqpy.core.sgk import sgio
    from seqpy.core.bioio import tabutils

    # if provided, read posfile
    posdf = None
    if args.posfile:
        posdf = posutils.read_posfile(args=args)

    # load dataset
    ds = sgio.load_dataset(args.infile, fields=['INFO/*', 'FORMAT/GT', 'FORMAT/AD', 'FORMAT/DP'],
                           max_alt_alleles=args.max_alt_alleles,
                           ploidy=args.ploidy)

    # XXX: check for duplicate position
    cerr('[Checking for duplicate positions]')
    if not pass_for_duplicate_positions(ds):
        cexit('ERROR: duplicate positions found! '
              'Please normalize or filter the VCF/Zarr data first.')

    # select SNPs
    if posdf is not None:
        ds = posdf.pos.sel_dataset(ds)

    # if need to select samples, performed here
    if args.samplefile:
        orig_N = ds.dims['samples']
        sample_df = pd.read_table(args.samplefile, header=None)
        ds = ds.sel(samples=ds.sample_id.isin(sample_df.iloc[:, 0].to_list()))
        curr_N = ds.dims['samples']
        if curr_N != len(sample_df):
            curr_samples = set(ds.sample_id.values)
            sel_samples = set(sample_df.iloc[:, 0])
            diff_samples = sel_samples - curr_samples
            raise ValueError(f'Samples not found: {diff_samples}')
        cerr(f'[Subsetting the samples from {orig_N} to {curr_N}]')

    # prepare metadata to be added to output file
    if args.metafile:
        samples = ds.sample_id.values
        sample_df, errs = tabutils.join_metafile(samples, args.metafile, percenttag=True)
    else:
        sample_df = None

    # convert using hetratio, adjust params

    if args.major:
        args.hetratio = -1
        args.minaltdepth = -1

    return ds, sample_df


def write_tabular(variants, ds, sample_df, args):

    from seqpy.core.bioio import tabutils

    cerr('[Transposing dataarray containing alleles..]')
    variants = variants.transpose().compute()
    tabular_df = tabutils.dataframe_from_variants(ds, variants)

    N, L = tabular_df.shape

    # add additional metadata, which already sorted based on samples
    if sample_df is not None:

        # assume 1st columns of tabular_df is SAMPLE, wchich can be replaced by sample_df
        cerr('[Concatenating tabular alleles with additional sample information')
        tabular_df = pd.concat([sample_df, tabular_df.iloc[:, 1:]], axis=1)

    match ext := os.path.splitext(args.outfile)[1].lower():

        case '.txt':
            with open(args.outfile, 'w') as fout:
                fout.write('#BAR\tSAMPLE\n')
                for idx, row in tabular_df.iterrows():
                    fout.write(''.join(row[1:]) + '\t' + row[0] + '\n')

        case _:
            tabutils.write_file(args.outfile, tabular_df)

    cerr(f'[Genotype data (L={L-1};N={N}) is written to {args.outfile}]')

    if args.outtarget:
        target_posdf = posutils.posframe_from_dataset(ds)
        target_posdf.to_csv(args.outtarget, sep='\t', index=False)
        cerr(f'[Target position is written to {args.outtarget}]')


def zarr2tabular(args):

    import dask
    dask.config.set(scheduler='synchronous')

    from seqpy.core.sgk import sgutils2 as sgutils

    ds, sample_df = prepare_dataset(args)

    #import IPython; IPython.embed()

    cerr('[Converting alleles to tabular format...]')
    variants = sgutils.get_alleles(
        sgutils._allele_for_barcode,
        ds,
        hetratio=args.hetratio,
        mindepth=args.mindepth,
        minaltdepth=args.minaltdepth,
        useGT=args.useGT,
        threads=args.threads
    )

    # saving to ouftile
    write_tabular(variants, ds, sample_df, args)


def main(args):
    zarr2tabular(args)

# EOF
