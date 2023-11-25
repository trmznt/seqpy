#!/usr/bin/env spcli

import argparse
import os

from seqpy import cerr, cexit


def init_argparser(p=None):
    p = p if p else argparse.ArgumentParser()

    p.add_argument('-d', '--mindepth', type=int, default=5,
                   help='minimum depth [5]')
    p.add_argument('--max_alt_alleles', type=int, default=8,
                   help='maximum number of altenate alleles [8]')
    p.add_argument('--fraction', action='store_true', default=False,
                   help="add F_MISS column for the calculated fraction")
    p.add_argument('--ploidy', type=int, default=2,
                   help='sample species ploidy')
    p.add_argument('--outsample', default='outsample.tsv',
                   help='output filename for sample missingness')
    p.add_argument('--outvariant', default='outvariant.tsv',
                   help='output filename for variant missingness')
    p.add_argument('infile',
                   help='zarr or VCF file to be qc-ed')

    return p


def zarr2qc(args):

    import dask
    dask.config.set(scheduler='synchronous')

    from seqpy.core.sgk import sgio
    import pandas as pd

    ds = sgio.load_dataset(args.infile, fields=['INFO/*', 'FORMAT/AD'],
                           max_alt_alleles=args.max_alt_alleles,
                           ploidy=args.ploidy)

    depths = ds.call_AD[:, :, :2].sum(axis=2).compute()
    missingness = depths < args.mindepth

    # sample missingness
    sample_miss = missingness.sum(axis=0)
    variant_miss = missingness.sum(axis=1)

    df_samples = pd.DataFrame({
        'SAMPLE': ds.sample_id.values,
        'MISS': sample_miss,
        'VARIANTS': len(missingness.variants)
    })

    df_variants = pd.DataFrame({
        'CHROM': ds.contig_id[ds.variant_contig].values,
        'POS': ds.variant_position.values,
        'MISS': variant_miss,
        'SAMPLES': len(missingness.samples)
    })

    if args.fraction:
        df_variants['F_MISS'] = df_variants.MISS / df_variants.SAMPLES
        df_samples['F_MISS'] = df_samples.MISS / df_samples.VARIANTS

    df_samples.to_csv(args.outsample, sep='\t', index=False)
    df_variants.to_csv(args.outvariant, sep='\t', index=False)

    cerr(f'[Sample missingness written to {args.outsample}]')
    cerr(f'[Variant missingness written to {args.outvariant}]')


def main(args):
    zarr2qc(args)

# EOF
