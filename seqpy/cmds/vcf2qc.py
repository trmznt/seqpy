#!/usr/bin/env spcli

import argparse
import os

from seqpy import cerr, cexit

# this utility calculate stats using GT field of VCF files

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
                   help='output filename for sample quality stats')
    p.add_argument('--outvariant', default='outvariant.tsv',
                   help='output filename for variant quality stats')
    p.add_argument('infile',
                   help='zarr or VCF file to be qc-ed')

    return p


def vcf2qc(args):

    import sys
    import numpy as np
    import pandas as pd
    from cyvcf2 import VCF, Writer

    if not args.infile:
        args.infile = '-'

    vcf = VCF(args.infile)

    samples = vcf.samples
    sample_len = len(samples)
    sample_gt_types = np.zeros((sample_len, 4), dtype=int)

    # for variants
    CHROM = []
    POS = []
    MISS = []
    SAMPLES = []
    HETS = []
    HOMREF = []
    HOMALT = []

    for v in vcf:

        GT = v.format('GT')
        genotypes = v.genotypes
        # 0: ref, 1: alt, -1: missing

        gt_types = v.gt_types
        # 0: homo ref, 1: het 2: missing, 3: homo alt

        # for variant quality, store the stats directly
        CHROM.append(v.CHROM)
        POS.append(v.POS)
        #SAMPLES.append(v.num_called)
        HOMREF.append(v.num_hom_ref)
        HOMALT.append(v.num_hom_alt)
        HETS.append(v.num_het)
        MISS.append(v.num_unknown)

        # for sample quality, store stats to sample holders
        # this can rewritten in numpy-style code
        for s_idx in range(len(v.gt_types)):
            sample_gt_types[s_idx, int(v.gt_types[s_idx])] += 1

    # generate data frames

    df_variants = pd.DataFrame({
        'CHROM': CHROM,
        'POS': POS,
        'MISS': MISS,
        'SAMPLES': sample_len,
        'HOMREF': HOMREF,
        'HOMALT': HOMALT,
        'HETS': HETS,
    })

    variant_len = len(df_variants)
    df_samples = pd.DataFrame({
        'SAMPLE': samples,
        'HOMREF': sample_gt_types[:, 0],
        'HETS': sample_gt_types[:, 1],
        'HOMALT': sample_gt_types[:, 3],
        'MISS': sample_gt_types[:, 2],
        'VARIANTS': variant_len,
    })

    if args.fraction:

        df_variants['F_MISS'] = df_variants.MISS / df_variants.SAMPLES
        df_variants['F_HETS'] = df_variants.HETS / df_variants.SAMPLES
        df_variants['F_HOMREF'] = df_variants.HOMREF / df_variants.SAMPLES
        df_variants['F_HOMALT'] = df_variants.HOMALT / df_variants.SAMPLES

        df_samples['F_MISS'] = df_samples.MISS / df_samples.VARIANTS
        df_samples['F_HETS'] = df_samples.HETS / df_samples.VARIANTS
        df_samples['F_HOMREF'] = df_samples.HOMREF / df_samples.VARIANTS
        df_samples['F_HOMALT'] = df_samples.HOMALT / df_samples.VARIANTS

    df_samples.to_csv(args.outsample, sep='\t', index=False)
    df_variants.to_csv(args.outvariant, sep='\t', index=False)

    cerr(f'[Sample quality stats written to {args.outsample}]')
    cerr(f'[Variant quality stats written to {args.outvariant}]')


def main(args):
    vcf2qc(args)


# EOF
