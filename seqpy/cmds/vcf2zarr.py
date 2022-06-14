#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser

from seqpy.core.sgk.sgio import read_vcf, save_dataset


def init_argparser():
    p = arg_parser("Convert VCF to ZARR datastore")
    p.add_argument('-o', '--outfile')
    p.add_argument('--max_alt_alleles', type=int, default=3)
    p.add_argument('infile')
    return p


def vcf2zarr(args):

    cerr(f'I: reading VCF {args.infile}')
    ds = read_vcf(args.infile, max_alt_alleles=args.max_alt_alleles)
    save_dataset(ds, args.outfile)
    cerr(f'I:: datastore written to {args.outfile}')


def main(args):
    vcf2zarr(args)

# EOF
