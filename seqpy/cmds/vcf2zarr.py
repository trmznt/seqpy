#!/usr/bin/env spcli

from seqpy import cerr
from seqpy.cmds import arg_parser


def init_argparser():
    p = arg_parser("Convert VCF to ZARR datastore")
    p.add_argument('-o', '--outfile')
    p.add_argument('--ploidy', type=int, default=2)
    p.add_argument('--max_alt_alleles', type=int, default=3)
    p.add_argument('--fields', default='FORMAT/GT,FORMAT/AD',
                   help='Extra fields to be extracted from VCF file, eg: INFO/*, '
                   'FORMAT/AD, FORMAT/GT, FORMAT/*, etc. [FORMAT/GT,FORMAT/AD]')
    p.add_argument('infiles', nargs='+')
    return p


def vcf2zarr(args):

    from seqpy.core.sgk.sgio import read_vcf, save_dataset
    from sgkit.io.vcf import vcf_to_zarr

    fields = args.fields.split(',')

    if len(args.infiles) > 1:
        from dask.distributed import Client

        client = Client(n_workers=16, threads_per_worker=1)
        cerr('[Working in parallel]')
        vcf_to_zarr(args.infiles, args.outfile, tempdir=args.outfile + '-temp',
                    max_alt_alleles=args.max_alt_alleles,
                    fields=fields,
                    ploidy=args.ploidy)

    #elif not args.outfile.endswith('.zip'):
    #    vcf_to_zarr(args.infiles[0], args.outfile,
    #                fields=fields,
    #                max_alt_alleles=args.max_alt_alleles,
    #                ploidy=args.ploidy)

    else:
        cerr(f'[Reading VCF {args.infiles[0]}]')
        ds = read_vcf(args.infiles[0],
                      max_alt_alleles=args.max_alt_alleles,
                      fields=fields,
                      ploidy=args.ploidy)
        cerr(f'[Writing to {args.outfile}]')
        save_dataset(ds, args.outfile)

    cerr(f'[Datastore written to {args.outfile}]')


def main(args):
    vcf2zarr(args)

# EOF
