#!/usr/bin/env spcli

from seqpy import cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import posutils


def init_argparser(p=None):
    p = p if p else arg_parser("Filtering zarr/VCF to another zarr file")
    p = posutils.init_argparser(p)

    p.add_argument('-s', '--samplefile', default=None)
    p.add_argument('--keep-only-complete-samples', action='store_true', default=False)
    p.add_argument('--keep-only-complete-variants', action='store_true', default=False)

    p.add_argument('--mindepth', type=int, default=5)
    p.add_argument('--ploidy', type=int, default=2,
                   help='Ploidy of samples (in VCF file), default = 2')
    p.add_argument('--max_alt_alleles', type=int, default=8,
                   help='Maximum number of alternate alleles (in VCF file), default = 8')

    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('infile')

    return p


def zarr2zarr(args):

    from seqpy.core.sgk import sgio, sgutils2 as sgutils
    import pandas as pd

    # check sanity of arguments
    if args.keep_only_complete_samples and args.keep_only_complete_variants:
        cexit("ERR: --keep-only-complete-samples and --keep-only-complete-variants can't be used together")

    # if provided, read posfile
    var_df = posutils.read_posfile(args=args) if args.posfile else None

    # load dataset
    ds = sgio.load_dataset(args.infile, fields=['INFO/*', 'FORMAT/*'],
                           max_alt_alleles=args.max_alt_alleles,
                           ploidy=args.ploidy)

    if var_df:
        cerr('[Checking for duplicate positions]')
        if sgutils.has_duplicate_positions(ds):
            cexit('ERROR: duplicate positions found! '
                  'Please normalize or filter the VCF/Zarr data first.')

    # select variants
    ds = var_df.pos.sel_dataset(ds) if var_df else ds

    # if need to select samples, performed here
    if args.samplefile:
        ds = sgutils.select_samples(ds, samplefile=args.samplefile)

    # further filtering down
    if args.keep_only_complete_samples or args.keep_only_complete_variants:

        depths = ds.call_AD[:, :, :2].sum(axis=2).compute()
        missingness = depths < args.mindepth

        if args.keep_only_complete_samples:
            # sample missingness
            sample_miss = missingness.sum(axis=0)

            sample_df = pd.DataFrame({
                'SAMPLE': ds.sample_id.values,
                'MISS': sample_miss,
            })

            # get only sample with MISS == 0
            complete_sample_df = sample_df.loc[sample_df.MISS == 0]
            complete_samples = sorted(list(complete_sample_df.SAMPLE))

            ds = sgutils.select_samples(ds, samples=complete_samples)

        elif args.keep_only_complete_variants:

            # variant missingness
            variant_miss = missingness.sum(axis=1)

            variant_df = pd.DataFrame({
                'CHROM': ds.contig_id[ds.variant_contig].values,
                'POS': ds.variant_position.values,
                'MISS': variant_miss,
            })

            complete_variant_df = variant_df.loc[variant_df.MISS == 0].reset_index(drop=True)
            ds = complete_variant_df.pos.sel_dataset(ds)

        else:

            raise NotImplementedError()

    sgio.save_dataset(ds, args.outfile, auto_rechunk=True)


def main(args):
    zarr2zarr(args)


# EOF
