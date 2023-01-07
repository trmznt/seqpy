#!/usr/bin/env spcli

import argparse

from seqpy import cerr, cexit


# this little utility helps in manipulate position files (BED, pos file)

def init_argparser():
    p = argparse.ArgumentParser(
        description='Small utility to manipulate position files (BED or POS format)')

    p.add_argument('--union', default=False, action='store_true')
    p.add_argument('--intersection', default=False, action='store_true')
    p.add_argument('--difference', default=False, action='store_true')
    p.add_argument('--symdiff', default=False, action='store_true')
    p.add_argument('-o', '--outfile', default=None,
                   help="output filename, if contains .bed then the format will be BED, "
                   "otherwise the format will be POS file (a tab-delimited file with CHROM and POS header")
    p.add_argument('infiles', nargs='+',
                   help="list of position files in either BED or POS format")

    return p


def postopos(args):

    from seqpy.core.bioio import posutils
    pos_dfs = [posutils.read_posfile(infile=infile) for infile in args.infiles]

    if args.union:

        # for now, use posutils combine method, in future use union method

        pos_df = pos_dfs[0]
        for pos2_df in pos_dfs[1:]:
            pos_df = pos_df.pos.combine(pos2_df)

    elif args.intersection:

        pos_df = pos_dfs[0]
        for pos2_df in pos_dfs[1:]:
            pos_df = pos_df.pos.intersection(pos2_df)

    elif args.difference:

        pos_df = pos_dfs[0]
        for pos2_df in pos_dfs[1:]:
            pos_df = pos_df.pos.difference(pos2_df)

    elif args.symdiff:

        pos_df = pos_dfs[0]
        for pos2_df in pos_dfs[1:]:
            pos_df = pos_df.pos.symdiff(pos2_df)

    else:

        cexit('ERR: need either --union, --intersection, --difference or --symdiff')

    if args.outfile:
        if '.bed' in args.outfile:
            # use BED format
            pos_df.pos.to_bed(args.outfile)
        else:
            pos_df.pos.to_pos(args.outfile)

        cerr(f'[Positioni written to {args.outfile}]')


def main(args):
    postopos(args)

# EOF
