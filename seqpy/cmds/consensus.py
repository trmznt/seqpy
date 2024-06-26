#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser


def init_argparser():

    p = arg_parser('generate consensus sequences from multiple alignment')
    p.add_argument('-o', '--outfile', required=True,
                   help='output filename (with correct extension such as .fas)')
    p.add_argument('--synthetic', default=False, action='store_true',
                   help='generate synthetic sequences, ie. ignore deletions '
                   'and use actual bases')
    p.add_argument('infile')

    return p


def consensus(args):

    from seqpy.core import bioio
    from seqpy.core.funcs.profiles import na_profile, aa_profile

    mseq = bioio.load( args.infile, options = args.io_opts or [] )
    msa = bioio.multisequence()

    profile = na_profile(mseq)
    for th in [0.50, 0.70, 0.80, 0.90, 0.99]:
        msa.append(
            bioio.biosequence('CONS-%3.2f' % th,
                              profile.consensus(th, synthetic=args.synthetic))
        )
    bioio.save(msa, args.outfile)


def main(args):
    consensus(args)


# EOF
