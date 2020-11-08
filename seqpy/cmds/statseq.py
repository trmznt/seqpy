#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

from hashlib import sha256

def init_argparser():

    p = arg_parser('perform some simple statistics on a multiple sequence alignment')
    p.add_argument('infile')

    return p


def main( args ):

    statseq( args )

def statseq( args ):

    mseq = bioio.load( args.infile, options = args.io_opts or [] )

    for s in mseq:
        seq = s.seq.upper()
        A_ = seq.count(b'A')
        C_ = seq.count(b'C')
        G_ = seq.count(b'G')
        T_ = seq.count(b'T')
        N_ = seq.count(b'N')
        d_ = seq.count(b'-')
        L = A_ + C_ + G_ + T_ + N_ + d_

        cout('A: %3d  C: %3d  G: %3d  T: %3d  N: %3d  -: %3d  L: %3d  |  \t%s' % (A_, C_, G_, T_, N_, d_, L, s.label))
