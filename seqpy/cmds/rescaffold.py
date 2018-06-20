#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.core import bioio, funcs
from seqpy.cmds import arg_parser



def init_argparser():

    p = arg_parser("rescaffold dna sequence based on available contigs and a ref sequence")
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('-r', '--reffile', required=True)
    p.add_argument('--max_mismatch', default=0.1, type=float)
    p.add_argument('contigsfile')
    return p


def main( args ):

    pass


def map_sequences():


    contigs = bioio.load( contigsfile )
    rseq = bioio.load( args.reffile )

    for contig in contigs:

        # map contig to ref sequence
        start, end, mismatch, _, _ = map_sequence( contig, ref, max_mismatch )
        if start < 0:
            contig = funcs.reverse_complemented( contig )
            start, end, mismatch, _, _ = map_sequence( contig, ref, max_mismatch )
            if start < 0:
                continue

        
        
