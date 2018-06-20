
from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

from hashlib import sha256

def init_argparser():

    p = arg_parser('perform some simple statistics on a multiple sequence alignment')
    p.add_argument('infile')

    return p


def main( args ):

    mseq = bioio.load( args.infile, options = args.io_opts or [] )

    print('Number of seqs: %d' % len(mseq))

    # get unique haplotype and sample cluster

    haplotypes = {}
    for seq in mseq:
        seq_hash = sha256(seq.seq)
        try:
            haplotypes[seq_hash].append( seq.label )
        except KeyError:
            haplotypes[seq_hash] = [ seq.label ]

    print('Number of unique haplotypes: %d' % len(haplotypes))

    for (idx, item) in enumerate( haplotypes.items() ):
        k, v = item
        print('Haplo %d =>' % idx)
        for label in v:
            print('  %s' % label)


