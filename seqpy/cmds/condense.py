

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.core.funcs import funcs
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("translate dna sequences")
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('infile')
    return p


def main( args ):

    mseq = bioio.load( args.infile, options = args.io_opts )
    cout('reading %d sequences from %s' % (len(mseq), args.infile))
    bioio.save( funcs.condensed( mseq ), args.outfile )

