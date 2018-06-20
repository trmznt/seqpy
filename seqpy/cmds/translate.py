

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.core.funcs import funcs
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("translate dna sequences")
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('--start_codon', default=1, type=int)
    p.add_argument('files', nargs='+')
    return p


def main( args ):

    aaseqs = bioio.multisequence()

    for infile in args.files:

        mseq = bioio.load( infile, options = args.io_opts )
        cout('reading %d sequences from %s' % (len(mseq), infile))

        for seq in mseq:
            aaseq = seq.clone()
            aaseq.set_sequence( funcs.translated(seq, start_pos = args.start_codon ) )
            aaseqs.append( aaseq )

    bioio.save( aaseqs, args.outfile )

