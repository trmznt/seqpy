

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.core.funcs import funcs
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("translate dna sequences")
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('--start_codon', default=1, type=int)
    p.add_argument('--start_sequence', default=None)
    p.add_argument('files', nargs='+')
    return p


def main( args ):

    aaseqs = bioio.multisequence()

    if args.start_sequence:
        args.start_sequence = args.start_sequence.upper().encode('ASCII')

    for infile in args.files:

        mseq = bioio.load( infile, options = args.io_opts )
        cout('reading %d sequences from %s' % (len(mseq), infile))

        for seq in mseq:
            aaseq = seq.clone()
            if args.start_sequence:
                # we use search restriction pattern function to locate
                # the position
                target_seq = funcs.uppercased(funcs.degapped(seq))
                res = funcs.search_restriction_site(target_seq, args.start_sequence)
                if len(res) != 1:
                    continue
                print(target_seq[res[0][0]:res[0][0]+30])
                aaseq.set_sequence( funcs.translated(target_seq, start_pos = res[0][0]+1) )
            else:
                aaseq.set_sequence( funcs.translated(seq, start_pos = args.start_codon ) )
            aaseqs.append( aaseq )

    bioio.save( aaseqs, args.outfile )

