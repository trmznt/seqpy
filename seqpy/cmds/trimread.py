
from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser
from seqpy.core.traceio import traceutils

def init_argparser():

    p = arg_parser('trimming capillary read and save as fasta file (fastq in future)')
    p.add_argument('-o', '--outfile')
    p.add_argument('files', nargs='+')
    return p


def main(args):

    mseq = bioio.multisequence()

    for infile in args.files:
        trace = bioio.load( infile )
        result = traceutils.trim( trace )
        if not result:
            continue

        bases, quals, upstream_trim, downstream_trim = result
        seq = bioio.biosequence( infile, bases )
        seq.add_attr( 'upstream_trim', str(upstream_trim) )
        seq.add_attr( 'downstream_trim', str(downstream_trim) )

        mseq.append( seq )

    bioio.save( mseq, args.outfile )
        
