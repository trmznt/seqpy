
from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser
from seqpy.core.traceio import traceutils

def init_argparser():

    p = arg_parser('deduplicate multiple sequences by name or content')
    p.add_argument('-o', '--outfile')
    #p.add_argument('--byname', action='store_const', value=True, default=False)
    p.add_argument('files', nargs='+')

    return p


def main(args):

    mseq = bioio.multisequence()
    labels = {}
    idx = 0

    for infile in args.files:
        mseqs = bioio.load( infile )
        for s in mseqs:
            if s.label in labels:
                continue
            mseq.append( s )
            labels[s.label] = idx
            idx += 1

    bioio.save( mseq, args.outfile )

        

