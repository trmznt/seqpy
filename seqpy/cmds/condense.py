

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.core.funcs import funcs
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("translate dna sequences")
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('-r', '--report', default=None)
    p.add_argument('infile')
    return p

def write_report(cmseq, outreport):

    ref = cmseq[0]
    pos = cmseq.get_control('position').seq.split(b',')
    print(pos)

    for i in range(1, len(cmseq)):
        seq = cmseq[i]
        snps = []
        for r, b, p in zip(ref.seq, seq.seq, pos):
            if b == 78 or b == 110: continue
            if r != b:
                snps.append( p.decode('ASCII') + chr(b) )

        print(seq.label, len(snps), ' '.join(snps))



def main( args ):

    mseq = bioio.load( args.infile, options = args.io_opts )
    cout('reading %d sequences from %s' % (len(mseq), args.infile))
    c_mseq = funcs.condensed( mseq )
    bioio.save( c_mseq, args.outfile )

    if args.report:
        write_report(c_mseq, args.report)

