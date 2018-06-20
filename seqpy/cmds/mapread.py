
from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser
from seqpy.core import funcs
from seqpy.core.traceio import traceutils
from seqpy.core.pwaligner import calign

def init_argparser():

    p = arg_parser("mapping traces to a reference sequence")
    p.add_argument('-o', '--outfile')
    p.add_argument('--ref')
    p.add_argument('--msaout', action='store_true')
    p.add_argument('files', nargs='+')
    return p


def main( args ):

    if not args.ref:
        cerr('reference sequence is needed')
    reference = bioio.load( args.ref, options = args.io_opts )

    traces = [ bioio.load( x ) for x in args.files ]

    if args.msaout:
        map_results = map_to_ref( reference[0], traces, msaout=True)
    else:
        map_results = map_to_ref( reference[0], traces )

    if not args.msaout:

        for res in map_results:
            r = res[2]
            cout( '>%s\n%s\t%d\t%d\t%d\t%d\t%d\t%s' % 
                (res[0].filename, res[1], r[0], r[1], r[2], r[3], r[4], alncode(r[5]))
            )

    else:

        for res in map_results:
            r = res[2]
            cout( '>%s\n%s\n' % (res[0].filename, res[1]) )
            cout( alnprint( r[:-1] ).decode('UTF-16') )

def alncode( listmap ):
    
    sym = {1:'I', 2:'D', 3:'M'}

    return ''.join( '%d%s' % (x[1], sym[x[0]]) for x in listmap )

def alnprint( seqs ):

    max_len = max( [ len(x) for x in seqs ] )

    lines = []

    line_len = 75
    x = 0
    while x < max_len:
        for s in seqs:
            lines.append(s[x:x+line_len])
        lines.append(b"\n")
        x += line_len

    return b'\n'.join( lines )
        



def map_to_ref(ref, traces, msaout=False):

    maps = []
    if not msaout:
        fmt = 'code'
    else:
        fmt = 'alignment'

    for trace in traces:
        # trim first
        upstream, downstream = traceutils.trim( trace, fmt='code' )

        # fwd strand
        bases = trace.edit_bases[upstream:downstream]
        print(bases)
        map_fwd = calign.aligner( ref.seq, bases, -14, -14, method='glocal', matrix='DNAREAD', fmt=fmt )
        print(map_fwd)
        # rev strand
        bases = funcs.reverse_complemented( bases )
        print(bases)
        map_rev = calign.aligner( ref.seq, bases, -14, -14, method='glocal', matrix='DNAREAD', fmt=fmt )
        print(map_rev)
        if fmt == 'code':
            cerr('%s\n%f\t%f' % (trace.filename, map_fwd[4], map_rev[4]))
            if map_fwd[4] > map_rev[4]:
                maps.append( (trace, '+', map_fwd) )
            else:
                maps.append( (trace, '-', map_rev) )
        else:
            cerr('%s\n%f\t%f' % (trace.filename, map_fwd[2], map_rev[2]))
            if map_fwd[2] > map_rev[2]:
                maps.append( (trace, '+', map_fwd) )
            else:
                maps.append( (trace, '-', map_rev) )

    return maps



