import seqpy
import numpy
import array
from subprocess import call
from .traceio import *

def ttuner_processing(filename):
    """ return (basecalls, qualities, matrix [numpy]) """

    # run ttuner
    retcode = call([ seqpy.path_binary('ttuner'), '-p', filename])
    if retcode != 0 and not os.path.exists( filename + '.phd.1' ):
        raise RuntimeError('ERROR: tracetuner failed!!')

    basecalls, bases, qualities = read_phd( filename + '.phd.1' )

    t = traceio.read_trace(filename)

    return ( basecalls, bases, qualities,
                numpy.array( [ t.trace_A, t.trace_C, t.trace_G, t.trace_T ], numpy.float ) )


def ttuner_align(filename, ref_file):
    """ return (basecalls, qualities, matrix [numpy], alignment_codes) """
    pass


def cap_align(filename):
    pass




def read_phd(filename):
    phd_file = open(filename)

    basecalls = array.array('i')
    bases = array.array('c')
    qualities = array.array('i')

    for line in phd_file:
        if line.startswith('BEGIN_DNA'):
            break

    for line in phd_file:
        if line.startswith('END_DNA'):
            break
        base, qual, call = line.split()
        bases.append(base)
        qualities.append(int(qual))
        basecalls.append(int(call))

    return ( numpy.array( basecalls, numpy.int ), bases, numpy.array( qualities, numpy.int ) )


def trim( trace, winsize=4, quality_threshold=20, preserve_case=False, fmt='sequence' ):

    bases = trace.edit_bases
    qualities = trace.edit_qualities

    new_qualities = [ 1 if x > quality_threshold else 0 for x in qualities ]

    # trim upstream & downstream

    upstream_trim = -1
    downstream_trim = -1

    for i in range( len(new_qualities) - winsize ):
        s = sum( new_qualities[i:i+winsize] )
        if s >= winsize:
            upstream_trim = i
            break

    for i in range( len(new_qualities), winsize, -1):
        s = sum( new_qualities[i-winsize:i] )
        if s >= winsize:
            downstream_trim = i
            break

    if fmt == 'code':
        return (upstream_trim, downstream_trim)

    if upstream_trim >= 0 and downstream_trim >= 0:
        b = bases[upstream_trim:downstream_trim]
        q = qualities[upstream_trim:downstream_trim]
        if not preserve_case:
            s = bytearray( len(b) )
            for i in range(len(s)):
                s[i] = b[i] if q[i] > 25 else b[i] + 32
        return (s, q, upstream_trim, downstream_trim)

    return None

        

