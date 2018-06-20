##
## Copyright (c) 2014 Hidayat Trimarsanto <anto@eijkman.go.id> <trimarsanto@gmail.com>
##
## This script runs under seqpy.
## This license of this script is similar to seqpy
##

from seqpy import cout, cerr
from seqpy.core import bioio, funcs
from seqpy.cmds import arg_parser

import sys

def init_argparser():

    p = arg_parser("recircularized dna sequence based on a reference sequence")
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('-r', '--reffile', required=True)
    p.add_argument('--max_mismatch', default=0.05, type=float)
    p.add_argument('--minlen', default=-1, type=int)
    p.add_argument('infile')
    return p



def main( args ):

    circseqs = bioio.multisequence()

    mseq = bioio.load( args.infile, options = args.io_opts )
    rseq = bioio.load( args.reffile )

    for seq in mseq:
        if seq.label != 'NODE_2_length_4501_cov_41.785': continue
        if len(seq) < len(rseq[0]):
            cerr('WARNING: %s is shorter than reference' % seq.label)
        circseq = seq.clone()
        if args.minlen > 0 and len(seq) > args.minlen:
            print('seq:', circseq.label)
            circseq.set_sequence( recircularize_sequence( seq.seq, rseq[0].seq,
                                    max_mismatch = args.max_mismatch ) )
        else:
            circseq.set_sequence( seq.seq )
        circseqs.append( circseq )

    bioio.save( circseqs, args.outfile )


def start_end_pos(seq):

    gap = ord('-')

    start_pos = -1
    for i in range(len(seq)):
        if seq[i] != gap:
            start_pos = i
            break

    end_pos = -1
    for i in range(len(seq)-1, 0, -1):
        if seq[i] != gap:
            end_pos = i
            break

    return (start_pos, end_pos)


def align_ref( ref, query, max_mismatch = 0.5 ):
    """ return:
            ref_start -> position where the reference start aligned
            ref_end -> position where the reference stop aligned

    """

    (rseq, qseq, score) = funcs.align( [ref, query] )

    ref_start, ref_end = start_end_pos( rseq )
    query_start, query_end = start_end_pos( qseq )

    return ref_start, ref_end, query_start, query_end, rseq, qseq, score/min(ref_end-query_start, query_end-ref_start)


def recircularize_sequence( seq, ref, max_mismatch ):

    revcomp = funcs.reverse_complemented(seq)
    ref_start, ref_end, query_start, query_end, arseq, aqseq, score = align_ref( ref, seq, max_mismatch )
    #print(arseq)
    #print(aqseq)
    ref_start2, ref_end2, query_start2, query_end2, arseq2, aqseq2, score2 = align_ref( ref, revcomp, max_mismatch )
    #print(arseq2)
    #print(aqseq2)
    if score2 > score:
        cerr('-> use reverse complement seq')
        ref_start, ref_end = ref_start2, ref_end2
        query_start, query_end = query_start2, query_end2
        seq = revcomp
        arseq = arseq2
        aqseq = aqseq2

    # starting from here, seq = original query seq, qseq, aqseq = aligned query seq

    circularized_seq = aqseq[ref_start:ref_end+1]

    upstream = downstream = ''

    print(arseq)
    print(aqseq)

    if ref_start > query_start:
        # sequence has headings
        upstream = aqseq[0:ref_start]

    if ref_end  < query_end:
        # sequence has downstream / tail
        downstream = aqseq[ref_end+1:]
        #print('boundary ->', aqseq[ref_end-5:ref_end+5])

    print('ustream ->', upstream)
    print('dstream ->', downstream)

    arseq = arseq[ref_start:ref_end+1]

    if len(upstream) > 15:
        merged_1, merged_2, _ = funcs.align( [ arseq, upstream ], degap=False )
        print('merged_1 >< merged_2 >< circularized_seq')
        print(merged_1)
        print(merged_2)
        print(circularized_seq)
        circularized_seq = funcs.merged( [circularized_seq, merged_2] )

    if len(downstream) > 15:
        #print('ref >< circularized_seq')
        #print(arseq)
        #print(circularized_seq)
        merged_3, merged_4, _ = funcs.align( [ funcs.degapped(arseq), downstream ], degap=False )
        print('merged_3 >< merged_4')
        print(merged_3)
        print(merged_4)
        merged_5, merged_6, _ = funcs.align( [merged_3, funcs.degapped(circularized_seq)], degap=False)
        print('merged_5 >< merged_6')
        print(merged_5)
        print(merged_6)
        print('circ >< merged_4 >< merged_6')
        print(circularized_seq)
        print(merged_4)
        print(merged_6)
        circularized_seq = funcs.merged( [merged_4, merged_6] )
        #print('ref >< circularized_seq')
        #print(ref)
        #print(circularized_seq)

    return circularized_seq

