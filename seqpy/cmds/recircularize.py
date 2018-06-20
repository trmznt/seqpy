##
## Copyright (c) 2014 Hidayat Trimarsanto <anto@eijkman.go.id> <trimarsanto@gmail.com>
##
## This script runs under seqpy.
## This license of this script is similar to seqpy
##

from seqpy import cout, cerr
from seqpy.core import bioio, funcs
from seqpy.cmds import arg_parser



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
        circseq = seq.clone()
        if args.minlen > 0 and len(seq) > args.minlen:
            print('seq:', circseq.label)
            circseq.set_sequence( recircularize_sequence( seq.seq, rseq[0].seq,
                                    max_mismatch = args.max_mismatch ) )
        else:
            circseq.set_sequence( seq.seq )
        circseqs.append( circseq )

    bioio.save( circseqs, args.outfile )



def map_sequence( short_seq, long_seq, max_mismatch = 0.9 ):
    """ map the short_seq to long_seq,
        return (start_pos, end_pos, inserts, deletions) """

    min_score = len(short_seq) * max_mismatch * 0.75
    (sseq, lseq, score) = funcs.align( [short_seq, long_seq] )

    gap = ord('-')

    start_pos = -1
    for i in range(len(sseq)):
        if sseq[i] != gap:
            start_pos = i
            break

    end_pos = -1
    for i in range(len(sseq)-1, 0, -1):
        if sseq[i] != gap:
            end_pos = i
            break

    if start_pos == end_pos:
        raise RuntimeError()

    inserts = lseq[start_pos:end_pos].count( gap )
    dels = sseq[start_pos:end_pos].count( gap )
    mismatch = 0
    for i in range(start_pos, end_pos):
        if lseq[i] != sseq[i]:
            mismatch += 1

    if mismatch/(end_pos - start_pos) > max_mismatch:
        return (-1, -1, mismatch, -1, -1)

    return (start_pos, end_pos, mismatch, inserts, dels)


def recircularize_sequence( seq, ref, match_len=30, max_mismatch = 0.9 ):
    """ recircularize a circular DNA, based on ref
        return the recircularized sequence
    """

    # first, map ref to seq to find head position(s)

    head_start, _, _, _, _ = map_sequence( ref[:match_len], seq, max_mismatch )
    if head_start < 0:

        seq = funcs.reverse_complemented( seq )
        head_start, _, _, _, _ = map_sequence( ref[:match_len], seq, max_mismatch )
        if head_start < 0:
            cerr('>>> head not found!')
            return seq

    head_start2 = -1
    if len(seq) - head_start >= len(ref):
        # remaining seq is longer than ref, possibly 2nd head pos exists
        # just pass half of ref
        offset = head_start + match_len
        head_start2, _, _, _, _ = map_sequence(
                                ref[:match_len], seq[offset:], max_mismatch )
        if head_start2 >= 0:
            head_start2 += offset

    #print(head_start, head_start2)

    if head_start2 > head_start:
        # this part deals in case where the full sequence appear in the seq

        circularized_seq = seq[head_start:head_start2]
        # create merged with seq[:head_pos] and seq[head_pos2:]

        upstream_part = seq[head_start2:]
        downstream_part = seq[:head_start]
        if len(upstream_part) > 0:
            merged_1, merged_2, _ = funcs.align( [ circularized_seq, upstream_part ] )
            circularized_seq = funcs.merged( [merged_1, merged_2] )
        if len(downstream_part) > 0:
            merged_3, merged_4, _ = funcs.align( [ circularized_seq, downstream_part ] )
            circularized_seq = funcs.merged( [merged_3, merged_4] )

    else:

        upstream_part = seq[head_start:]
        downstream_part = seq[:head_start]
        overlap_len = len(seq) - len(ref) + 15

        overlap_len = min( overlap_len, len(upstream_part), len(downstream_part) )

        if overlap_len == 0:
            return upstream_part + downstream_part

        # map the end of upstream_part to ref
        #print( upstream_part[-overlap_len:] )
        #print( downstream_part[:overlap_len] )
        overlap_start1, overlap_end1, mismatch1, ins1, dels1 = map_sequence( upstream_part[-overlap_len:], ref )
        overlap_start2, overlap_end2, mismatch2, ins2, dels2 = map_sequence( downstream_part[:overlap_len], ref )
        #print(overlap_start1, overlap_end1)
        #print(overlap_start2, overlap_end2)
        if overlap_start2 >= overlap_end1:
            cerr('>> algorithm problem for this sample!')
            return seq
        upstream_tobe_merged, downstream_tobe_merged, score = funcs.align(
                    [ upstream_part[overlap_start2: overlap_end1],
                      downstream_part[:overlap_end1 - overlap_start2] ] )
        merged_segment = funcs.merged( [upstream_tobe_merged, downstream_tobe_merged] )
        circularized_seq = upstream_part[:overlap_start2] + merged_segment + downstream_part[overlap_end1 - overlap_start2:]

    return circularized_seq

