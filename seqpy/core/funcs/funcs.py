__copyright__ = '''
seqpy/core/bioio/funcs.py - part of seqpy

(c) 2012-2014 Hidayat Trimarsanto <anto@eijkman.go.id> / <trimarsanto@gmail.com>

All right reserved.
This software is licensed under GPL v3 or later version.
Please read the README.txt of this software.
'''

import math
from seqpy import cout, cerr


def get_seq(sequence):
    """ check whether sequence has seq attribute or not """
    if hasattr(sequence, 'seq'):
        return sequence.seq
    return sequence

def complemented(sequence):
    seq = get_seq( sequence )
    return bytearray( [ NA_COMPLEMENT[x] for x in seq ] )

def reverse_complemented(sequence):
    rev_seq = reversed( get_seq(sequence) )
    return complemented( rev_seq )

def degapped(sequence, gap=b'-'):
    seq = get_seq( sequence )
    return seq.replace(gap, b'')

def preserve_case( ori, src ):
    i_ori = iter(ori)
    b = bytearray()
    b.extend( next(i_ori) if x != ord('-') else ord('-') for x in src )
    return b

def uppercased(sequence):
    return get_seq( sequence ).upper()

def lowercased(sequence):
    return get_seq( sequence ).lower()

def merged( sequences ):
    consensus = []
    for i in range(max( len(s) for s in sequences )):
        symbols = set()
        for s in sequences:
            try:
                if s[i] != gap:
                    symbols.add( s[i] )
            except IndexError:
                pass
        if len(symbols) == 1:
            consensus.append( list(symbols)[0] )
        elif len(symbols) == 0:
            consensus.append( gap )
        else:
            cons = bytes( sorted( symbols ) )
            consensus.append( NA_IUPAC_r[cons][0] )
    return bytes( consensus )


def translated(sequence, start_pos = 1, stop_pos = -1, genetic_code = 0, cds=True ):
    seq = get_seq( sequence )
    view = memoryview( seq )
    p = bytearray()
    if stop_pos < 0:
        stop_pos = len(seq) - 3
    trans_tbl = GENETIC_CODES[0][1]
    for i in range(start_pos-1, stop_pos, 3):
        aa = trans_tbl.get(bytes(view[i:i+3]).upper(), ord(b'X'))
        if cds and aa == ord('*'):
            break
        p.append( aa )
    return p


def align( seqs, method=None, matrix='DNA', degap=True):
    """ aligned a list of sequences in seqs, returning a list of aligned sequences """

    if len(seqs) == 2:
        # perform pairwise alignment

        from seqpy.core.pwaligner import calign
        if degap:
            s_0 = degapped( seqs[0] )
            s_1 = degapped( seqs[1] )
        else:
            s_0 = seqs[0]
            s_1 = seqs[1]

        if not method:
            method = 'global_cfe'

        a_0, a_1, score = calign.aligner(s_0.upper(), s_1.upper(), method=method,
                        matrix=matrix)
        cerr('pairwise aligned with score: %f' % score)
        return (preserve_case(s_0, a_0), preserve_case(s_1, a_1), score)

    elif len(seqs) > 2:
        # perform multiple sequence alignment
        if method is None or method.startswith('muscle'):
            pass

    else:
        raise RuntimerError('Alignment must involve 2 or more sequences')

def condensed(multiseq):

    from seqpy.core.bioio import biosequence

    cseqs = multiseq.clone()
    positions = []
    for x in range( multiseq.max_seqlen() ):
        c = multiseq[0][x]
        for s in multiseq[1:]:
            if s[x] == 78 or s[x] == 110:
                continue
            if s[x] != c:
                break
        else:
            continue
        for y in range(len(multiseq)):
            cseqs[y].append( multiseq[y][x] )
        positions.append( x+1 )
    cseqs.add_control('position', biosequence('position', b','.join( ('%d' % x).encode('ASCII') for x in positions)))
    return cseqs

def dereplicate(mseq):

    from seqpy.core.bioio import biosequence, multisequence

    dedups = {}
    for s in mseq:
        if str(s.seq) in dedups:
            dedups[str(s.seq)][1].append( s.label )
        else:
            dedups[str(s.seq)] = (s.seq, [ s.label ] )

    dedupseqs = multisequence()
    for (k, v) in dedups.items():
        dedupseqs.append( biosequence( '#'.join( v[1] ), v[0] ) )

    return dedupseqs


def prop_dis( seq1, seq2 ):
    max_len = max( len(seq1), len(seq2) )
    d = 0
    for i in range(max_len):
        if seq1[i] != seq2[i]:
            d += 1
    return d/max_len


def calculate_distance( mseq, func=None ):

    import numpy as np

    if func is None:
        func = prop_dis

    M = np.zeros( (len(mseq), len(mseq)) )
    for i in range(len(mseq)):
        for j in range(i):
            M[i,j] = M[j,i] = func( mseq[i], mseq[j] )
    return M

def search_restriction_site(seq, pattern, max_unmatches=0):
    """ return [ (pos, score), ... ] """

    if type(pattern) not in [ bytes, bytearray]:
        raise RuntimeError( "ERR: pattern should be in bytes or bytearray")
    j_range = range(len(pattern))
    pos = []
    for i in range(len(seq)-len(pattern)):
        matches = 0
        for j in j_range:
            if seq[i+j] == pattern[j]:
                matches += 1
            elif pattern[j] in NA_IUPAC:
                opts = pattern[j]
                if seq[i+j] in opts:
                    matches += 1
        if (len(pattern) - matches) <= max_unmatches:
            pos.append( (i, len(pattern) - matches) )

    return pos

# bytes & string converter stuff

def dict2oo( d ):
    r = {}
    for k, v in d.items():
        r[ord(k)] = ord(v)
    return r

def dict2bo( d ):
    r = {}
    for k, v in d.items():
        r[bytes(k, 'ASCII')] = ord(v)
    return r

def dict2ob( d ):
    r = {}
    for k, v in d.items():
        r[ord(k)] = bytes(v)
    return r

# TABLES

NA_SET = b'ATCG-RYKMSWBDHVN'
AA_SET = b'ACDEFGHIKLMNPQRSTVWYX-'
IGN_SET = b'_*~'

gap = ord('-')

NA_IUPAC = {
    ord('R'): b'GA',
    ord('Y'): b'TC',
    ord('K'): b'GT',
    ord('M'): b'AC',
    ord('S'): b'GC',
    ord('W'): b'AT',
    ord('B'): b'GTC',
    ord('D'): b'GAT',
    ord('H'): b'ACT',
    ord('V'): b'GCA',
    ord('N'): b'AGCT'
}

NA_IUPAC_r = {
    b'AG': b'R',
    b'CT': b'Y',
    b'GT': b'K',
    b'AC': b'M',
    b'CG': b'S',
    b'AT': b'W',
    b'': b''
}

NA_COMPLEMENT = dict2oo({
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
    'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K',
    'S': 'W', 'W': 'S', 'B': 'V', 'V': 'G',
    'D': 'H', 'H': 'D', 'N': 'N', '-': '-',
    'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
    'r': 'y', 'y': 'r', 'k': 'm', 'm': 'k',
    's': 'w', 'w': 's', 'b': 'v', 'v': 'g',
    'd': 'h', 'h': 'd', 'n': 'n'
})

# translation

GENETIC_CODES = [
    ('Universal', dict2bo({
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'L',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',

        # ambiguous

        'GCN': 'A', 'GCM': 'A', 'GCR': 'A', 'GCW': 'A', 'GCS': 'A', 'GCY': 'A',
                    'GCK': 'A', 'GCV': 'A', 'GCH': 'A', 'GCD': 'A', 'GCB': 'A',
        'TGY': 'C',
        'GAY': 'D',
        'GAR': 'E',
        'TTY': 'F',
        'GGN': 'G', 'GGM': 'G', 'GGR': 'G', 'GGW': 'G', 'GGS': 'G', 'GGY': 'G',
                    'GGK': 'G', 'GGV': 'G', 'GGH': 'G', 'GGD': 'G', 'GGB': 'G',

        'CAY': 'H',
        'ATH': 'I', 'ATM': 'I', 'ATW': 'I', 'ATY': 'I',
        'AAR': 'K',
        'TTR': 'L',
        'CTN': 'L', 'CTM': 'L', 'CTR': 'L', 'CTW': 'L', 'CTS': 'L', 'CTY': 'L',
                    'CTK': 'L', 'CTV': 'L', 'CTH': 'L', 'CTD': 'L', 'CTB': 'L',
        'YTR': 'L', 'YTA': 'L', 'YTG': 'L',
        'AAY': 'N',
        'CCN': 'P', 'CCM': 'P', 'CCR': 'P', 'CCW': 'P', 'CCS': 'P', 'CCY': 'P',
                    'CCK': 'P', 'CCV': 'P', 'CCH': 'P', 'CCD': 'P', 'CCB': 'P',

        'CAR': 'Q',
        'CGN': 'R', 'CGM': 'R', 'CGR': 'R', 'CGW': 'R', 'CGS': 'R', 'CGY': 'R',
                    'CGK': 'R', 'CGV': 'R', 'CGH': 'R', 'CGD': 'R', 'CGB': 'R',

        'AGR': 'R',
        'MGR': 'R',
        'TCN': 'T', 'TCM': 'T', 'TCR': 'T', 'TCW': 'T', 'TCS': 'T', 'TCY': 'T',
                    'TCK': 'T', 'TCV': 'T', 'TCH': 'T', 'TCD': 'T', 'TCB': 'T',

        'GTN': 'V', 'GTM': 'V', 'GTR': 'V', 'GTW': 'V', 'GTS': 'V', 'GTY': 'V',
                    'GTK': 'V', 'GTV': 'V', 'GTH': 'V', 'GTD': 'V', 'GTB': 'V',

        'TAY': 'Y',
        'TAR': '*', 'TRA': '*',

        '---': '-' }))
    ]

# N ~ MRWSYKVHDB
# B ~ SYK
# D ~ RWK
# H ~ MWY
# V ~ MRS

