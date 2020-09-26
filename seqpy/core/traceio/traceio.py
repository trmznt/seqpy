# traceio.py
'''
Copyright (C) 2004, 2007, 2008 Hidayat Trimarsanto <anto@eijkman.go.id>
Eijkman Institute for Molecular Biology

This module is part of seqpy, a set of python scripts for sequence processing
and is released under GNU Lesser General Public License version 3 or later:
http://www.gnu.org/licenses/gpl.html
'''

__version__ = '20081006a'

import struct
import bisect
import sys, os
import numpy as np
from io import BytesIO
import re

try:
    from seqpy import cerr, cout
except ImportError:
    def cerr(msg):
        print(msg)

    def cout(msg):
        print(msg)

DEBUG = False

def D( text ):
    if DEBUG:
        cerr(text)

class SeqTrace(object):

    def __init__(self):
        self.original_trace = None
        self.num = 0
        self.bases = None
        self.trace_A = None
        self.trace_C = None
        self.trace_G = None
        self.trace_T = None
        self.prob_A = None
        self.prob_C = None
        self.prob_G = None
        self.prob_T = None
        self.basecalls = None
        self.quality = None
        self.revcomp = False
        self.metadata = {}

    def to_seqtrace(self):
        return self

    def get_basecall_index(self, tracepos):
        return bisect.bisect(self.basecalls, tracepos)

    def get_bases(self):
        return self.bases

    def __len__(self):
        return len(self.trace_A)

    def max_length(self):
        return max( len(self.trace_A), len(self.trace_C),
                len(self.trace_G), len(self.trace_T) )

    def max_height(self):
        return max( max(self.trace_A), max(self.trace_C),
                max(self.trace_G), max(self.trace_T) )

    def reverse_complement(self):
        base_dict = dict( a='t', g='c', t='a', c='g', A='T', G='C', T='A', C='G' )
        base_dict['-'] = '-'
        base_list = list(self.bases)
        base_list.reverse()
        rc_bases = ''.join( [ base_dict.get(x, 'N') for x in base_list ] )
        self.bases = rc_bases
        self.trace_A, self.trace_T = list(self.trace_T), list(self.trace_A)
        self.trace_C, self.trace_G = list(self.trace_G), list(self.trace_C)
        if self.prob_A:
            self.prob_A, self.prob_T = list(self.prob_T), list(self.prob_A)
            self.prob_C, self.prob_G = list(self.prob_G), list(self.prob_C)
            for t in [ self.prob_A, self.prob_C, self.prob_G, self.prob_T ]:
                t.reverse()

        for t in [ self.trace_A, self.trace_C, self.trace_G, self.trace_T ]:
            t.reverse()

        # re-set the quality
        if self.prob_A:
            d = {}
            d['a'] = d['A'] = self.prob_A
            d['c'] = d['C'] = self.prob_C
            d['g'] = d['G'] = self.prob_G
            d['t'] = d['T'] = self.prob_T
            d['-'] = d['n'] = d['N'] = self.prob_A
            self.quality = []
            for i in range(0, len(self.bases)):
                try:
                    self.quality.append( d[self.bases[i]][i] )
                except KeyError:
                    self.quality.append( d['a'][i] )
        elif self.quality:
            self.quality = list(self.quality)
            self.quality.reverse()

        # reverse basecalls
        max_tp = self.num
        self.basecalls = list(self.basecalls)
        self.basecalls.reverse()
        rc_basecalls = [ (max_tp - x) for x in self.basecalls ]
        self.basecalls = rc_basecalls

        self.revcomp = True

    def stretch(self):
        avg_spacing = 0.0

    def tracedata(self):
        """ return as TraceData
            trace is arranged as A, C, G, T
        """

        td = _SeqTraceCls()
        td.edit_bases = bytearray(self.bases)
        td.edit_qualities = np.array( self.quality )
        td.edit_basecalling = np.array( self.basecalls )
        td.traces = np.array( (self.trace_A, self.trace_C, self.trace_G, self.trace_T) )
        td.metadata = self.metadata
        td._xscale = 0.95
        td._yscale = 0.10
        td.set_filename( self.filename )

        return td


scf_header_fmt = '>4sLLLLLLLL4sLLLL'


def read_scf_stream(istream):

    bdata = istream.read()

    fmt = scf_header_fmt
    block_size = struct.calcsize(fmt)

    (magic, sample_num, sample_off, base_num, base_leftclip, base_rightclip, base_off,
        comment_size, comment_off, version, sample_size, code_set,
        private_size, private_off) = struct.unpack( fmt, bdata[:block_size] )

    # print('Version:', version)
    # print('Sample size:', sample_size)

    if magic != b'.scf':
        cerr( "Warning: Not an SCF file" )
        sys.exit(1)
    
    trace = SeqTrace()
    trace.num = sample_num
    if sample_size == 1:
        fmt = ">%db" % (sample_num*4)
    else:
        fmt = ">%dh" % (sample_num*4)
    trace_samples = list(struct.unpack( fmt,
        bdata[sample_off : sample_off + struct.calcsize(fmt)] ))

    if version[0] == b'3'[0]:

        trace.trace_A, trace.trace_C, trace.trace_G, trace.trace_T = (
            trace_samples[:sample_num], trace_samples[sample_num:2*sample_num],
            trace_samples[2*sample_num:3*sample_num], trace_samples[3*sample_num:4*sample_num] )

        # get the bases

        offset = base_off
        fmt = ">%dL" % base_num
        datasize = struct.calcsize(fmt)
        trace.basecalls = struct.unpack(fmt, bdata[ offset : offset + datasize ])
        offset += datasize

        fmt = ">%db" % base_num
        datasize = struct.calcsize(fmt)
        trace.prob_A = struct.unpack( fmt, bdata[ offset : offset + datasize ] )
        offset += datasize
        trace.prob_C = struct.unpack( fmt, bdata[ offset : offset + datasize ] )
        offset += datasize
        trace.prob_G = struct.unpack( fmt, bdata[ offset : offset + datasize ] )
        offset += datasize
        trace.prob_T = struct.unpack( fmt, bdata[ offset : offset + datasize ] )
        offset += datasize

        trace.bases = struct.unpack( "%ds" % base_num, bdata[ offset : offset + base_num] )[0]

        def normalize_trace(traces, sample_num):
            p_sample = 0
            for i in range(0, sample_num):
                traces[i] = traces[i] + p_sample
                p_sample = traces[i]
            for i in range(0, sample_num):
                traces[i] = traces[i] + p_sample
                p_sample = traces[i]

        for traces in (trace.trace_A, trace.trace_C, trace.trace_G, trace.trace_T):
            normalize_trace(traces, sample_num)
            min_value = min(traces)
            delta = 0 - min_value
            for i in range(0, sample_num):
                traces[i] += delta

    elif version[0] == b'2'[0]:

        trace.trace_A = trace_samples[0::4]
        trace.trace_C = trace_samples[1::4]
        trace.trace_G = trace_samples[2::4]
        trace.trace_T = trace_samples[3::4]

        fmt = '>L4Bs3B'
        datasize = struct.calcsize(fmt)

        trace.basecalls = []
        trace.prob_A = []
        trace.prob_C = []
        trace.prob_G = []
        trace.prob_T = []
        bases = []

        for i in range(base_num):
            values = struct.unpack( fmt, bdata[base_off + datasize*i: base_off + datasize*(i+1)])
            trace.basecalls.append( values[0] )
            trace.prob_A.append( values[1] )
            trace.prob_C.append( values[2] )
            trace.prob_G.append( values[3] )
            trace.prob_T.append( values[4] )
            bases.append( values[5] )

        trace.bases = b''.join(bases)

    # set the quality
    d = {}
    d[ord('a')] = d[ord('A')] = trace.prob_A
    d[ord('c')] = d[ord('C')] = trace.prob_C
    d[ord('g')] = d[ord('G')] = trace.prob_G
    d[ord('t')] = d[ord('T')] = trace.prob_T
    d[ord('-')] = d[ord('n')] = d[ord('N')] = trace.prob_A
    # print trace.bases
    trace.quality = []
    for i in range(0, base_num):
        try:
            trace.quality.append( d[trace.bases[i]][i] )
        except KeyError:
            trace.quality.append( d[ord('a')][i] )

    comments = bdata[comment_off:comment_off + comment_size]
    if comments:
        for commline in comments.split(b'\n'):
            if not commline: continue
            tag, data = commline.split(b'=', 1)
            trace.metadata[tag] = data

    trace.filename = istream.name
    return trace


def write_scf_stream(trace, ostream):

    ostream.seek(0)

    sample_off = 128
    sample_num = trace.num
    base_off = struct.calcsize('>%dh' % (sample_num*4)) + sample_off
    code_set = 0
    sample_size=2
    version = '3.10'

    def denormalize_trace(traces, sample_num):
        p_delta = 0
        for i in range(sample_num):
            p_sample = traces[i]
            traces[i] = traces[i] - p_delta
            p_delta = p_sample
        p_delta = 0
        for i in range(sample_num):
            p_sample = traces[i]
            traces[i] = traces[i] - p_delta
            p_delta = p_sample

    denorm_traces = []
    for traces in (trace.trace_A, trace.trace_C, trace.trace_G, trace.trace_T):
        _traces = list(traces)
        denormalize_trace(_traces, sample_num)
        denorm_traces += _traces

    #traces = trace.trace_A + trace.trace_C + trace.trace_G + trace.trace_T
    #denormalize_trace(traces, sample_num)

    fmt = ">%dh" % (sample_num*4)
    print( 'sample_num:', sample_num)
    print( 'trace_len:', len(trace.trace_A))
    print( 'traces_len:', len(traces))
    buffers = struct.pack(fmt, *denorm_traces)
    base_off = len(buffers) + sample_off
    
    # base_leftclip & base_rightclip
    # is unused
    base_leftclip = base_rightclip = 0

    base_num = len(trace.basecalls)
    fmt = ">%dL" % base_num
    basecalls = struct.pack(fmt, *trace.basecalls)

    fmt = ">%db" % (base_num * 4)
    basecalls += struct.pack(fmt, *(trace.prob_A + trace.prob_C + trace.prob_G + trace.prob_T))

    basecalls += struct.pack('%ds' % base_num, trace.bases)

    basecalls += struct.pack(">%db" % (base_num * 3), *([0] * (base_num * 3)))

    comment_off = base_off + len(basecalls)
    buffers += basecalls

    comments = '\n'.join( [ "%s=%s" % (k,v) for (k,v) in trace.metadata.items() ] )
    comment_size = len(comments)
    buffers += comments

    private_size = 0
    private_off = 0


    # write header
    ostream.write( struct.pack( scf_header_fmt, '.scf', sample_num, sample_off, 
        base_num, base_leftclip, base_rightclip, base_off,
        comment_size, comment_off, version, sample_size, code_set,
        private_size, private_off) )
    ostream.write('\0' * 72)     # unused/reserved

    # write content
    ostream.write( buffers )


#---------------------------- ABIF format -------------------------------#

abitypes = {
    1: '%dB',
    18: 'x%ds',
    19: '%ds',
    2: '%ds',
    201: '%db',
    4: '>%dh',
    5: '>%df',
    7: '>%df',
    10: '4s',
    11: '4s',
    13: '%db',
    1024: '%ds'
    }

abitags = {
    b'GELP1': 2,
    b'PCON1': 201
}


abif_direntry = '>4slhhll4sl'

class ABIF_DirEntry(object):

    def __init__(self, tag, no, etype, esize, num, dsize, drec, dhdl):
        self.tag = tag
        self.no = no
        self.etype = etype
        self.esize = esize
        self.num = num
        self.dsize = dsize
        self.drec = drec
        self.data = None

    def set_data(self, data):
        self.data = data
        self.num = len(data)

    def get_data(self):
        return self.data


class ABIF(object):

    def __init__(self):
        self.dir_entries = {}
        self.version = None

    def get_entry(self, tagno):
        tag, no = tagno[:4], int(tagno[4:])
        return self.dir_entries[tag][no]

    def get_data(self, tagno):
        return self.get_entry(tagno).get_data()

    def to_seqtrace(self):
        t = SeqTrace()
        t.original_trace = self
        t.num = self.get_entry(b'DATA12').num
        t.bases = self.get_data(b'PBAS1')
        t.basecalls = self.get_data(b'PLOC1')
        try:
            t.quality = self.get_data(b'PCON1')
            t.prob_A = t.prob_C = t.prob_G = t.prob_T = t.quality
        except KeyError:
            pass
        order = self.get_data(b'FWO_1')
        order.upper()
        t.trace_A = self.get_data(('DATA%d' % (9 + order.index(b'A'))).encode('ASCII'))
        t.trace_C = self.get_data(('DATA%d' % (9 + order.index(b'C'))).encode('ASCII'))
        t.trace_G = self.get_data(('DATA%d' % (9 + order.index(b'G'))).encode('ASCII'))
        t.trace_T = self.get_data(('DATA%d' % (9 + order.index(b'T'))).encode('ASCII'))

        return t


def read_abif_stream(istream):

    bdata = istream.read()

    if not bdata.startswith(b'ABIF'):
        raise RuntimeError("Warning: not an ABIF file")

    t = ABIF()
    t.version = struct.unpack('>h', bdata[4:6])[0]

    dir_entry_size = struct.calcsize( abif_direntry )
    header = struct.unpack(abif_direntry, bdata[6: 6 + dir_entry_size])
    dir_entry_num = header[4]
    dir_entry_off = struct.unpack('>l', header[6])[0]

    # read dir_entry and its associated data

    for i in range(0, dir_entry_num):
        offset = dir_entry_off + 28 * i
        elems = struct.unpack( abif_direntry, bdata[ offset : offset + 28 ] )
        de = ABIF_DirEntry( *elems )
        if de.tag in t.dir_entries:
            t.dir_entries[de.tag][de.no] = de
        else:
            t.dir_entries[de.tag] = { de.no : de }

        alt_type = abitags.get(de.tag + str(de.no).encode('ASCII'), de.etype)
        if alt_type != de.etype:
            D( "Warning: inconsistent element type for %s" % de.tag )
        if alt_type == 18: de.num -= 1
        etype_fmt = abitypes.get( alt_type )
        if alt_type not in (10, 11):
            etype_fmt = etype_fmt % de.num
        #D( etype_fmt, alt_type )
        if de.dsize <= 4:
            de.data = struct.unpack( etype_fmt, de.drec[:de.dsize] )
        else:
            offset = struct.unpack('>l', de.drec)[0]
            buf = bdata[offset : offset + de.dsize]
            #print de.tag, de.no, de.etype, de.esize, etype_fmt, de.dsize
            de.data = struct.unpack( etype_fmt, buf )
        if de.num == 1 or alt_type in (18, 19, 2):
            de.data = de.data[0]

    _t = t.to_seqtrace()
    _t.filename = istream.name
    return _t



def read_trace(filename, options={}):

    f = open(filename, 'rb')
    magic = f.read(132)
    f.seek(0)
    if magic.startswith(b'\037\213'):
        import gzip
        f.close()
        f = gzip.open(filename)
        magic = f.read(132)
        f.seek(0)
    if magic.startswith( b'.scf' ):
        return read_scf_stream( f ).tracedata()
    elif magic.startswith( b'ABIF' ):
        return read_abif_stream( f ).tracedata()
    elif magic[128:131] == b'ABIF':
        f.seek(128)
        D( "Warning: MacJunk detected!" )
        return read_abif_stream( f ).tracedata()
    else:
        raise RuntimeError('Supported format: SCF, ABI')


class trace(object):

    def __init__(self):
        self.filename = None
        self.edit_bases = None
        self.edit_qualities = []
        self.edit_basecalling = []
        self.traces = None
        self.metadata = {}

    def set_filename(self, filename):
        self.filename = filename

    def __len__(self):
        return len(self.traces[0])

    def get_basecall_index(self, tracepos):
        return bisect.bisect(self.edit_basecalling, tracepos)

    def as_scf(self):
        ostream = BytesIO()

        sample_off = 128
        sample_num = len(self)
        base_off = struct.calcsize('>%dh' % (sample_num*4)) + sample_off
        code_set = 0
        sample_size=2
        version = b'3.10'

        def denormalize_trace(traces, sample_num):
            p_delta = 0
            for i in range(sample_num):
                p_sample = traces[i]
                traces[i] = traces[i] - p_delta
                p_delta = p_sample
            p_delta = 0
            for i in range(sample_num):
                p_sample = traces[i]
                traces[i] = traces[i] - p_delta
                p_delta = p_sample

        denorm_traces = []
        for i in range(4):
            traces = self.traces[i]
            _traces = list(traces)
            denormalize_trace(_traces, sample_num)
            denorm_traces += _traces

        fmt = ">%dh" % (sample_num*4)
        print('sample_num:', sample_num)
        print('trace_len:', len(self.traces[0]))
        print('traces_len:', len(traces))
        buffers = struct.pack(fmt, *denorm_traces)
        base_off = len(buffers) + sample_off
    
        # base_leftclip & base_rightclip
        # is unused
        base_leftclip = base_rightclip = 0

        base_num = len(self.edit_basecalling)
        fmt = ">%dL" % base_num
        basecall_fmt = [ int(x) for x in self.edit_basecalling ]
        print(fmt)
        print(basecall_fmt)
        basecalls = struct.pack(fmt, *basecall_fmt)

        fmt = ">%db" % (base_num * 4)
        qualities = list(self.edit_qualities) * 4
        basecalls += struct.pack(fmt, *qualities)

        basecalls += struct.pack('%ds' % base_num, self.edit_bases)

        basecalls += struct.pack(">%db" % (base_num * 3), *([0] * (base_num * 3)))

        comment_off = base_off + len(basecalls)
        buffers += basecalls

        comments = b'\n'.join( [ (k + b'=' + v) for (k,v) in self.metadata.items() ] )
        comment_size = len(comments)
        buffers += comments

        private_size = 0
        private_off = 0


        # write header
        ostream.write( struct.pack( scf_header_fmt, b'.scf', sample_num, sample_off, 
            base_num, base_leftclip, base_rightclip, base_off,
            comment_size, comment_off, version, sample_size, code_set,
            private_size, private_off) )
        ostream.write(b'\0' * 72)     # unused/reserved

        # write content
        ostream.write( buffers )

        return ostream.getvalue()

    def boundaries(self, pos):
        midpoint = self.edit_basecalling[pos]
        if pos > 0:
            left = (self.edit_basecalling[pos-1] + midpoint) / 2
        else:
            left = 0

        if pos >= len(self.edit_bases) - 1:
            right = ( len(self) + midpoint ) / 2
        else:
            right = (self.edit_basecalling[pos+1] + midpoint) / 2

        return (left, right)


    def search_re(self, regexp):

        res = []
        for m in re.finditer(regexp, self.edit_bases, re.I):
            res.append( (m.start(), m.end()) )
        return res

    def insert_base(self, pos, base=b'n'):
        # add bases, qualities, basecalls
        midpoint = (self.edit_basecalling[pos] + self.edit_basecalling[pos+1])/2
        pos += 1
        self.edit_bases[pos:pos] = base
        self.edit_qualities = np.insert( self.edit_qualities, pos, 0)
        self.edit_basecalling = np.insert( self.edit_basecalling, pos, midpoint)

    def modify_base(self, pos, base):
        self.edit_bases[pos] = ord(base)

    def name(self):
        print(self.filename)
        if self.filename:
        #    return self.filename
            return os.path.splitext( os.path.basename( self.filename ) )[0]
        return 'untitled-trace'

    def get_sequence(self):
        from seqpy.core.bioio import biosequence
        return biosequence( self.name(), self.edit_bases )

_SeqTraceCls = trace

def set_SeqTraceCls(cls):
    global _SeqTraceCls
    _SeqTraceCls = cls

