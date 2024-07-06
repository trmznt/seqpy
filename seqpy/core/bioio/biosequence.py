
from seqpy.core import funcs

## Design Consideration:
##
## biosequence label, definition and attr should be Unicode strings (UTF-8 or ASCII encoding)
## biosequence.seq should be bytearray (or bytes)

class biosequence(object):

    __slots__ = [ 'label', 'seq', 'attr', 'definition', '_reverse_complement',
                    '_gap_pos', '_gap_ext', '_extdata', 'rec' ]

    def __init__(self, label, seq=b''):
        self.label = label              # string
        self.seq = bytearray(seq)       # bytearray
        self.attr = None                # string
        self.definition = None          # string
        self._reverse_complement = False
        self._gap_pos = self._gap_ext = None
        self._extdata = None
        self.rec = None

    def set_label(self, label):
        self.label = label
        return self

    def set_sequence(self, seq):
        assert type(seq) in [ bytes, bytearray ], "Sequence must be either bytes or bytearray"
            
        self.seq = seq
        if not self._gap_pos is None:
            self.init_listgap()
        return self

    def append(self, symbol):
        """ appending one symbol """
        self.seq.append( symbol )

    def extend(self, seq):
        """ extending the sequence """
        self.seq.extend( seq )

    def count(self, symbol):
        self.seq.count( symbol )

    def insert(self, idx, seq):
        """ inserting sequence """
        self.seq[idx:idx] = seq

    def add_attr(self, key, value):
        if value.startswith('"') and value.endswith('"'):
            value = value[1:-1]
        if self.attr is None:
            self.attr = { key: value }
        else:
            self.attr[key] = value

    def __delitem__(self, key):
        """ deleting part of sequence """
        del self.seq[key]

    def __repr__(self):
        return '<biosequence: (%d) %s - [%s...]>' % (len(self.seq), self.label, self.seq[:25].decode('ASCII'))

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]

    def __setitem__(self, key, value):
        self.seq[key] = value

    # FUNCTIONS

    def reverse_complement(self):
        self.seq = funcs.reverse_complemented(self.seq)
        self._reverse_complemented = False if self._reverse_complemented else True
        return self

    def reverse_complemented(self):
        return self.clone().set_sequence( funcs.reverse_complemented(self.seq) )

    def degap(self):
        self.seq = self.seq.replace(b'-', b'')
        return self

    def degapped(self):
        return self.clone().set_sequence( self.seq.replace(b'-', b'') )

    def clone(self, items=None):
        s = biosequence( self.label )
        if s.attr:
            s.attr = copy(self.attr)
        return s

    def init_listgap(self):
        self._gap_pos, self._gap_ext = list_gap(self.seq)

    def relative_pos(self, pos):
        if self._gap_pos == None:
            self.init_listgap()

        if not self._gap_pos:
            return pos

        t = 0
        for g, e in zip(self._gap_pos, self._gap_ext):
            if pos < g:
                break
            if pos < g + e:
                if g == 0:
                    return -1
                #if pos < g:
                #    break
                return g - t - 1
            t += e
        return pos - t

    def absolute_pos(self, pos):
        if self._gap_pos == None:
            self.init_listgap()

        if not self._gap_pos:
            return pos

        t = 0
        for g, e in zip( self._gap_pos, self._gap_ext ):
            if pos < g - t:
                break
            t += e
        return pos + t

    def upper(self):
        self.seq = self.seq.upper()
        return self


def list_gap(seqstr):
    """ return (gap_pos, gap_ext) """
    gap_pos = []
    gap_ext = []
    if not seqstr:
        return (gap_pos, gap_ext)
    next_pos = 0
    while True:
        gap = seqstr.find(b'-', next_pos)
        if gap < 0:
            break
        extend = 1
        try:
            while seqstr[gap+extend] == ord('-'):
                extend += 1
        except IndexError:
            pass
        gap_pos.append( gap )
        gap_ext.append( extend )
        next_pos = gap + extend
    return (gap_pos, gap_ext)

