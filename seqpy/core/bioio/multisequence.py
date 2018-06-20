
from .seqtype import DNA, RNA, PROTEIN
from .biosequence import biosequence


class multisequence(object):


    def __init__(self, label=None):
        self.label = label
        self.seqs = []
        self._type = None
        self.filename = ''
        self._controls = None


    def __getitem__(self, key):
        return self.seqs[key]


    def clear(self):
        del self.seqs[:]
        self.filename = ''
        self._controls = None
        self._type = None


    def append(self, seq):
        """ append a new sequence """
        self.seqs.append( seq )


    def extend(self, multiseq):
        """ extend with multiseq """
        self.seqs.extend( multiseq )


    def __len__(self):
        return len(self.seqs)


    def __delete__(self, key):
        del self.seqs[key]


    def delete(self, key):
        if type(key) == list:
            key.sort()
            key.reverse()
            for k in key:
                del self.seqs[k]
        else:
            del self.seqs[key]


    def max_seqlen(self):
        if self.seqs:
            return max( [ len(s) for s in self.seqs ] )
        return 0


    def insert(self, key, seq):
        self.seqs.insert(key, seq)


    def type(self, nocache=False):
        if self._type is not None and not nocache:
            return self._type

        if len(self) == 0:
            self._type = DNA        # default type is DNA
            return self._type

        p_dna = p_rna = 0.0
        for s in (self[0], self[-1], self[len(self)//2]):
            seq = s.seq.lower()
            a = seq.count(b'a')
            c = seq.count(b'c')
            g = seq.count(b'g')
            t = seq.count(b't')
            u = seq.count(b'u')
            gaps = seq.count(b'-') + seq.count(b'~')
            p_dna += (a + c + g + t) / (len(seq) - gaps)
            p_rna += (a + c + g + u) / (len(seq) - gaps)
        avg_p_dna = p_dna / 3   # len(self)
        avg_p_rna = p_rna / 3
        if avg_p_dna > 0.5 or avg_p_rna > 0.5:
            if avg_p_dna > avg_p_rna:
                self._type = DNA
            else:
                self._type = RNA
        else:
            self._type = PROTEIN

        return self._type


    def set_type(self, seqtype):
        self._type = seqtype


    def summary(self, key):
        pass


    def set_filename(self, filename):
        self.filename = filename


    def align(self, indexes, method=None, matrix=None):
        if matrix is None:
            if self.type() in [ DNA, RNA ]:
                matrix = 'DNA'
            else:
                matrix = 'BLOSUM62'
        src = [ self[idx] for idx in indexes ]

        from seqpy.core import funcs
        results = funcs.align( src, method, matrix )
        for idx, r in zip(indexes, results):
            self[idx].set_sequence( r )


    def clone(self):
        new_mseq = multisequence()
        for s in self:
            new_mseq.append( s.clone() )
        return new_mseq


    def pop(self, key):
        return self.seqs.pop( key )


    def __add__(self, obj):
        self.extend( obj )


    def add_control(self, label, values):
        if self._controls is None:
            self._controls = { label: values }
        else:
            self._controls[label] = values


    def get_control(self, label):
        return self._controls[label]


    def deflate(self):
        """ remove empty sequences (or sequences containing only gaps '-') from list """
        for idx in reversed(range(len(self))):
            seq = self.seqs[idx].seq
            if seq.count(b'-') == len(seq):
                del self.seqs[idx]


    def get_by_label(self, label):
        for s in self:
            if s.label == label:
                return s
        return None


    def sort(self, key=None, reverse=False):
        self.seqs.sort(key=key, reverse=reverse)


    def is_unilen(self):
        l = len(self[0])
        for s in self[1:]:
            if len(s) != l:
                return False
        return True

