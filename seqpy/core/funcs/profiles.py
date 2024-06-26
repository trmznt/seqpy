
from seqpy.core.funcs.funcs import NA_SET, AA_SET, IGN_SET, NA_IUPAC
import numpy

class SeqProfile(object):

    def __init__(self, char_set, ign_set, iupac_set=None, mat=None, apriori=1e-10):
        self.char_set = char_set
        self.ign_set = ign_set
        self.iupac_set = iupac_set
        self.mat = mat
        self.apriori = apriori

    def update(self, mseq):

        char_set = { self.char_set[i]: i for i in range(len(self.char_set)) }

        if self.mat is None:
            import numpy
            charset_len = len(self.char_set)
            max_seqlen = mseq.max_seqlen()
            self.mat = numpy.empty( (max_seqlen, charset_len), numpy.float )
            self.mat.fill(self.apriori)

        cons = self.mat
        for s in mseq:
            seq = s.seq.upper()
            for i in range(0, len(seq)):
                idx = char_set.get(seq[i], -1)
                if idx >= 0: cons[i, idx] += 1
                elif seq[i] in self.ign_set: pass
                elif self.iupac_set and seq[i] in self.iupac_set:
                    chars = self.iupac_set[seq[i]]
                    ratio = 1.0 / len(chars)
                    for c in chars:
                        idx = char_set[c]
                        cons[i, idx] += ratio
                else: raise RuntimeError("Unknown char: %s in sample %s" % (seq[i], s.label))

    def update(self, mseq, weights = None):

        char_set = { self.char_set[i]: i for i in range(len(self.char_set)) }
        if weights == None:
            weights = [1.0] * len(mseq)

        if self.mat is None:
            import numpy
            charset_len = len(self.char_set)
            max_seqlen = mseq.max_seqlen()
            self.mat = numpy.empty( (max_seqlen, charset_len), float )
            self.mat.fill(self.apriori)

        cons = self.mat
        for s, w in zip(mseq, weights):
            seq = s.seq.upper()
            for i in range(0, len(seq)):
                idx = char_set.get(seq[i], -1)
                if idx >= 0: cons[i, idx] += 1*w
                elif seq[i] in self.ign_set: pass
                elif self.iupac_set and seq[i] in self.iupac_set:
                    chars = self.iupac_set[seq[i]]
                    ratio = 1.0 / len(chars)
                    for c in chars:
                        idx = char_set[c]
                        cons[i, idx] += ratio*w
                else: raise RuntimeError("Unknown char: %s in sample %s" % (seq[i], s.label))

    def normalize(self):
        cons = self.mat
        for i in range(0, len(cons)):
            sum_pos = sum( cons[i] )
            if sum_pos < 1.0: continue    # no characters here
            cons[i] = cons[i]/sum_pos

    def consensus(self, threshold = 0.5, ref=None, non_consensus = None,
                  synthetic=False):
        """ return a consensus sequence, set synthetic=True for generating
            synthetic reference sequence (ie. replace deletion/dash to the
            highest base)
        """
        ord_dash = ord('-')

        cons_seq = bytearray()
        cons = self.mat
        char_set = self.char_set
        for i in range(0, len(cons)):
            max_pos = cons[i].argmax()
            if cons[i, max_pos] < threshold:
                if non_consensus:
                    c = non_consensus
                else:
                    c = char_set[max_pos] + 32 if 64 < char_set[max_pos] < 91 else char_set[max_pos]
            else:
                c = char_set[max_pos]
            if synthetic and c == ord_dash:
                # in synthetic mode, use 2nd highest in case this is dash
                c = char_set[cons[i].argsort()[-2]]
            if ref and ref[i] == c:
                c = ord('.')
            cons_seq.append( c )

        return cons_seq

    def log_lk(self, seq):
        """ return log L(this profile | seq) or log P(seq | this profile) """
        # assume seq is aligned with profile !
        log_lk = 0.0
        seq = seq.upper()
        mat = self.mat
        for i in range(0, len(mat)):
            idx = self.char_set.find(seq[i])
            if idx >= 0:
                log_lk += math.log( mat[i, idx] )
        return log_lk


def seq_profile(char_set, ign_set, iupac_set, mseq, apriori=1e-10, normalize=True):
    prof = SeqProfile(char_set, ign_set, iupac_set, None, apriori)
    if mseq:
        prof.update(mseq)
        if normalize:
            prof.normalize()
        assert prof.mat is not None
    return prof

def na_profile(mseq, apriori=1e-10, normalize=True, additional_set = b'', additional_ignore = b''):
    return seq_profile(NA_SET + additional_set, IGN_SET + additional_ignore, NA_IUPAC, mseq, apriori, normalize)

def aa_profile(mseq, apriori=1e-10, normalize=True, additional_set = b'', additional_ignore = b''):
    return seq_profile(AA_SET + additional_set, IGN_SET + additional_ignore, None, mseq, apriori, normalize)


# EOF
