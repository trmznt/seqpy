
import numpy as np


def list2dictidx(l):
    d = {}
    for idx, i in enumerate(l):
        d[i] = idx
    return d


class SNPProfile(object):
    """ this is a saveable object to keep the SNP profile, which is decoupled from
        processing methods so that the processing methods could be modified without
        affecting the stored data
    """

    def __init__(self, code, positions, groups):
        """ code: code name for this profile
            positions: [ (chr, pos, ref, alt), (chr, pos, ref, alt), (<str>, <int>, <str>, <str>), ...]
            groups: [ grp1, grp2, grp3, <str>, ...]
        """

        self.code = code
        self.positions = positions
        self.groups = groups
        self.group2idx = list2dictidx(groups)
        self.pos2idx = list2dictidx(positions)
        # matrix keeps count of both ref and alt alleles
        self.M = np.full(shape=( len(groups), len(positions), 2 ), fill_value=0.1)
        self._logM = None


    def group_index(self):
        if len(self.group2idx) == 0:
            for idx, group in enumerate(self.groups):
                self.group2idx[group] = idx
        return self.group2idx


    def add_data(self, X, Y):

        # filling-up matrices
        # we are iterating row wise
        for i, y in enumerate(Y):
            y_idx = self.group2idx[y]
            n_alt = X[i]
            M_y = self.M[y_idx]

            M_y[n_alt == 0, 0] += 2
            M_y[n_alt == 1] += 1
            M_y[n_alt == 2, 1] += 2

        return self

        # convert M into fractions
        for y in M:
            M_y = M[y]
            total = np.sum(M_y, axis=1)
            M[y] = np.log(M_y / total[:, None])

        self.matrices = M

        return self


    @property
    def logM(self):
        if self._logM is None:
            total = np.sum(self.M, axis=2)
            self._logM = np.log(self.M / total[:,:,None])
        return self._logM


    def group_size(self):
        return len(self.groups)


    def to_pickle(self, outfile):
        """ export this object to a pickle
        """

        d = self.to_dict()
        # save as a pickle

        pickle.save(d, outfile)


    def to_dict(self):
        """ export this object to a dictionary
        """

        d = { 'code': self.code, 'positions': self.positions, 'groups': self.groups}

        # convert M to text and store as d['M']

        return d


    def to_yaml(self, outfile):
        """ export this object to YAML file """

        d = self.to_dict()
        yaml.dump(d, outfile)


class HaplotypeBaseModel(object):

    def __init__(self, model_id, seed=None):
        self.model_id = model_id
        self.logs = []

    def log(self, text):
        self.logs.append( text )

    def flush_log(self):
        log = self.logs
        self.logs = []
        return log

    def select(self, haplotypes, groups):
        pass

    def fit_and_predict(self, X_train, y_train, X_test, k):
        """ return prediction from likelihood and prediction from original method
            (if available) as (lk_prediction, haplotypeprofile)
        """
        pass


    def score(self, genotype_train, group_train, genotype_test, group_test, simid, k_fold):
        """ return dataframe containing score(s)
        """
        pass


class MultiSegmentHaplotypeModel(object):

    def fit_and_predict(self, X_train, y_train, X_test, k):
        """ assume X_train contains a list of segment
        """

        # prepare profile

        #


class HaplotypeRegionProfile(object):

    def __init__(code, groups, positions=None, haplotypes=None):
        self.code = code
        self.positions = positions
        self.haplotypes2idx = {}
        self.haplotypes = [ {} for g in groups ]
        # matrix only keeps count of each haplotype
        #self.M = np.fill( shape=( len(groups), len(haplotypes) ) )
        self._M = None
        self._logM = None


    def add_haplotype(self, haplotype, grp_idx):
        self._M = None
        # XXX: this code contain errors!
        h = haplotype.tobytes()
        try:
            idx = self.haplotype2idx[h]
        except KeyError:
            idx = self.haplotype2idx[h] = len(self.haplotype2idx)
        try:
            self.haplotypes[grp_idx][idx] += 1
        except KeyError:
            self.haplotypes[grp_idx][idx] = 1


    def to_dict(self):

        d = { 'code': self.code, 'positions': self.positions
            , 'haplotypes': self.haplotypes, 'groups': self.groups
            , 'haplotype2idx': self.haplotype2idx
            }
        return d


    @property
    def M(self):
        if self._M is None:
            self._M = np.zeros( shape=( len(haplotypes), len(haplotype2idx)))
            for i in self.haplotypes:
                for j in self.haplotypes[i]:
                    self._M[i,j] = self.haplotypes[i][j]
        return self._M


    @property
    def logM(self):
        if self._logM is None:
            self._logM = np.log(self.M)
        return self._logM


class HaplotypeProfile(object):

    def __init__(self, code, groups = None, profiles = None):
        self.code = code
        self.profiles = profiles
        self.groups = groups
        self.group2idx = list2dictidx(groups)
        self.M = None


    def prepare_profiles(self, segments):
        self.profiles = []
        for segment in segments:
            self.profiles.add( HaplotypeRegionProfile(segment, self.groups.keys()) )


    def add_data(self, haplotypes, group):
        grp_idx = self.group2idx[group]
        for h, profile in zip(haplotypes, self.profiles):
            profile.add_haplotype(h, grp_idx)


    def to_dict(self):

        return    { 'code': self.code
                , 'profiles': [ p.to_dict() for p in self.profiles]
                }

    @classmethod
    def from_dict(cls, d):
        profiles = [ HaplotypeProfile.from_dict(p) for p in d['profiles'] ]
        return cls(d['code'], profiles)


class LikelihoodEstimator(object):

    def __init__(self, profile):
        self.profile = profile
        # convert to log

    def fit(self, X_train, y_train):
        pass

    def predict(self, X_test):
        pass

    def calculate_loglks(self, X):
        """ given a matrix X, will return a matrix M (len(X), len(groups))
            containing the log likelihood
        """
        loglk = np.zero( ( len(X), self.profile.group_size()) )
        for i in range(len(X)):
            for j in range( self.profile.group_size() ):
                loglk[i,j] = self.loglk(X[i], j)

        return loglk


class SNPLikelihoodEstimator(LikelihoodEstimator):

    def __init__(self, profile=None):
        super().__init__(profile)

    def fit(self, X, Y):
        """
        """

        if self.profile != None:
            raise RuntimeError(
                'WARN: cannot execute fit() since profile is already exists, use refit() instead!')

        # check how many groups
        uniq_Y = sorted(set(Y))
        L = len(X[0])

        self.profile = SNPProfile(code='#', positions=list(range(L)), groups=uniq_Y)
        self.profile.add_data(X, Y)

        return self


    def predict(self, X):
        loglks = self.calculate_loglks(X)
        grp_indexes = np.argmax(loglks, axis=1)
        return [ self.profile.groups[idx] for idx in grp_indexes ]


    def calculate_loglks(self, X):
        """ given a matrix X, will return a matrix M (len(X), len(groups))
            containing the log likelihood
        """
        log_lk = np.zeros( ( len(X), self.profile.group_size()) )
        for i in range(len(X)):
            n_alt_data = X[i]
            for j in range( self.profile.group_size() ):
                M = self.profile.logM[j]

                n_alt_0 = 2 * np.sum(M[np.where(n_alt_data == 0), 0])

                # optimize for scarce amount of 1 alt data
                index = np.where(n_alt_data == 1)
                if M[index].size:
                    n_alt_1 = np.sum(M[index])
                else:
                    n_alt_1 = 0

                n_alt_2 = 2 * np.sum(M[np.where(n_alt_data == 2), 1])

                log_lk[i, j] = n_alt_0 + n_alt_1 + n_alt_2

        return log_lk


class HaplotypeLikelihoodEstimator(LikelihoodEstimator):

    def __init__(self, profile = None):
        super().__init__()
        self.profile = profile

    def fit(self, X, y):
        """ X is haplotype matrix as:
            [    ['x0-0', 'x1-0', 'x2-0'],
                ['x0-1', 'x1-1', 'x2-1']
            ]
        """

        if self.profile is not None:
            raise RuntimeError('ERR: cannot use fit() if profile already exists, use refit() instead!')

        uniq_Y = sorted( set(Y) )
        L = len(X[0])

        self.profile = HaplotypeProfile(code=l, groups=uniq_Y)
        self.profile.prepare_profiles([ 'X%d' for d in range(L) ])
        self.profile.add_data(X, Y)

        return self


    def predict(self, X):

        log_lk = self.calculate_loglks(X)
        total_lk = log_lk.sum(axis=1)
        preds = []
        import IPython; IPython.embed()
        for i in range(len(X)):
            preds.append

    def calculate_loglks(self, X):
        """ X is [ [x0-0, x0-1, x0-2], ...] where x0-0, x0-1, x0-2 ... are haplotypes for different haplotyped regions
        """

        log_lk = np.zero( ( len(X), self.profile.group_size(), len(self.profile.profiles)) )
        for feat_idx, (Xn, Pn) in enumerate(zip(X, self.profile.profiles)):
            for i in range(N):
                x = Pn.haplotype2idx[Xn[i]]
                for j in range(self.profile.group_size()):
                    log_lk[i, j, feat_idx] = Pn.logM[j, x]

        return log_lk