# Likelihood evaluator

import numpy as np
import pandas as pd

import itertools, time, math, os, pickle, datetime

import numpy as np, pandas as pd, pyparsing as pp
from sklearn.model_selection import StratifiedKFold
from sklearn.tree import DecisionTreeClassifier
import allel

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser



class BaseSelector(object):

    code = 'base'

    def __init__(self, model_id, k=None, snpindex=None, iteration=1, seed=None):
        self.model_id = model_id
        self.k_list = k or [0]
        self.iteration = iteration
        self.snpindex = snpindex
        self.logs = []
        self.reseed( seed )


    def reseed(self, seed):
        self.randomstate = np.random.RandomState( seed )


    def select(self):
        """ return list of SNP list """
        pass


    def fit_and_predict(self, X_train, y_train, X_test, k):
        """ return prediction from likelihood and prediction from original method (if available) as
            (lk_prediction, snp_list, model_prediction)
        """

        # select SNPs based on model
        L, orig_predictions, params = self.select(X_train, y_train, X_test, k)
        if len(L) <= 0:
            return None, None, None, None
        X_train_0 = X_train[:,L]
        X_test_0 = X_test[:,L]

        # create, train and predict using LikelihoodProfile

        lk_est = SNPLikelihoodEstimator()
        lk_est.fit(X_train_0, y_train)

        return lk_est.predict(X_test_0), L, orig_predictions, params


    def score(self, genotype_train, group_train, genotype_test, group_test, simid, k_fold):
        """ return a dataframe containing scores, dict of snps, logs, and prediction result """

        results = []
        snps = {}
        log = []
        preds = []

        for k in self.k_list:

            # best score containe (F, score_dataframe, orig_score_dataframe, snplist)
            best_score = (-1, None, None, None)

            for i in range(self.iteration):

                lk_pred, snplist, orig_pred, params = self.fit_and_predict( genotype_train
                        , group_train, genotype_test, k)
                if lk_pred is None:
                    continue

                if len(snplist) > 0:
                    begin_snp = min(snplist)
                    end_snp = max(snplist)
                else:
                    begin_snp = end_snp = -1

                scores = calculate_scores(group_test, lk_pred
                            , EST = 'lk'
                            , k = len(snplist), _k = k, SELECTOR = self.code
                            , MODELID = self.model_id, SIMID = simid, FOLD = k_fold
                            , BEGIN_SNP = begin_snp, END_SNP = end_snp
                            , **params)

                orig_scores = None
                if orig_pred is not None:
                        orig_scores = calculate_scores(group_test, orig_pred
                            , EST = self.code
                            , k = len(snplist), _k = k, SELECTOR = self.code
                            , MODELID = self.model_id, SIMID = simid, FOLD = k_fold
                            , BEGIN_SNP = begin_snp, END_SNP = end_snp
                            , **params)

                f_score = scores.loc[ scores['REG'] == 'MIN', 'MCC'].values[0]
                if f_score > best_score[0]:
                    best_score = (f_score, scores, orig_scores, snplist.tolist(), lk_pred)

            if best_score[0] < 0:
                continue

            results.append( best_score[1] )
            preds.append( (k, best_score[4]) )
            if best_score[2] is not None:
                results.append( best_score[2] )
            snps['%s/%d/%d/%d/%d'
                % (self.model_id, simid, k_fold, k, len(best_score[3]))] = best_score[3]

        log += [ '[I - {%d|%s}: %s]' % (simid, self.model_id, line)
                        for line in self.flush_log() ]

        if len(results) <= 0:
            return (pd.DataFrame(), snps, log, preds)
        return (pd.concat( results, sort=False ), snps, log, preds)


    def log(self, text):
        self.logs.append( text )

    def flush_log(self):
        log = self.logs
        self.logs = []
        return log

class AllSelector(BaseSelector):

    code = 'all'

    def select(self, haplotypes_train, groups_train, haplotypes_test, k=None):
        return (np.arange(len(haplotypes_train[0])), None, {})


class RandomSelector(BaseSelector):

    code = 'rand'

    def __init__(self, model_id, k, snpindex=None, iteration=1, seed=None):
        super().__init__(model_id, k, snpindex, iteration, seed)

    def select(self, haplotypes_train, groups_train, haplotypes_test, k):
        if self.snpindex is not None:
            return (np.random.choice(self.snpindex, k, replace=False), None, {})
        return (self.randomstate.randint(0, len(haplotypes_train[0]), k), None, {})


class FixSNPSelector(BaseSelector):

    code = 'fix'

    def __init__(self, model_id, k=None, snpindex=None, snpfile=None, iteration=1, seed=None):
        if snpindex is not None:
            self.L = snpindex = np.array(snpindex)
        elif snpfile:
            #self.L = np.array( [ int(x) for x in open(snpfile).read().split('\n')] )
            self.L = snpindex = np.loadtxt( snpfile, dtype=int)
        else:
            cerr('[Err - need snpindex or snpfile]')
        super().__init__(model_id, k, snpindex, iteration, seed)

    def select(self, haplotypes, groups, haplotest, k=None):
        return (self.L, None, {})


class DecisionTreeSelector(BaseSelector):

    code = 'dt'

    def __init__(self, model_id, k=None, max_depth=None, min_samples_leaf=2
            , snpindex=None, iteration=1, seed=None):

        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        super().__init__(model_id, k, snpindex, iteration, seed)


    def select(self, haplotypes_train, groups_train, haplotypes_test, k=None):

        if k <= 0:
            k = None

        classifier = DecisionTreeClassifier(class_weight='balanced', max_features=k
                , random_state = self.randomstate
                , max_depth = self.max_depth
                , min_samples_leaf = self.min_samples_leaf
                )
        classifier = classifier.fit(haplotypes_train, groups_train)
        features = classifier.tree_.feature
        d = { 'MAX_DEPTH': classifier.tree_.max_depth }

        return (np.unique(features[ features >= 0]), classifier.predict(haplotypes_test), d)


def parse_guide_tree( treefile ):

    pp.ParserElement.inlineLiteralsUsing(pp.Suppress)
    identity = pp.Word(pp.alphas + ' ')
    element = pp.Group(identity)
    parser = pp.OneOrMore( pp.nestedExpr(content=pp.delimitedList(element) + pp.Optional(','))
                      | pp.delimitedList(element))

    return parser.parseString( treefile.read() ).asList()


def flatten(lst):
    return sum( ([x] if not isinstance(x, list) else flatten(x)
             for x in lst), [] )


def traverse(tree, level=0):

    #print(tree)
    if len(tree) != 2:
        raise RuntimeError('[E - FATAL PROG ERR: misformat tree structure: child size=%d]' % len(tree))

    if isinstance(tree[0], list) and len(tree[0]) > 1:
        pop1 = flatten(tree[0])
        next1 = traverse(tree[0], level+1)
    else:
        pop1 = tree[0]
        next1 = []

    if isinstance(tree[1], list) and len(tree[1]) > 1:
        pop2 = flatten(tree[1])
        next2 = traverse(tree[1], level+1)
    else:
        pop2 = tree[1]
        next2 = []

    return [ (level, pop1, pop2) ] + next1 + next2


def count_allele(haplotypes):
    """ return numpy array of [ [c1, c2], [c1, c2], [c1, c2], .... ] based on haplotypes """

    n_alt = haplotypes.sum(axis=0)
    ac = np.empty( (len(n_alt), 2))
    ac[:,0] = len(haplotypes)*2 - n_alt
    ac[:,1] = n_alt
    return ac


class HierarchicalFSTSelector(BaseSelector):

    code = 'hfst'

    def __init__(self, model_id, k, guide_tree = None, min_fst = 0.9
                    , priority = None, max_leaf_snp = 0
                    , snpindex = None, iteration = 1, seed = None):
        super().__init__(model_id, k, snpindex, iteration, seed)

        if guide_tree is None:
            cexit('[E - HierarchicalFSTSelector requires guide tree]')
        self.guide_tree = guide_tree
        self.min_fst = min_fst
        self.priority = priority
        self.max_leaf_snp = max_leaf_snp


    def select(self, haplotypes, groups, haplotest, k=None):

        # we use k for redundancy parameters
        if k == 0 or k is None:
            k = 1

        candidate_L = []     # [ (pos, rank, no_actual_pops)]
        # we traverse through the tree
        for (level, pop1, pop2) in traverse(self.guide_tree):

            n_pops = len(pop1) + len(pop2)
            haplotypes1 = haplotypes[ np.isin(groups, pop1) ]
            haplotypes2 = haplotypes[ np.isin(groups, pop2) ]

            if len(haplotypes1) < 4:
                cerr('[I - insufficient population size for %s]' % pop1)
            if len(haplotypes2) < 4:
                cerr('[I - insufficient population size for %s]' % pop2)

            # convert haplotypes to allele counts
            ac1 = count_allele(haplotypes1)
            ac2 = count_allele(haplotypes2)

            # calculate highest FST
            FST = []
            num, den = allel.hudson_fst(ac1, ac2)

            # NOTE: the line below avoids warning (invalid value in true_divide)
            # when den == 0, which should be perfectly ok for FST calculation
            den[ den == 0 ] = -1
            fst = num/den

            # check for FST == 1.0
            ultimate_fst_pos = np.nonzero( fst == 1.0 )[0]
            if len(ultimate_fst_pos) > 0:
                self.log('FST: 1.0 at %s for pop %s <> %s' % (str(ultimate_fst_pos), pop1, pop2))

            if len(ultimate_fst_pos) > k and self.priority is not None:
                # get ultimate_fst based on priority

                ultimate_priority = self.priority[ ultimate_fst_pos ]
                sortidx = ultimate_fst_pos[ np.argsort( ultimate_priority ) ]

                #import IPython; IPython.embed()

            else:
                #fst[ np.isnan(fst) ] = 0
                sortidx = np.argsort( fst )

            # get highest FST
            #highest_fst_pos = sortidx[-(k+1):-1]
            #highest_fst_pos = list(reversed(sortidx))[:k]
            highest_fst_pos = sortidx[-k:]
            highest_fst_val = fst[ highest_fst_pos ]
            #self.log('highest FST: %5.4f at %d for pops %s <> %s' % (highest_fst_val, highest_fst_pos, pop1, pop2))
            if len(ultimate_fst_pos) > 0 and highest_fst_pos not in ultimate_fst_pos:
                pass
                #import IPython; IPython.embed()

            # check suitability of SNPs
            snplist, F = None, -1
            if highest_fst_val.max() < self.min_fst:

                if self.max_leaf_snp > k:

                    X_train =  np.append(haplotypes1, haplotypes2, axis=0)
                    y_train =  np.array( [1] * len(haplotypes1) + [2] * len(haplotypes2) )

                    best_iteration = (-1, None)
                    for i in range(k, self.max_leaf_snp):
                        features = sortidx[-(i+1):-1]

                        model = FixSNPSelector('dummy', snpindex=features)
                        lk_predictions, snplist, _, params = model.fit_and_predict(X_train
                                                            , y_train, X_train, len(features))
                        scores = calculate_scores(y_train,  lk_predictions)

                        F = scores.loc[ scores['REG'] == 'MIN', 'F'].values[0]
                        if best_iteration[0] < F:
                            best_iteration = (F, snplist)

                    snplist, F = best_iteration[1], best_iteration[0]

                snplist_2, F_2 = self.select_2(haplotypes1, haplotypes2)
                if F_2 > F:
                    snplist, F = snplist_2, F_2

                if snplist is not None:
                    self.log('F: %5.4f SNP: %d for pop %s <> %s => %s'
                        % (F, len(snplist), pop1, pop2, snplist))

                    for p in snplist:
                        candidate_L.append( (p, level, n_pops) )
                    continue

                # TODO: 2nd approach: find 2 SNPs with highest r^2(st) eg r^2 subpopulation vs r^2 total population


                # if snplist is None, just provide warning notice and skip this node!
                else:
                    self.log('low FST = %5.4f < %5.4f for %s vs %s; skipping'
                        % ( highest_fst_val.max(), self.min_fst, pop1, pop2))
                    continue

            # append to candidate_L
            for p in highest_fst_pos:
                candidate_L.append( (p, level, n_pops) )
            self.log('FST: %s SNP: %d for pop %s <> %s => %s'
                        % (str(highest_fst_val), len(highest_fst_pos), pop1, pop2, str(highest_fst_pos)))

        # process candidate_L
        L = np.unique( np.array( sorted( [ x[0] for x in candidate_L] ) ) )

        # return snp position
        return (L, None, {})


    def select_2(self, haplotypes1, haplotypes2):
        return None, -1


class HHFSTDTSelector(HierarchicalFSTSelector):

    code = 'hfst+dt'


    def __init__(self, model_id, k, guide_tree = None, min_fst = 0.9
                    , max_depth=None, min_samples_leaf=2
                    , priority = None, max_leaf_snp = 0
                    , snpindex = None, iteration = 1, seed = None):

        super().__init__(model_id, k, guide_tree, min_fst, priority, max_leaf_snp
                            , snpindex, iteration, seed)
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf

    def select_2(self, haplotypes1, haplotypes2):
        """ return (snplist, F):
            snplist - a list of SNP positions after further selection
            F = F score for these particular SNP set """

        X_train =  np.append(haplotypes1, haplotypes2, axis=0)
        y_train =  np.array( [1] * len(haplotypes1) + [2] * len(haplotypes2) )

        best_score = (-1, None, None, None)
        for i in range(3):

            classifier = DecisionTreeClassifier(class_weight='balanced'
                , random_state = self.randomstate
                , max_depth = self.max_depth
                , min_samples_leaf = self.min_samples_leaf)
            classifier = classifier.fit(X_train, y_train)
            features = classifier.tree_.feature

            # remove features with negative position and redundant
            features = np.unique(features[ features >= 0])

            model = FixSNPSelector('dummy', snpindex=features)
            lk_predictions, snplist, _, params = model.fit_and_predict(X_train
                                                , y_train, X_train, len(features))
            if lk_predictions is None:
                continue

            scores = calculate_scores(y_train,  lk_predictions)

            f_score = scores.loc[ scores['REG'] == 'MIN', 'F'].values[0]
            if f_score > best_score[0]:
                best_score = (f_score, scores, None, features.tolist())

        return best_score[3], best_score[0]




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
            positions: [ (chr, pos), (chr, pos), (<str>, <int>), ...]
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

def calculate_scores(target, prediction, **kwargs):

    assert len(target) == len(prediction)

    accs = []
    grp = {}
    for g,e in zip(target, prediction):
        # g: group, e: estimation/prediction
        if g not in grp:
            g_grp = grp[g] = [0, 0, 0, 0]   # TP, TN, FP, FN
        else:
            g_grp = grp[g]

        if e not in grp:
            e_grp = grp[e] = [0, 0, 0, 0]
        else:
            e_grp = grp[e]

        if g==e:
            g_grp[0] += 1   # TP
        else:
            g_grp[3] += 1   # FN should be g but instead e
            e_grp[2] += 1   # FP

    for g in grp:
        g_grp = grp[g]
        TP = g_grp[0]; FP = g_grp[2]; FN = g_grp[3]
        TN = len(prediction) - (TP+FP+FN)
        TP_FN = (TP + FN) if (TP+FN) > 0 else 1e-5
        TN_FP = (TN + FP) if (TN+FP) > 0 else 1e-5
        TP_FP = (TP + FP) if (TP+FP) > 0 else 1e-5
        TN_FN = (TN + FN) if (TN+FN) > 0 else 1e-5
        d = {
                'TN': TN, 'FP': FP, 'TP': TP, 'FN': FN,
                'ACC': (TP+TN)/(TP+TN+FP+FN),
                'TPR': TP/TP_FN,
                'TNR': TN/TN_FP,
                'PPV': TP/TP_FP,
                'NPV': TN/TN_FN,
                'MCC': (TP*TN - FP*FN)/np.sqrt(TP_FP * TP_FN * TN_FP * TN_FN),
                'REG': g
        }
        d.update( kwargs )
        d['BACC'] = (d['TPR'] + d['TNR'])/2
        d['F'] = 2 * d['PPV'] * d['TPR']/( d['PPV'] + d['TPR'] ) if (d['PPV'] + d['TPR']) > 0 else 0

        accs.append( d )

    # calculate median, mean, min
    scores = {}
    for score in [ 'ACC', 'TPR', 'TNR', 'PPV', 'NPV', 'BACC', 'F', 'MCC']:
        values = [ x[score] for x in accs ]
        scores[score] = ( np.median(values), np.mean(values), np.min(values)
                        , np.percentile(values, 25) )

    for i, g in enumerate(['MEDIAN', 'MEAN', 'MIN', 'Q1']):
        d = { 'REG': g, 'TN': -1, 'FP': -1, 'TP': -1, 'FN': -1 }
        d.update( kwargs )
        for score in [ 'ACC', 'TPR', 'TNR', 'PPV', 'NPV', 'BACC', 'F', 'MCC']:
            d[score] = scores[score][i]
        accs.append( d )

    return pd.DataFrame(accs)
