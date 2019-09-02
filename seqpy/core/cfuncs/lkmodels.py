# Likelihood evaluator

import numpy as np
import pandas as pd

import itertools, time, math, os, pickle, datetime

import numpy as np, pandas as pd, pyparsing as pp
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import MultinomialNB, GaussianNB, ComplementNB, BernoulliNB
import allel

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import lkest


class ClassifierLK(object):

    classifier = lambda x: lkest.SNPLikelihoodEstimator()
    classifier_label = 'LK'

class ClassifierLR(object):

    classifier = lambda x: LogisticRegression()
    classifier_label = 'LR'


class ClassifierMNB(object):

    classifier = lambda x: MultinomialNB(fit_prior=False)
    classifier_label = 'MNB'


class ClassifierCNB(object):

    classifier = lambda x: ComplementNB()
    classifier_label = 'CNB'


class ClassifierGNB(object):

    classifier = lambda x: GaussianNB()
    classifier_label = 'GNB'

class ClassifierBNB(object):

	classifier = lambda x: BernoulliNB(binarize=1.0, fit_prior=False)
	classifier_label = 'BNB'


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
        if len(L) == len(X_train[0]):
            X_train_0 = X_train
            X_test_0 = X_test
        else:
            X_train_0 = X_train[:,L]
            X_test_0 = X_test[:,L]

        # create, train and predict using LikelihoodProfile

        lk_est = self.classifier()
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
                            , EST = self.classifier_label
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


class AllSelectorLR(AllSelector):

    classifier = lambda x: LogisticRegression()


class AllSelectorMNB(AllSelector):

    classifier = lambda x: MultinomialNB()


class AllSelectorGNB(AllSelector):

    classifier = lambda x: GaussianNB()


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


class FixSNPSelectorLK(FixSNPSelector, ClassifierLK):
	pass


class FixSNPSelectorLR(FixSNPSelector, ClassifierLR):
	pass


class FixSNPSelectorMNB(FixSNPSelector, ClassifierMNB):
	pass


class FixSNPSelectorCNB(FixSNPSelector, ClassifierCNB):
	pass


class FixSNPSelectorGNB(FixSNPSelector, ClassifierGNB):
	pass


class FixSNPSelectorBNB(FixSNPSelector, ClassifierBNB):
	pass


class DecisionTreeSelector(BaseSelector):

    code = 'DT'

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


class RandomForestSelector(DecisionTreeSelector):

    def select(self, haplotypes_train, groups_train, haplotypes_test, k=None):

        if k <= 0:
            k = None

        classifier = RandomForestClassifier(class_weight='balanced', max_features=k
                , random_state = self.randomstate
                , max_depth = self.max_depth
                , min_samples_leaf = self.min_samples_leaf
                , n_estimator = 50
                )
        classifier = classifier.fit(haplotypes_train, groups_train)
        importances = classifier.importances_        
        features = np.where( importances > 0.15 )
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

                        F = scores.loc[ scores['REG'] == 'MIN', 'MCC'].values[0]
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

            f_score = scores.loc[ scores['REG'] == 'MIN', 'MCC'].values[0]
            if f_score > best_score[0]:
                best_score = (f_score, scores, None, features.tolist())

        return best_score[3], best_score[0]



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
