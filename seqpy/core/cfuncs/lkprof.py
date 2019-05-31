# Likelihood evaluator

import numpy as np
import pandas as pd

import itertools, time, math, os, pickle, datetime
from multiprocessing import Pool, RawArray

import numpy as np, pandas as pd, pyparsing as pp
from sklearn.model_selection import StratifiedKFold
from sklearn.tree import DecisionTreeClassifier
import allel

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import naltparser
from seqpy.core.cfuncs import lkprof


class BaseSelector(object):

    code = 'base'

    def __init__(self, model_id, k, snpindex=None, iteration=1, seed=None):
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
            (lk_prediction, actual_k, model_prediction)
        """

        # select SNPs based on model
        L, orig_predictions, params = self.select(X_train, y_train, X_test, k)
        X_train_0 = X_train[:,L]
        X_test_0 = X_test[:,L]

        # create, train and predict using LikelihoodProfile

        lkp = lkprof.LikelihoodProfile()
        lkp.fit(X_train_0, y_train)

        return lkp.predict(X_test_0), L, orig_predictions, params


    def score(self, genotype_train, group_train, genotype_test, group_test, simid, k_fold):
        """ return a dataframe containing scores and dict of snps """

        results = []
        snps = {}
        log = []

        for k in self.k_list:

            # best score containe (F, score_dataframe, orig_score_dataframe, snplist)
            best_score = (-1, None, None, None)

            for i in range(self.iteration):

                lk_pred, snplist, orig_pred, params = self.fit_and_predict( genotype_train
                        , group_train, genotype_test, k)

                scores = lkprof.calculate_scores(group_test, lk_pred
                            , EST = 'lk'
                            , k = len(snplist), _k = k, SELECTOR = self.code
                            , MODELID = self.model_id, SIMID = simid, FOLD = k_fold
                            , **params)

                orig_scores = None
                if orig_pred is not None:
                        orig_scores = lkprof.calculate_scores(group_test, orig_pred
                            , EST = self.code
                            , k = len(snplist), _k = k, SELECTOR = self.code
                            , MODELID = self.model_id, SIMID = simid, FOLD = k_fold
                            , **params)

                f_score = scores.loc[ scores['REG'] == 'MIN', 'F'].values[0]
                if f_score > best_score[0]:
                    best_score = (f_score, scores, orig_scores, snplist.tolist())

            results.append( best_score[1] )
            if best_score[2] is not None:
                results.append( best_score[2] )
            snps['%s/%d/%d/%d/%d'
                % (self.model_id, simid, k_fold, k, len(best_score[3]))] = best_score[3]

        log += [ '[I - {%d|%s}: %s]' % (simid, self.model_id, line)
                        for line in self.flush_log() ]

        return (pd.concat( results, sort=False ), snps, log)


    def log(self, text):
        self.logs.append( text )

    def flush_log(self):
        log = self.logs
        self.logs = []
        return log


class RandomSelector(BaseSelector):

    code = 'rand'

    def __init__(self, model_id, k, snpindex=None, iteration=1, seed=None):
        super().__init__(model_id, k, snpindex, iteration, seed)

    def select(self, haplotypes_train, groups_train, haplotypes_test, k):
        if self.snpindex:
            return (np.random.choice(snpindex, k, replace=False), None, {})
        return (self.randomstate.randint(0, len(haplotypes_train[0]), k), None, {})


class FixSNPSelector(BaseSelector):

    code = 'fix'

    def __init__(self, model_id, k=None, snpindex=None, snpfile=None, iteration=1, seed=None):
        if snpindex is not None:
            self.L = snpindex
        elif snpfile:
            #self.L = np.array( [ int(x) for x in open(snpfile).read().split('\n')] )
            self.L = snpindex = np.loadtxt( snpfile, dtype=int)
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
                    , snpindex = None, iteration = 1, seed = None):
        super().__init__(model_id, k, snpindex, iteration, seed)

        if guide_tree is None:
            cexit('[E - HierarchicalFSTSelector requires guide tree]')
        self.guide_tree = guide_tree
        self.min_fst = min_fst


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

            # NOTE: the line below might produce warning (invalid value in true_divide)
            # if den == 0, which should be perfectly ok for FST calculation
            fst = num/den

            fst[ np.isnan(fst) ] = 0
            sortidx = np.argsort( fst )

            # get highest FST
            highest_fst_pos = sortidx[-(k+1):-1]
            highest_fst_val = fst[ highest_fst_pos ]
            #cerr('[I - highest FST: %5.4f at %d for pops %s and %s' % (highest_fst_val, highest_fst_pos, pop1, pop2))

            # check suitability of SNPs
            if highest_fst_val.max() < self.min_fst:

                snplist, F = self.select_2(haplotypes1, haplotypes2)
                if snplist:
                    self.log('F: %5.4f SNP: %d for pop %s <> %s => %s'
                        % (F, len(snplist), pop1, pop2, snplist))

                    for p in snplist:
                        candidate_L.append( (p, level, n_pops) )
                    continue

                # TODO: 2nd approach: find 2 SNPs with highest r^2(st) eg r^2 subpopulation vs r^2 total population


                # if snplist is None, just provide warning notice !
                else:
                    self.log('low FST = %5.4f for %s vs %s' % ( highest_fst_val.max(), pop1, pop2))

            # append to candidate_L
            for p in highest_fst_pos:
                candidate_L.append( (p, level, n_pops) )

        # process candidate_L
        L = np.unique( np.array( sorted( [ x[0] for x in candidate_L] ) ) )

        # return snp position
        return (L, None, {})


    def select_2(self, haplotypes1, haplotypes2):
        return None, None


class HHFSTDTSelector(HierarchicalFSTSelector):

    code = 'hfst+dt'


    def __init__(self, model_id, k, guide_tree = None, min_fst = 0.9
                    , max_depth=None, min_samples_leaf=2
                    , snpindex = None, iteration = 1, seed = None):

        super().__init__(model_id, k, guide_tree, min_fst, snpindex, iteration, seed)
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
                , random_state = self.randomstate, min_samples_leaf = self.min_samples_leaf)
            classifier = classifier.fit(X_train, y_train)
            features = classifier.tree_.feature

            # remove features with negative position and redundant
            features = np.unique(features[ features >= 0])

            model = FixSNPSelector('dummy', snpindex=features)
            lk_predictions, snplist, _, params = model.fit_and_predict(X_train
                                                , y_train, X_train, len(features))
            scores = lkprof.calculate_scores(y_train,  lk_predictions)

            f_score = scores.loc[ scores['REG'] == 'MIN', 'F'].values[0]
            if f_score > best_score[0]:
                best_score = (f_score, scores, None, features.tolist())

        return best_score[3], best_score[0]


def prepare_stratified_samples(haplotypes, group_keys, k_fold, haplotype_func=None):
    """ check the suitability of sample sets and modify haplotypes and group_keys properly """

    groups = []
    for group_key, count in zip( *np.unique(group_keys, return_counts=True)):
        # we make sure that every group has at least 2 * k_fold member
        if count < k_fold * 2:
            groups.append( (group_key, math.ceil(k_fold*2 / count)) )

    if len(groups) == 0:
        # nothing to modify
        return (haplotypes, group_keys)

    cerr('[I - prepare_stratified_sample() replicated group: %s]'
            % ' '.join( x[0] for x in groups ))
    #import IPython; IPython.embed()
    new_haplotypes = [ haplotypes ]
    new_group_keys = [ group_keys ]
    for group_key, m_factor in groups:
        indexes = np.where( group_keys == group_key )
        for i in range(m_factor):
            new_haplotypes.append( haplotypes[indexes] )
            new_group_keys.append( group_keys[indexes] )

    haplotypes = np.concatenate( new_haplotypes, axis=0 )
    group_keys = np.concatenate( new_group_keys, axis=0 )

    return (haplotypes, group_keys)


def cross_validate_worker( args ):

    pid = os.getpid()
    y, fold, simid = args
    models = var_dict['models']

    cerr('[I - pid %d: validator_worker() started]' % pid)

    np.random.seed(simid % pid)

    # reseed all models
    for model in models: model.reseed(np.random.randint(1e8))

    # obtain haplotype matrix
    if var_dict['X_shape'] == None:
        X = var_dict['X']
    else:
        X = np.frombuffer( var_dict['X'], dtype=np.int8 ).reshape( var_dict['X_shape'])

    results = []
    snps = {}
    k_fold = -1
    log = []

    if fold <= 0:

        for m in models:

            scores, snplist, mlog = m.score(X, y, X, y, simid, k_fold)

            results.append( scores )
            snps.update( snplist )
            log += mlog

        return (simid, pd.concat(results, sort=False), snps, log)

    X, y = prepare_stratified_samples(X, y, fold)

    skf = StratifiedKFold(n_splits = fold, shuffle=True, random_state = np.random.randint(1e8))

    for train_index, test_index in skf.split(X, y):

        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        k_fold += 1

        for m in models:

            scores, snplist, mlog = m.score(X_train, y_train, X_test, y_test, simid, k_fold)

            results.append( scores )
            snps.update( snplist )
            log += mlog

    return (simid, pd.concat(results, sort=False), snps, log)


# global variable for multiprocessing
var_dict = {}

def init_worker(X, X_shape, models):
    var_dict['X'] = X
    var_dict['X_shape'] = X_shape
    var_dict['models'] = models


def cross_validate( models, haplotypes, group_keys, repeats, fold
        , outfile, procs=1, outsnp=None, logfile=None):

    start_time = time.monotonic()
    cerr('[I - cross-validation for %d model(s) with repeats: %d %d-fold]'
        % (len(models), repeats, fold))

    logf = None
    if logfile:
        logf = open(logfile, 'w')

    group_keys = np.array(group_keys) if type(group_keys) != np.ndarray else group_keys

    seed = np.random.randint(1e7)
    arguments = [ (group_keys, fold, seed+n) for n in range(repeats) ]


    if procs > 1:
        # perform multiprocessing

        # create a shared-memory for numpy array
        cerr('[I - preparing for shared-memory numpy array]')
        X_shape = haplotypes.shape
        X = RawArray('b', X_shape[0]*X_shape[1])
        X_np = np.frombuffer(X, dtype=np.int8).reshape(X_shape)

        # fill share-memory with haplotypes
        np.copyto(X_np, haplotypes)

        with Pool(procs, initializer=init_worker, initargs=(X, X_shape, models)) as pool:
            c = 0
            for (n, result, snps, log) in pool.imap_unordered(cross_validate_worker, arguments):
                c += 1
                cerr('[I - receiving result from repeat #%d (%d/%d) with %d results]'
                        % (n-seed+1, c, repeats, len(result)))
                # write to temporary files
                if outfile:
                    with open('%s.%d' % (outfile, n), 'wb') as fout:
                        pickle.dump(result, fout, pickle.HIGHEST_PROTOCOL)
                if outsnp:
                    with open('%s.%d' % (outsnp, n), 'wb') as fout:
                        pickle.dump(snps, fout, pickle.HIGHEST_PROTOCOL)

                # write to log
                if logf and log:
                    logf.write( '\n'.join( log ) )
                    logf.write( '\n' )

    else:

        init_worker( haplotypes, None, models )
        c = 0
        for (n, result, snps, log) in map(cross_validate_worker, arguments ):
            c += 1
            cerr('[I - receiving result from repeat #%d (%d/%d) with %d results]'
                    % (n-seed+1, c, repeats, len(result)))

            # write to temporary files
            if outfile:
                with open('%s.%d' % (outfile, n), 'wb') as fout:
                    pickle.dump(result, fout, pickle.HIGHEST_PROTOCOL)
            if outsnp:
                with open('%s.%d' % (outsnp, n), 'wb') as fout:
                    pickle.dump(snps, fout, pickle.HIGHEST_PROTOCOL)

            # write to log
            if logf and log:
                logf.write( '\n'.join( log ) )
                logf.write( '\n' )

    cerr('[I - combining output files]')
    results = []
    snp_tables = {}

    for n in range(repeats):

        if outfile:
            filename = '%s.%d' % (outfile, n+seed)
            with open(filename, 'rb') as fin:
                results.append( pickle.load(fin) )
            os.remove(filename)

        if outsnp:
            filename = '%s.%d' % (outsnp, n+seed)
            with open(filename, 'rb') as fin:
                snp_tables.update( pickle.load(fin) )
            os.remove(filename)

    if outfile:
        df = pd.concat( results )
        df.to_csv(outfile, sep='\t', index=False)
        cerr('[I - writing scores to %s]' % outfile)

    if outsnp:
        pickle.dump(snp_tables, open(outsnp, 'wb'))
        cerr('[I - writing SNP table to %s]' % outsnp )

    cerr('[I - cross-validate() finished in %6.2f minute(s) at %s]'
            % ((time.monotonic() - start_time)/60, datetime.datetime.now()))


class LikelihoodProfile(object):


    def __init__(self, matrices=None):
        self.matrices = matrices


    def fit(self, X_train, y_train):
        """ create a set of matrices for each group """

        # create empty matrices
        M = {}
        l = len(X_train[0])
        for y in set(y_train):
            M[y] = np.full(shape=(l, 2), fill_value=0.001)

        # filling-up matrices
        # we are iterating row wise
        for i, y in enumerate(y_train):
            n_alt = X_train[i]
            M_y = M[y]

            M_y[n_alt == 0, 0] += 2
            M_y[n_alt == 1] += 1
            M_y[n_alt == 2, 1] += 2

        # convert M into fractions
        for y in M:
            M_y = M[y]
            total = np.sum(M_y, axis=1)
            M[y] = np.log(M_y / total[:, None])

        self.matrices = M

        return self


    def predict(self, X_test, with_confidences=False):

        predictions = []

        for n_alt_data in X_test:
            lkscore = self.calculate_likelihoods(n_alt_data)
            lkscore.sort()
            predictions.append( lkscore[-1][1] )

        return predictions


    def predict_with_confidence(self, X_test):
        """ not: condifende score need to be rethought!! """

        predictions = []
        confidences = []

        for n_alt_data in X_test:
            lkscore = self.calculate_likelihoods(n_alt_data)
            lkscore.sort()
            predictions.append( lks[-1][1] )
            confidences.append( lkscore[-1][0]/lkscore[-2][0] )

        return predictions, confidences


    def predict_with_proba(self, X_test):

        predictions = []
        probabilities = []

        for n_alt_data in X_test:
            lkscore = self.calculate_likelihoods(n_alt_data)
            lkscore.sort()
            predictions.append( lks[-1][1] )

            lkhoods = [ math.e*x[0] for x in lkscore ]
            lk_ratio = lkhoods[0] / sum(lkhoods[1:])
            proba = lk_ratio/(1+lk_ratio)
            probabilities.append( proba )

        return predictions, probabilities


    def calculate_likelihoods(self, n_alt_data):
        """ calculate log likelihood of each group """

        R = []
        for g in self.matrices:
            R.append(
                ( self.calculate_likelihood(self.matrices[g], n_alt_data), g)
            )

        return R


    def calculate_likelihood(self, M, n_alt_data):
        n_alt_0 = 2 * np.sum(M[np.where(n_alt_data == 0), 0])

        # optimize for scarce amount of 1 alt data
        index = np.where(n_alt_data == 1)
        if M[index].size:
            n_alt_1 = np.sum(M[index])
        else:
            n_alt_1 = 0

        n_alt_2 = 2 * np.sum(M[np.where(n_alt_data == 2), 1])

        log_lk = n_alt_0 + n_alt_1 + n_alt_2

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
        d = {
                'ACC': (TP+TN)/(TP+TN+FP+FN),
                'TPR': TP/(TP + FN) if (TP+FN) > 0 else 1e-5,
                'TNR': TN/(TN + FP) if (TN+FP) > 0 else 1e-5,
                'PPV': TP/(TP+FP) if (TP+FP) > 0 else 1e-5,
                'NPV': TN/(TN+FN) if (TN+FN) > 0 else 1e-5,
                'REG': g
        }
        d.update( kwargs )
        d['BACC'] = (d['TPR'] + d['TNR'])/2
        d['F'] = 2 * d['PPV'] * d['TPR']/( d['PPV'] + d['TPR'] ) if (d['PPV'] + d['TPR']) > 0 else 0

        accs.append( d )

    # calculate median, mean, min
    scores = {}
    for score in [ 'ACC', 'TPR', 'TNR', 'PPV', 'NPV', 'BACC', 'F']:
        values = [ x[score] for x in accs ]
        scores[score] = ( np.median(values), np.mean(values), np.min(values) )

    for i, g in enumerate(['MEDIAN', 'MEAN', 'MIN']):
        d = { 'REG': g }
        d.update( kwargs )
        for score in [ 'ACC', 'TPR', 'TNR', 'PPV', 'NPV', 'BACC', 'F']:
            d[score] = scores[score][i]
        accs.append( d )

    return pd.DataFrame(accs)


    #    props.append( match/len(paired_results) )

    #print('Prop of correct assignment - avg = %f, med = %f. max = %f, min = %f' %
    #        ( np.mean(props), np.median(props), np.max(props), np.min(props) )
    #)
