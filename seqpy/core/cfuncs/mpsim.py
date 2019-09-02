# single/multi-procs simulator

import itertools, time, os, pickle, datetime, math
from multiprocessing import Pool, RawArray
from collections import defaultdict

import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold

from seqpy import cout, cerr, cexit



class Segment(object):
    """ Segment instance will be passed to other process during multi-process,
        so this class must contain as minimum data as possible, without
        referencing other classes/instances that might needed to be passed
        as well.
    """

    def __init__(self, segment, begin, end):
        self.segment = segment
        self.begin = begin
        self.end = end

    def get_haplotypes(self, haplotypes):
        #import IPython; IPython.embed()
        return haplotypes[:, self.begin:self.end + 1]

    def add_annotation(self, dataframe):
        dataframe['SEGMENT'] = self.segment
        return dataframe

    def __repr__(self):
        return '<Segment: %s [%d-%d]>' % (self.segment, self.begin, self.end)


def create_gene_segments(region, group_keys, start_count=None):

    current_segment = None
    counter = itertools.count(start_count or np.random.randint(1e8))
    for idx, posline in enumerate(region.P):
        if type(posline[4]) is not str:
            continue
        if current_segment is None:
            current_segment = Segment(posline[4], idx, idx)
            continue
        if current_segment.segment != posline[4]:
            if current_segment.end - current_segment.begin > 3:
                yield next(counter), current_segment, group_keys
            current_segment = Segment(posline[4], idx, idx)
            continue
        current_segment.end = idx
    if current_segment.end - current_segment.begin > 3:
        yield next(counter), current_segment, group_keys


def create_window_segments(region, group_keys, window_size = 250, min_snp = 3, start_count=None):

    curr_segment = None
    counter = itertools.count(start_count or np.random.randint(1e8))
    cur_chrom = ''
    end_idx = last_end_idx = 0
    max_idx = len(region.P) - 1
    for idx, posline in enumerate(region.P):
        #print(idx, posline)
        chrom, pos = posline[0], posline[1]

        # move end_idx until length is longer than window_size
        end_posline = region.P[end_idx]
        end_chrom, end_pos = end_posline[0], end_posline[1]

        while end_idx < max_idx and end_pos - pos <= window_size and end_chrom == chrom:
            end_idx += 1
            end_posline = region.P[end_idx]
            end_chrom, end_pos = end_posline[0], end_posline[1]

        if last_end_idx >= end_idx - 1:
            continue

        if end_idx - idx < min_snp:
            continue

        last_end_idx = end_idx - 1
        end_posline = region.P[last_end_idx]
        yield   ( next(counter)
                , Segment('%s:%d/%d/%d' % (chrom, pos, end_posline[1] - pos, last_end_idx - idx+1), idx, last_end_idx)
                , group_keys
                )


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


def scan_segment_worker( args ):

    pid = os.getpid()
    simid, segment, y = args
    models = var_dict['models']

    cerr('[I - pid %d: scan_segment__worker() started for %s]' % (pid, str(segment)) )

    np.random.seed(simid % pid)

    # reseed all models
    for model in models: model.reseed(np.random.randint(1e8))

    # obtain haplotype matrix
    if var_dict['X_shape'] == None:
        X = var_dict['X']
    else:
        X = np.frombuffer( var_dict['X'], dtype=np.int8 ).reshape( var_dict['X_shape'])

    X_ = segment.get_haplotypes(X)
    results = []
    snps = {}
    log = []

    for m in models:

        try:
            cerr('[I - pid %d: scoring model %s]' % (pid, m.model_id))
            scores, snplist, mlog, preds = m.score(X_, y, X_, y, simid, -1)

        except:
            cerr('ERR - pid %d model %s segment %s' % (pid, m.model_id, str(segment)))
            raise

        results.append( scores )
        snps.update( snplist )
        log += mlog

    return (simid, segment.add_annotation(pd.concat(results, sort=False)), snps, log, [])


def cross_validate_worker( args ):

    pid = os.getpid()
    y, fold, simid = args
    models = var_dict['models']

    cerr('[I - pid %d: cross_validate_worker() started]' % pid)

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
    predictions = []

    if fold <= 0:

        # using whole data set, no cross-validation

        for m in models:

            cerr('[I - pid %d: scoring model %s]' % (pid, m.model_id))
            scores, snplist, mlog, preds = m.score(X, y, X, y, simid, k_fold)

            results.append( scores )
            snps.update( snplist )
            log += mlog
            predictions.append( (m.model_id, preds) )

        return (simid, pd.concat(results, sort=False), snps, log, predictions)

    orig_y_size = len(y)
    X, y = prepare_stratified_samples(X, y, fold)

    skf = StratifiedKFold(n_splits = fold, shuffle=True, random_state = np.random.randint(1e8))

    for train_index, test_index in skf.split(X, y):

        # using stratified k-fold cross validation data set

        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        k_fold += 1

        for m in models:

            cerr('[I - pid %d: scoring model %s]' % (pid, m.model_id))
            scores, snplist, mlog, preds = m.score(X_train, y_train, X_test, y_test, simid, k_fold)

            results.append( scores )
            snps.update( snplist )
            log += mlog

            predictions.append( (m.model_id, test_index, preds) )

        # gather prediction results from all fold and remove the duplicated samples

    predictions_by_models = {}
    for model_id, indexes, preds in predictions:
        try:
            model_preds = predictions_by_models[model_id]
        except KeyError:
            model_preds = predictions_by_models[model_id] = {}
        for (k, values) in preds:
            try:
                pred_results = model_preds[k]
            except KeyError:
                pred_results = model_preds[k] = [''] * orig_y_size
            for i, p in zip(indexes, values):
                if i < orig_y_size:
                    pred_results[i] = p

    aggregate_predictions = []
    for k in predictions_by_models:
        aggregate_predictions.append( (k, list(predictions_by_models[k].items())) )

    #import IPython; IPython.embed()
    return (simid, pd.concat(results, sort=False), snps, log, aggregate_predictions)


# global variable for multiprocessing
var_dict = {}

def init_worker(X, X_shape, models):
    var_dict['X'] = X
    var_dict['X_shape'] = X_shape
    var_dict['models'] = models


def run_worker(models, haplotypes, arguments, worker_func, procs, outfile, outsnp, logfile, outpred=None):

    logf = None
    if logfile:
        logf = open(logfile, 'w')

    simids = []
    pred_results = []

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
            for (n, result, snps, log, m_preds) in pool.imap_unordered(worker_func, arguments):
                c += 1
                cerr('[I - receiving result from simid #%d (%d) with %d results]'
                        % (n, c, len(result)))
                simids.append(n)

                # write to temporary files
                if outfile:
                    with open('%s.%d' % (outfile, n), 'wb') as fout:
                        pickle.dump(result, fout, pickle.HIGHEST_PROTOCOL)
                if outsnp:
                    with open('%s.%d' % (outsnp, n), 'wb') as fout:
                        pickle.dump(snps, fout, pickle.HIGHEST_PROTOCOL)
                if outpred:
                    with open('%s.%d' % (outpred, n), 'wb') as fout:
                        pickle.dump(m_preds, fout, pickle.HIGHEST_PROTOCOL)

                # write to log
                if logf and log:
                    logf.write( '\n'.join( log ) )
                    logf.write( '\n' )

    else:

        init_worker( haplotypes, None, models )
        c = 0
        for (n, result, snps, log, m_preds) in map(worker_func, arguments ):
            c += 1
            cerr('[I - receiving result from simid #%d (%d) with %d results]'
                    % (n, c, len(result)))
            simids.append(n)

            if outfile:
                with open('%s.%d' % (outfile, n), 'wb') as fout:
                    pickle.dump(result, fout, pickle.HIGHEST_PROTOCOL)
            if outsnp:
                with open('%s.%d' % (outsnp, n), 'wb') as fout:
                    pickle.dump(snps, fout, pickle.HIGHEST_PROTOCOL)
            if outpred:
                with open('%s.%d' % (outpred, n), 'wb') as fout:
                    pickle.dump(m_preds, fout, pickle.HIGHEST_PROTOCOL)

            # write to log
            if logf and log:
                logf.write( '\n'.join( log ) )
                logf.write( '\n' )



    cerr('[I - combining output files]')
    results = []
    snp_tables = {}
    predictions = []

    for n in simids:

        if outfile:
            filename = '%s.%d' % (outfile, n)
            with open(filename, 'rb') as fin:
                results.append( pickle.load(fin) )
            os.remove(filename)

        if outsnp:
            filename = '%s.%d' % (outsnp, n)
            with open(filename, 'rb') as fin:
                snp_tables.update( pickle.load(fin) )
            os.remove(filename)

        if outpred:
            filename = '%s.%d' % (outpred, n)
            with open(filename, 'rb') as fin:
                predictions.extend( pickle.load(fin) )
            os.remove(filename)


    if outfile:
        df = pd.concat( results, sort=False )
        df.to_csv(outfile, sep='\t', index=False)
        cerr('[I - writing scores to %s]' % outfile)

    if outsnp:
        pickle.dump(snp_tables, open(outsnp, 'wb'))
        cerr('[I - writing SNP table to %s]' % outsnp )

    if outpred:
        all_preds = consolidate_predictions( predictions )
        pickle.dump(all_preds, open(outpred, 'wb'))
        cerr('[I - writing aggregate predictions to %s]' % outpred)


def scan_segment(models, haplotypes, group_keys, arguments
        , outfile, outsnp=None, logfile=None, procs=1):
    """ distribute the arguments over multi processes
        models: list of model instance, with suitable score method
        haplotypes: list of haploytpe
        arguments:
    """

    start_time = time.monotonic()
    cerr('[I - scan_segment() for %d model(s)]'
        % (len(models)))

    worker_func = scan_segment_worker
    run_worker(models, haplotypes, arguments, worker_func, procs, outfile, outsnp, logfile)

    cerr('[I - scan_segment() finished in %6.2f minute(s) at %s]'
            % ((time.monotonic() - start_time)/60, datetime.datetime.now()))


def cross_validate(models, haplotypes, group_keys, repeats, fold
        , outfile, outsnp=None, logfile=None, outpred=None, procs=1):
    """ distribute the repeats over multi process
    """

    start_time = time.monotonic()
    cerr('[I - cross_validate() for %d model(s)]'
        % (len(models)))

    seed = np.random.randint(1e7)
    group_keys = np.array(group_keys) if type(group_keys) != np.ndarray else group_keys

    arguments = [ (group_keys, fold, seed+n) for n in range(repeats) ]

    worker_func = cross_validate_worker
    run_worker(models, haplotypes, arguments, worker_func, procs, outfile, outsnp, logfile, outpred)

    cerr('[I - cross_validate() finished in %6.2f minute(s) at %s]'
            % ((time.monotonic() - start_time)/60, datetime.datetime.now()))


def consolidate_predictions(m_preds):

    # m_preds = [('ALL',  [(0, [ group1, group2, ...]), ... ]) ]

    model_predictions = {}

    for model_id, predictions in m_preds:
        try:
            current_predictions = model_predictions[model_id]
        except KeyError:
            current_predictions = model_predictions[model_id] = {}

        for k, result_list in predictions:
            try:
                current_predictions[k].append(result_list)
            except KeyError:
                current_predictions[k] = [ result_list ]

    return model_predictions

