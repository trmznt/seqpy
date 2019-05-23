# Likelihood evaluator

import numpy as np
import pandas as pd


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


class LikelihoodEvaluator(object):

    def __init__(self, X_train, y_train, X_test, y_test):
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test
        self.groups = set(y_train)


    def estimate(self, n_alt_data):

        lks = self.calculate_likelihoods(n_alt_data)
        lks.sort()
        return lks[-1][1]


    def evaluate(self):
        """ return the accuracy of the region """

        estimates = []

        for grp, n_alt_data in zip(self.y_test, self.X_test):
            est_group = self.estimate(n_alt_data)
            estimates.append((grp, est_group))

        return estimates


    def prepare_matrix(self):
        """ create a set of matrices for each group """

        # create empty matrices
        M = {}
        l = len(self.X_train[0])
        for group in self.groups:
            M[group] = np.full(shape=(l, 2), fill_value=0.01)

        # filling-up matrices
        # we are iterating row wise
        for i, group in enumerate(self.y_train):
            n_alt = self.X_train[i]
            M_group = M[group]

            M_group[n_alt == 0, 0] += 2
            M_group[n_alt == 1] += 1
            M_group[n_alt == 2, 1] += 2

        # convert M into fractions
        for group in M:
            M_group = M[group]
            total = np.sum(M_group, axis=1)
            M[group] = np.log(M_group / total[:, None])

        self.matrices = M
        return M


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
