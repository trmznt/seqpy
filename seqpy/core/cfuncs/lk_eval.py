# Likelihood evaluator

import numpy as np

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



def calculate_accuracy(results, outfile):

    # for now, just print no of proportion with correct assignment
    props = []
    for paired_results in results:
        match = 0
        for g,e in paired_results:
            if g==e: match += 1

        props.append( match/len(paired_results) )

    print('Prop of correct assignment - avg = %f, med = %f. max = %f, min = %f' %
            ( np.mean(props), np.median(props), np.max(props), np.min(props) )
    )