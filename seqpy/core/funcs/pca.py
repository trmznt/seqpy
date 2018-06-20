#
# this is ripped from matplotlib, so that users with matplotlib < 1.2 can
# still use PCA module

import numpy as np

class PCA(object):

    def __init__(self, a, standardize=True):
        """
        compute the SVD of a and store data for PCA.  Use project to
        project the data onto a reduced set of dimensions

        Inputs:

          *a*: a numobservations x numdims array
          *standardize*: True if input data are to be standardized. If False,
          only centering will be carried out.

        Attrs:

          *a* a centered unit sigma version of input a

          *numrows*, *numcols*: the dimensions of a

          *mu*: a numdims array of means of a. This is the vector that points
          to the origin of PCA space.

          *sigma*: a numdims array of standard deviation of a

          *fracs*: the proportion of variance of each of the principal
          components

          *s*: the actual eigenvalues of the decomposition

          *Wt*: the weight vector for projecting a numdims point or array into
          PCA space

          *Y*: a projected into PCA space


        The factor loadings are in the Wt factor, i.e., the factor
        loadings for the 1st principal component are given by Wt[0].
        This row is also the 1st eigenvector.

        """
        n, m = a.shape
        if n < m:
            raise RuntimeError('we assume data in a is organized with '
                               'numrows>numcols')

        self.numrows, self.numcols = n, m
        self.mu = a.mean(axis=0)
        self.sigma = a.std(axis=0)
        self.standardize = standardize

        a = self.center(a)

        self.a = a

        U, s, Vh = np.linalg.svd(a, full_matrices=False)

        # Note: .H indicates the conjugate transposed / Hermitian.

        # The SVD is commonly written as a = U s V.H.
        # If U is a unitary matrix, it means that it satisfies U.H = inv(U).

        # The rows of Vh are the eigenvectors of a.H a.
        # The columns of U are the eigenvectors of a a.H.
        # For row i in Vh and column i in U, the corresponding eigenvalue is
        # s[i]**2.

        self.Wt = Vh

        # save the transposed coordinates
        Y = np.dot(Vh, a.T).T
        self.Y = Y

        # save the eigenvalues
        self.s = s**2

        # and now the contribution of the individual components
        vars = self.s/float(len(s))
        self.fracs = vars/vars.sum()

    def project(self, x, minfrac=0.):
        '''
        project x onto the principle axes, dropping any axes where fraction
        of variance<minfrac
        '''
        x = np.asarray(x)

        ndims = len(x.shape)

        if (x.shape[-1] != self.numcols):
            raise ValueError('Expected an array with dims[-1]==%d' %
                             self.numcols)

        Y = np.dot(self.Wt, self.center(x).T).T
        mask = self.fracs >= minfrac
        if ndims == 2:
            Yreduced = Y[:, mask]
        else:
            Yreduced = Y[mask]
        return Yreduced

    def center(self, x):
        '''
        center and optionally standardize the data using the mean and sigma
        from training set a
        '''
        if self.standardize:
            return (x - self.mu)/self.sigma
        else:
            return (x - self.mu)


def jitters( data ):
    """ returning normal distribution noise for data jittering"""
    # need to calculate s based on the smallest and longest distance of the points
    d_min = d_max = np.abs(data[0] - data[1])
    for i in range(len(data)):
        for j in range(len(data)):
            d = np.abs( data[i] - data[j] )
            if d > 0:
                if d > d_max:
                    d_max = d
                elif d < d_min:
                    d_min = d

    s = min(d_min / d_max * len(data) / 2, d_max / 75)
    return s * np.random.randn( len(data) )
