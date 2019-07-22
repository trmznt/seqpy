# naltparser

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

from seqpy.core.bioio import grpparser

import numpy as np
import pandas as pd
import attr
import io
import itertools
import functools

# requires scikit-allel
try:
    import allel
except:
    cexit('ERR: require properly installed scikit-allel!')



class Region(object):

    def __init__(self, name, P=None, M=None):
        self.name = name    # name of region
        self.P = P or []        # position

        self.M = M or []        # n_alt matrix (no of alternate allele)
        # self.M is structured as
        # [ [ snp1_smaple1 snp1_sample2 snp1_sample3 ...]
        #    [ snp2_sample1 snp2_sample2 snp2_sample3 ...]
        # ]

        self.H = None            # haploytpes as 2D numpy array of short
        # self.H is structured as
        # [    sample1_snp1 sample1_snp2 sample1_snp3 ...
        #    sample2_snp1 sample2_snp2 sample2_snp3 ...
        # ]

    def append(self, posinfo, n_alt):
        self.P.append(posinfo)
        self.M.append(n_alt)

    def samples(self):
        return list(self.df_M)

    def haplotypes(self):
        if self.H:
            return self.H

        self.H = np.transpose( np.array(self.M) )
        return self.H

    def ralt_to_nalt(self, hetratio=0.25):

        from genoutils import ralt_to_nalt

        if hetratio < 0:
            n_mdp = np.array([], dtype=np.intc)
        for i in range(len(self.M)):
            self.M[i] = ralt_to_nalt(self.M[i], n_mdp, hetratio)

    def parse_positions(self):
        for i in range(len(self.M)):
            yield (self.P[i], self.M[i])

    def filter_poslines(self, poslines, inplace=True):

        return self.filter_positions( self.get_position_indexes(poslines), inplace )

    def get_position_indexes(self, poslines):

        # create dictionary of [chr][pos] = index
        d = {}
        for i, line in enumerate(self.P):
            try:
                d[line[0]][line[1]] = i
            except KeyError:
                d[line[0]] = { line[1]: i}

        indexes = []
        counter = 0
        for line in poslines:
            if not line: continue
            try:
                counter += 1
                indexes.append( d[line[0]][int(line[1])])
            except KeyError:
                cerr('[I - warning: position not found: %s %s]' % (line[0], line[1]))

        cerr('[I - warning: only found %d out of %d positions]' % (len(indexes), counter))

        return indexes



    def filter_positions(self, posindex, inplace=True):
        posindex = np.sort(posindex)
        if inplace:
            self.df_M = self.df_M.iloc[posindex, :]
            self.M = self.df_M.values
            self.df_P = self.df_P.iloc[posindex, :]
            self.P = self.df_P.values
            return self

        new_reg = Region(self.name)
        new_reg.df_M = self.df_M.iloc[posindex, :]
        new_reg.M = self.df_M.values
        new_reg.df_P = self.df_P.iloc[posindex, :]
        new_reg.P = self.df_P.values
        return new_reg

    def filter_samples(self, indvindex):
        self.df_M = self.df_M.iloc[:, indvindex]
        self.M = self.df_M.values


    def filter_mac(self, mac = 1, inplace=True):

        # get posindex whose MAC >= mac
        snpindex = self.get_snpindex( mac = mac )
        cerr('[I - filtering MAC = %d from %d SNPs to %d SNPs]'
            % (mac, len(allele_mac), len(snpindex)))

        return self.filter_positions(snpindex, inplace)

    def get_mac(self):
        allele_0 = np.count_nonzero(self.M < 0.5, axis=1)
        allele_1 = len(self.M[0]) - allele_0
        return np.minimum(allele_0, allele_1)

    def get_snpindex(self, mac = 0):
        """ get SNP index with provided criteria, please extend as necessary """

        snpindexes = []
        if mac > 0:
            # get posindex whose MAC >= mac
            snpindexes.append( np.flatnonzero( self.get_mac() >= mac ) )

        return functools.reduce( np.intersect1d, snpindexes)

    def get_sample_indexes(self, samples):
        cols = self.df_M.columns.values
        sidx = np.argsort(cols)
        #import IPython; IPython.embed()
        return sidx[np.searchsorted(cols, samples, sorter=sidx)]


    def save(self, fmt, prefixname=None, autofilename=False, with_position=False):

        # we make assumption on the type of data
        if self.M.dtype == np.int8:
            datatype = 'nalt'
        else:
            datatype = 'ralt'

        if autofilename:
            prefixname = '%s-%d-%d' % (
                    'r' if datatype == 'ralt' else 'n',
                    len(samples), len(whole_region.M)
            )

        outmatrix = prefixname + ('.ralt' if datatype == 'ralt' else '.nalt')
        if fmt == 'pickle':
            outmatrix = outmatrix + '.pickle.gz'
            self.df_M.to_pickle(outmatrix)
        elif fmt == 'npy':
            outmatrix = outmatrix + '.npy.gz'
            with gzopen(outmatrix, 'wb') as f:
                a = np.array( [ np.array(self.df_M.columns), self.df_M.values] )
                np.save(f, a)
        else:
            outmatrix = outmatrix + '.txt.gz'
            self.df_M.to_csv(outmatrix, sep='\t', index=False)
        cerr('[I - writing genotype data to %s]' % outmatrix)

        if with_position:
            outpos = prefixname + '.pos.txt.gz'
            self.df_P.to_csv(outpos, sep='\t', index=False)

            cerr('[I - writing position data to %s]' % (outpos))



class PositionParser(object):

    def __init__(self, args):

        # positions
        if args.posfile is None:
            cexit('ERROR: required --posfile')
        self.posfilename = args.posfile
        self.posfile = None
        self.posfile_header = None
        self.positions = None
        self.header = None
        self.n = args.n
        self.df = None
        self.M = None

    def read_data(self):
        if self.M is None:
            self.df = pd.read_csv(self.posfilename, delimiter='\t',
                        nrows=self.n if self.n > 0 else None)
            self.M = self.df.values
            self.header = self.df.columns

    def get_df(self):
        self.read_data()
        return self.df

    def get_M(self):
        self.read_data()
        return self.M

    def get_positions(self, indexes):
        pass

    def get_indexes(self, positions):
        pass

    def get_posinfo(self):

        if not self.posfilename:
            c = 1
            while True:
                yield ('NA', 'NA', c)
                c += 1
            return

        if not self.posfile:
            self.posfile = gzopen(self.posfilename, 'rt')

        if not self.posfile_header:
            self.posfile.seek(0)
            self.posfile_header = next(self.posfile).strip().split('\t')

        for line in self.posfile:
            yield line.strip().split('\t')


def init_argparser(p=None, with_group=True, with_position=True):

    if p is None:
        p = arg_parser('Genotype file parser')

    if with_group:
        p = grpparser.init_argparser( p )

    if with_position:
        p.add_argument('--posfile', default=None)
        p.add_argument('--includepos', default='')

    p.add_argument('-n', type=int, default=-1)
    p.add_argument('--fmt', default='tab', choices = ['tab', 'pickle'])
    p.add_argument('infile')

    return p


class NAltLineParser(object):

    # n_alt file
    # ----> s1     s2 s3
    # SNP1    0    0    1
    # SNP2    0    2    1
    # SNP3    0    0    0
    #
    # there are 2 variant of data reader;
    # 1) using pandas.read_csv (fast) or
    # 2) using numpy.fromfile (memory )

    def __init__(self, args, datatype='nalt', with_group=True, with_position=True):

        self.group_parser = grpparser.GroupParser( args ) if with_group else None
        self.position_parser = PositionParser( args ) if with_position else None

        self.infile = args.infile
        self.fmt = args.fmt
        self.n = args.n

        self.dtype = np.int8 if datatype=='nalt' else np.float
        #self.convert_data = lambda line: np.loadtxt(io.StringIO(line),
        #                        dtype = self.dtype, delimiter='\t')

        #self.convert_data = lambda line: pd.read_csv(io.StringIO(line),
        #                        dtype = dtype, delimiter='\t', header=None).values

        self.convert_data = lambda line: np.fromfile(io.StringIO(line),
                                dtype = dtype, delimiter='\t')

        self.df = None
        self.M = None
        self.samples = None

        self.parse_samples()


    def read_data(self):
        cerr('[I - reading file %s]' % self.infile)
        if not self.M:

            if self.fmt == 'pickle':
                self.df = pd.read_pickle(self.infile)
            elif self.fmt == 'npy':
                with gzopen(self.infile, 'rb') as f:
                    a = np.load(f)
                    self.df = pd.DataFrame(columns = a[0], data = a[1])
            else:
                self.df = pd.read_csv(self.infile, dtype=self.dtype, delimiter='\t',
                            nrows=self.n if self.n > 0 else None)
            self.samples = self.df.columns
            self.M = self.df.values


    def parse_samples(self):
        if self.samples is None:
            self.read_data()
        return self.samples


    def get_sample_header(self):
        return '\t'.join(self.samples)


    def parse_grouping(self):

        # assign samples to group
        samples = self.parse_samples()
        groups = self.group_parser.assign_groups( samples )
        return groups


    def parse_whole(self, n=-1, mask=None):
        """ parse whole genome, return a Region """

        region = Region('whole')
        region.df_M = self.df
        region.M = self.M
        if self.position_parser:
            region.df_P = self.position_parser.get_df()
            region.P = self.position_parser.get_M()
        else:
            region.df_P = region.P = None

        return region


    def parse_lines(self, n=-1):
        """ parse whole file yielding line-by-line """
        pass

    def parse_chromosomes(self):
        pass


    def parse_genes(self):
        pass


    def subset(self, L, name='subset'):

        new_region = Region(name)
        for l in L:
            new_region.append( self.P[l], self.M[l] )

        return new_region


def open(filename, fmt='tab', datatype='ralt', n=-1):
    """ simple function to open a nalt or ralt file
    """

    from types import SimpleNamespace
    args = SimpleNamespace( infile=filename, fmt=fmt, datatype=datatype, n=n)
    return NAltLineParser(args, with_group=False, with_position=False)


def get_position_indexes(poslines, pos_dataframe):

        # create dictionary of [chr][pos] = index
        d = {}
        for i, line in enumerate(pos_dataframe.values):
            try:
                d[line[0]][line[1]] = i
            except KeyError:
                d[line[0]] = { line[1]: i}

        indexes = []
        counter = 0
        for line in poslines:
            if not line: continue
            try:
                counter += 1
                indexes.append( d[line[0]][int(line[1])])
            except KeyError:
                cerr('[I - warning: position not found: %s %s]' % (line[0], line[1]))

        cerr('[I - warning: only found %d out of %d positions]' % (len(indexes), counter))

        return indexes

