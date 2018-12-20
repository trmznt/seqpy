# naltparser

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

from seqpy.core.bioio import grpparser

import numpy as np
import pandas as pd
import attr
import io
import itertools

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

    def haplotypes(self):
        if self.H:
            return self.H

        self.H = np.transpose( np.array(self.M) )
        return self.H

    def ralt_to_nalt(self, majority=False, hetratio=0.25):

        from genoutils import ralt_to_nalt

        if majority:
            hetratio = -1
        for i in range(len(self.M)):
            self.M[i] = ralt_to_nalt(self.M[i], hetratio)

    def parse_positions(self):
        for i in range(len(self.M)):
            yield (self.P[i], self.M[i])


class PositionParser(object):

    def __init__(self, args):

        # positions
        self.posfilename = args.posfile
        self.posfile = None
        self.posfile_header = None
        self.positions = None
        self.header = None
        self.M = None

    def read_data(self):
        if self.M is None:
            df = pd.read_table(self.posfilename, delimiter='\t')
            self.M = df.values
            self.header = df.columns

    def get_M(self):
        self.read_data()
        return self.M

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


def init_argparser(p=None):

    if p is None:
        p = arg_parser('Genotype file parser')

    p = grpparser.init_argparser( p )

    p.add_argument('--posfile', default=None)
    p.add_argument('--includepos', default='')
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
    # 1) using pandas.read_table (fast) or
    # 2) using numpy.fromfile (memory )

    def __init__(self, args, datatype='nalt'):

        self.group_parser = grpparser.GroupParser( args )
        self.position_parser = PositionParser( args )

        self.infile = args.infile

        self.dtype = np.int if datatype=='nalt' else np.float
        #self.convert_data = lambda line: np.loadtxt(io.StringIO(line),
        #                        dtype = self.dtype, delimiter='\t')

        #self.convert_data = lambda line: pd.read_table(io.StringIO(line),
        #                        dtype = dtype, delimiter='\t', header=None).values

        self.convert_data = lambda line: np.fromfile(io.StringIO(line),
                                dtype = dtype, delimiter='\t')

        self.M = None
        self.samples = None

        self.parse_samples()


    def read_data(self):
        if not self.M:
            df = pd.read_table(self.infile, dtype=self.dtype, delimiter='\t')
            self.samples = df.columns
            self.M = df.values

    def parse_samples(self):
        if self.samples is None:
            self.read_data()
        return self.samples


    def parse_grouping(self):

        # assign samples to group
        samples = self.parse_samples()
        groups = self.group_parser.assign_groups( samples )
        return groups


    def parse_whole(self, n=-1, mask=None):
        """ parse whole genome, return a Region """

        region = Region('whole')
        region.M = self.M
        region.P = self.position_parser.get_M()

        return region


    def parse_chromosomes(self):
        pass


    def parse_genes(self):
        pass


    def subset(self, L, name='subset'):

        new_region = Region(name)
        for l in L:
            new_region.append( self.P[l], self.M[l] )

        return new_region

