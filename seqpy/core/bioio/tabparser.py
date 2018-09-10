# tabparser

from seqpy import cout, cerr
from seqpy.cmds import arg_parser

import numpy as numpy
import pandas

# requires scikit-allel
try:
    import allel
except:
    cexit('ERR: require properly installed scikit-allel!')


def init_argparser(p=None):

    if p is None:
        p = arg_parser('Genotype file parser')

    p.add_argument('--groupfile', default='')
    p.add_argument('--metafile', default='')
    p.add_argument('--column', default='1,2')
    p.add_argument('--posfile', required=True)
    p.add_argument('--includepos', default='')
    p.add_argument('--infile')

    return p

class GenotypeLineParser(object):

    # genotype file
    # ----> sample1 sample2 sample3 sample4
    # SNP1
    # SNP2
    # SNP3
    # etc

    diploid_translator = { '0': [0, 0], '2': [1, 1], '-1': [-1, -1], '1': [0, 1] }
    haploid_translator = { '0': [0], '2': [2], '-1': [-1], '1': [1] }
    diploid_na_translator = { '0': [0, 0], '2': [1, 1], '-1': [0, 1], '1': [0, 1] }
    haploid_na_translator = { '0': [0], '2': [2], '-1': [1], '1': [1] }


    def __init__(self, args):

        if args.infile.endswith('.gz'):
            import gzip
            self.infile = gzip.open( args.infile, 'rt' )
        else:
            self.infile = open(args.infile)
        self.posfile = open(args.posfile)

        if not (args.groupfile or args.metafile):
            cexit('E: required either --groupfile or --metafile')

        self.groupfile = open(args.groupfile) if args.groupfile else None
        self.metafile = open(args.metafile) if args.metafile else None
        self.delimiter = ( ',' if args.metafile[-3:].lower() == 'csv' else '\t' ) if args.metafile else None
        self.column = args.column
        self.groups = {}
        self.position = None
        self.posfile_header = None

        self.include_positions = {}
        if args.includepos:
            with open(args.includepos) as infile:
                next(infile)
                for line in infile:
                    tokens = line.split()
                    self.include_positions[ (tokens[0], tokens[1]) ] = True


    def parse_grouping(self):

        # create a dictionary of groups <> sample_idx
        if self.groupfile:
            # this is a YAML/JSON file, open with YAML
            import yaml
            grouping = yaml.load( self.groupfile )
            groups = {}
            for g in grouping:
                for s in grouping[g]:
                    groups[s] = g
            self.group_info = groups

        elif self.metafile:
            # this is a tab/comma delimited file
            metadf = pandas.read_table(self.metafile, sep=self.delimiter)
            sample_column, group_column = self.column.split(',')
            if sample_column.isdigit():
                sample_column = metadf.columns[ int(sample_column)-1 ]
            if group_column.isdigit():
                group_column = metadf.columns[ int(group_column)-1 ]

            sampledf = metadf.loc[:, [sample_column, group_column] ]
            groups = {}
            for i in range(len(sampledf)):
                r = sampledf.loc[i]
                groups[r[0]] = r[1]

            self.group_info = groups

        else:
            cexit('E: need groupfile or metafile')

        # read the header of genotype and posfile
        header = next(self.infile)
        self.posfile_header = next(self.posfile)
        samples = header.strip().split('\t')
        cerr('I: reading %s samples from genotype file')
        groups = {}
        sample_idx = []
        for idx, code in enumerate(samples):
            grp_key = self.group_info[code]
            if grp_key in groups:
                groups[grp_key].append(idx)
            else:
                groups[grp_key] = [ idx ]
            sample_idx.append( idx )

        self.samples = samples
        self.sample_idx = set(sample_idx)
        self.groups = groups
        return self.groups


    def parse(self):
        """ this is a generator, returning (pos_info, geno_array) """


        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile)):
            pos_line, geno_line = paired_line
            tokens = geno_line.split()
            pos_info = pos_line.split()

            if len(tokens) != len(self.samples):
                cexit('E: genotype file does not match sample number at line %d' % idx)

            g = self.translate(tokens)

            yield( (pos_info, g) )


    def set_translator(self, translator):
        self.translator = translator


    def translate(self, tokens):
        return [ self.translator[k] for k in tokens ]


    def parse_all(self, maxline=-1):
        """ this return a full array from the data
            as such, ensure that the memory is big enough before calling this method
        """
        M = []
        for (idx, line) in enumerate(self.infile):
            if maxline > 0 and idx >= maxline:
                break

            tokens = line.split()

            if len(tokens) != len(self.samples):
                cexit('E: genotype file does not match sample number at line %d' % idx)

            M.append( self.translate(tokens) )

        self.parse_position(maxline)
        return M


    def parse_raw_lines(self, maxline=-1):
        """ this is a generator, returning (posline, genoline)
        """

        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile)):
            if maxline > 0 and idx >= maxline:
                break

            yield( paired_line )


    def parse_position(self, maxline):

        self.position = []
        for line in self.posfile:
            self.position.append( line.strip().split() )

        return self.position

    def get_sample_header(self, bytestring=False):
        header = '\t'.join( self.samples )
        if bytestring:
            return header.encode('UTF-8')
        return header


    def parse_haplotypes(self, maxline=-1):
        """ this return a list like the following:
        [   '0000022020',
            '0002020000' ]
        """

        M = []
        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile) ):
            if maxline > 0 and idx >= maxline:
                break

            posline, genoline = paired_line
            if self.include_positions:
                posinfo = posline.split()
                if (posinfo[0], posinfo[1]) not in self.include_positions:
                    continue
            tokens = genoline.split()

            M.append( x[0] for x in tokens )

        cerr('I: reading %d SNP positions' % len(M))

        # do transpose
        M_t = [ *zip( *M) ]
        H = [ ''.join( x ).encode('UTF-8') for x in M_t ]
        return H
