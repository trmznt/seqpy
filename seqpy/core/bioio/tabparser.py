# tabparser

from seqpy import cout, cerr, cexit, gzopen
from seqpy.cmds import arg_parser

from seqpy.core.bioio import grpparser

import numpy as np
import pandas
import attr

# requires scikit-allel
try:
    import allel
except:
    cexit('ERR: require properly installed scikit-allel!')


def init_argparser(p=None):

    if p is None:
        p = arg_parser('Genotype file parser')

    p = grpparser.init_argparser( p )

    p.add_argument('--posfile', default=None)
    p.add_argument('--includepos', default='')
    p.add_argument('infile')

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


        self.group_parser = grpparser.GroupParser( args )
        self.infile = gzopen(args.infile, 'rt')
        self.posfilename = args.posfile
        self.position = None
        self.posfile_header = None
        self.posfile = None
        self.sample_header = None
        self.samples = None

        # read included positions
        self.include_positions = {}
        if args.includepos:
            with open(args.includepos) as infile:
                next(infile)
                for line in infile:
                    tokens = line.strip().split('\t')
                    self.include_positions[ (tokens[0], tokens[1]) ] = True

        # need to read header of genotype
        self.parse_sample()


    def parse_sample(self):
        if not self.samples:
            self.infile.seek(0)
            self.sample_header = next(self.infile).strip()
            self.samples = self.sample_header.strip().split('\t')
        return self.samples


    def parse_position_header(self):
        if not self.posfilename:
            cexit('E: need --posfile')
        self.posfile = gzopen(self.posfilename)
        self.posfile_header = next(self.posfile).strip()


    def parse_grouping(self):

        # assign samples to group
        samples = self.parse_sample()
        groups = self.group_parser.assign_groups( samples )
        return groups

    @property
    def sample_idx(self):
        return self.group_parser.sample_idx

    def parse_grouping_XXX(self):

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

        if not self.posfile:
            self.parse_position_header()

        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile)):
            pos_line, geno_line = paired_line
            tokens = geno_line.strip().split('\t')
            pos_info = pos_line.strip('\n').split('\t')

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

            tokens = line.strip().split('\t')

            if len(tokens) != len(self.samples):
                cexit('E: genotype file does not match sample number at line %d' % idx)

            M.append( self.translate(tokens) )

        self.parse_position(maxline)
        return M


    def parse_raw_lines(self, maxline=-1):
        """ this is a generator, returning (posline, genoline) including new lines
        """

        if not self.posfile:
            self.parse_position_header()

        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile)):
            if maxline > 0 and idx >= maxline:
                break

            yield( paired_line )


    def parse_position(self, maxline=-1):

        if not self.posfile_header:
            self.parse_position_header()

        self.position = []
        for idx, line in enumerate(self.posfile):
            if maxline > 0 and idx >= maxline:
                break
            self.position.append( line.strip('\n').split('\t') )

        return self.position

    def get_sample_header(self, bytestring=False):
        if not self.samples:
            cexit('E: need to parse sample header first')
        header = '\t'.join( self.samples )
        if bytestring:
            return header.encode('UTF-8')
        return header

    def get_position_header(self):
        if not self.posfile_header:
            self.parse_position_header()
        return self.posfile_header


    def parse_haplotypes(self, maxline=-1):
        """ this return a list like the following:
        [   '0000022020',
            '0002020000' ]
        """

        if not self.posfile:
            self.parse_position_header()

        M = []
        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile) ):
            if maxline > 0 and idx >= maxline:
                break

            posline, genoline = paired_line
            if self.include_positions:
                posinfo = posline.strip('\n').split('\t')
                if (posinfo[0], posinfo[1]) not in self.include_positions:
                    continue
            tokens = genoline.strip().split('\t')

            M.append( x[0] for x in tokens )

        cerr('I: haplotyping for %d SNP positions' % len(M))

        # do transpose
        M_t = [ *zip( *M) ]
        H = [ ''.join( x ).encode('UTF-8') for x in M_t ]
        return H


    def parse_np_haplotypes(self, maxline=-1):
        """ this return a numpy array haplotypes
        [   [0, 0, 0, 0, 2, 2, 0, 2, 0],
            [0, 0, 0, 2, 0, 0, -1, -0. -1] ]
        """

        token2value = { '0': 0, '1': 1, '2': 2, '-': -1}

        S = len(self.samples)

        if self.include_positions:
            M = np.zeros( (S, len(self.include_positions)), np.int8 )

            l = 0
            for (idx, paired_line) in enumerate( zip(self.posfile, self.infile) ):
                posline, genoline = paired_line
                posinfo = posline.strip('\n').split('\t')
                if (posinfo[0], posinfo[1]) in self.include_positions:
                    tokens = genoline.strip().split('\t')
                    if len(tokens) != S:
                        cexit('E: inconsistent number of samples!')

                    for i in range(S):
                        M[i, l] = token2value[ tokens[i][0] ]

                    l += 1

        else:
            # we need to parse positions first
            positions = self.parse_position()
            L = maxline if maxline > 0 else len(positions)
            M = np.zeros( (len(S), L), np.int8 )

            for (idx, genoline) in enumerate(self.infile):
                if idx >= L:
                    break

                tokens = genoline.strip().split('\t')
                if len(tokens) != S:
                    raise RuntimeError('E: inconsistent number of samples!')

                for i in range(S):
                    M[i, idx] = token2value[ tokens[i][0] ]

        return M


    def parse_genes(self):
        """ a generator that returns:
            ( (chrom, pos, gene), positions, genotype matrix )
        """

        if not self.posfile:
            self.parse_position_header()

        region = None

        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile) ):
            posline, genoline = paired_line
            posinfo = posline.strip('\n').split('\t')
            gene = posinfo[4]
            if not gene:
                if region:
                    yield region
                    region = None
                continue
            if region is None:
                region = Region(gene, [posinfo], [genoline.strip().split('\t')])
            elif gene == region.name:
                region.append(posinfo, genoline.strip().split('\t'))
            elif gene != region.name:
                yield region
                region = Region(gene, [posinfo], [genoline.strip().split('\t')])

        if region:
            yield region


    def parse_chromosomes(self):
        """ a generator that returns a region of whole chromosome """
        if not self.posfile:
            self.parse_position_header()

        region = None

        for (idx, paired_line) in enumerate( zip(self.posfile, self.infile) ):
            posline, genoline = paired_line
            posinfo = posline.strip('\n').split('\t')
            chrom = posinfo[0]
            if not chrom:
                if region:
                    yield region
                    region = None
                continue
            if region is None:
                region = Region(chrom, [posinfo], [genoline.strip().split('\t')])
            elif chrom == region.name:
                region.append(posinfo, genoline.strip().split('\t'))
            elif chrom != region.name:
                yield region
                region = Region(chrom, [posinfo], [genoline.strip().split('\t')])

        if region:
            yield region

    def parse_haplogenes(self):
        """ return a list of haplotypes corresponding on the gene,
            missing SNPs will be counted as another genotype, so make
            sure we do not have missing calls\
        """

        M = []
        for (idx, data) in enumerate(self.parse_genes()):
            geninfo, posinfo, genotypes = data


@attr.s
class Region(object):
    name = attr.ib()
    P = attr.ib()
    A = attr.ib()
    G = None
    H = None


    def append(self, posinfo, genotypes):
        self.P.append(posinfo)
        self.A.append(genotypes)


    def haplotypes(self):
        """ return haplotypes """

        if self.H:
            return self.H

        M = []
        for genotypes in self.A:
            M.append( x[0] for x in genotypes )

        #cerr('I: haplotyping for %d SNP positions' % len(M))

        # do transpose
        M_t = [ *zip( *M) ]
        self.H = [ ''.join( x ).encode('UTF-8') for x in M_t ]
        return self.H


    def encode_haplotypes(self, haplotypes=None):

        if not haplotypes:
            haplotypes = self.haplotypes()

        haplo_map = {}
        c = 0
        for k in haplotypes:
            if k in haplo_map: continue
            haplo_map[k]  = c
            c += 1

        en_haplo = np.zeros( len(haplotypes), dtype='i2' )
        for i in range( len(haplotypes) ):
            en_haplo[i] = haplo_map[ haplotypes[i] ]

        return en_haplo


    def genotypes(self):
        """ return genotype array """

        if self.G:
            return self.G

        M = []
        for ha in self.A:
            M.append( [GenotypeLineParser.haploid_na_translator[k] for k in ha] )

        self.G = M
        return self.G

    def positions(self):
        return self.P

