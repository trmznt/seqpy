# genoparser.py
#

from seqpy import cout, cerr
from seqpy.cmds import arg_parser

from collections import defaultdict

import numpy as np

# requires scikit-allel
try:
    import allel
except:
    cerr('ERR: require properly installed scikit-allel!')


def init_argparser(p=None):

    if p is None:
        p = arg_parser("Genotype file parser")
    p.add_argument('--samplefile', required=True)
    p.add_argument('--legendfile', required=True)
    p.add_argument('--gff3file', required=True)
    p.add_argument('infile')

    return p


class GenotypeParser(object):

    # genotype file
    # ----> sample1 sample2 sample3 sample4
    # SNP1
    # SNP2
    # SNP3
    # etc

    diploid_translator = { '0': [0, 0], '1': [1, 1], 'NA': [-1, -1], '0.5': [0, 1] }
    haploid_translator = { '0': [0], '1': [2], 'NA': [-1], '0.5': [1] }
    diploid_na_translator = { '0': [0, 0], '1': [1, 1], 'NA': [-1, -1], '0.5': [0, 1] }
    haploid_na_translator = { '0': [0], '1': [2], 'NA': [1], '0.5': [1] }


    def __init__(self, genofile, samplefile, legendfile, gff3file):
        self.genofile = genofile
        self.samplefile = samplefile
        self.legendfile = legendfile
        self.gff3file = gff3file
        self.groups = None
        self.gff3 = None
        self.sample_idx = []
        self.translator = self.diploid_na_translator


    def parse_samples(self):
        # sample file is comma or tab-delimited file
        if self.groups is not None:
            return self.groups
        groups = {}
        sample_idx = []
        with open(self.samplefile) as f:
            line = next(f)
            delim = ',' if line.count(',') > line.count('\t') else '\t'
            for (idx, line) in enumerate(f):
                tokens = line.strip().split(delim)
                group = tokens[1]
                if group in groups:
                    groups[group].append( idx )
                else:
                    groups[group] = [ idx ]
                sample_idx.append( idx )

            self.groups = groups
            self.sample_idx = set(sample_idx)
            self.no_of_samples = len(sample_idx)
        return self.groups

    def group_info(self):
        if self.groups is None:
            self.parse_samples()
        return 'Grouping:\n' + '\n'.join(
                '  %12s %3d' % (k, len(self.groups[k])) for k in self.groups
            ) + '\n'


    def parse_gff3(self):
        self.gff3 = allel.gff3_to_recarray(self.gff3file, [ 'ID', 'Name', 'alias'])
        return self.gff3


    def parse_genes(self):
        # this is a generator that yields (gene, genotype, seqid, positions) until StopIteration

        # prepare gene list
        gff3 = self.parse_gff3()
        genes = [ r for r in gff3 if r.type == 'gene' ]

        gene_iter = iter(genes)

        with open(self.genofile) as genofile, open(self.legendfile) as legendfile:

            # sanity check for sample number in genotype file
            header = next(genofile).strip()
            # check for delim
            delim = ',' if header.count(',') > header.count('\t') else '\t'
            samples = header.split(delim)
            if len(samples) != self.no_of_samples:
                cerr('ERR: no of sample in infile (%d) does not match the sample file (%d)'
                    % (len(samples), self.no_of_samples) )
                cerr('ERR: last sample in infile is %s' % samples[-1])
                print(samples)
                import IPython
                IPython.embed()

            next(legendfile)
            gene = next(gene_iter)
            genotype = []
            positions = []

            # generate genotype_line
            for (idx, line) in enumerate(zip(genofile, legendfile)):
                geno_line, legend_line = line
                tokens = geno_line.strip().split(delim)
                legends = legend_line.replace('"', '').strip().split()
                seqid = legends[0]
                pos = int(legends[1])

                if len(tokens) != self.no_of_samples:
                    cerr('ERR: line %d - no of sample in infile does not match the sample file' % (idx + 1))

                # check if we are within current gene

                if ((seqid == gene.seqid and pos > gene.end) or
                    seqid != gene.seqid):

                    # either we are already past a gene or change a chromosome
                    # yield (gene, genotype)
                    cerr('Changing from gene %s' % gene.ID)

                    if genotype:
                        assert len(genotype) == len(positions)
                        yield( (gene, genotype, seqid, positions) )

                    genotype = []
                    positions = []
                    gene = next(gene_iter)

                if legends[0] == gene.seqid and pos < gene.start:
                    continue

                # assign our data to diploid to facilitate ploidy

                genotype.append( self.translate(tokens) )
                positions.append( pos )


    def parse_snps(self):

        """ this is generator that yield (legends, genotype) """

        with open(self.genofile) as genofile, open(self.legendfile) as legendfile:

            # sanity check for sample number in genotype file
            header = next(genofile).strip()
            samples = header.split()
            if len(samples) != self.no_of_samples:
                cerr('ERR: no of sample in infile does not match the sample file')

            next(legendfile)
            gene = next(gene_iter)
            genotype = []

            # generate genotype_line
            for (idx, line) in enumerate(zip(genofile, legendfile)):
                geno_line, legend_line = line
                tokens = geno_line.strip().split()
                legends = legend_line.replace('"', '').strip().split()
                seqid = legends[0]
                pos = int(legends[1])

                if len(tokens) != self.no_of_samples:
                    cerr('ERR: line %d - no of sample in infile does not match the sample file' % (idx + 1))

                g = self.translate(tokens)

                yield( (legends, g))


    def set_translator(self, translator):
        self.translator = translator


    def translate(self, tokens):
        return [ self.translator[k] for k in tokens ]