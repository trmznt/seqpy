#!/usr/bin/env spcli

from seqpy import cout, cerr
from seqpy.cmds import arg_parser
from seqpy.core import bioio
from collections import defaultdict


def init_argparser():

    p = arg_parser('Convert VCF file to sequence file (eg. fasta)')
    p.add_argument('-o', '--outfile', default='-')
    p.add_argument('--heterozygosity', type=float, default=0.0)
    p.add_argument('--chr', default='')
    p.add_argument('--opts', default='')
    p.add_argument('vcffile')
    return p


def main( args ):

    vcf2seq( args )


class Filters(object):

    NoIndel = 'NoIndel'
    LowQual = 'LowQual'
    MissingThreshold = 'MissingThreshold'
    HetThreshold = 'HetThreshold'
    DP = 'DP'
    AD = 'AD'
    MAC = 'MAC'
    MAF = 'MAF'

    keywords = [ NoIndel, LowQual, MissingThreshold, HetThreshold, DP, AD, MAC, MAF ]

def vcf2seq( args ):

    vcf2seqhelper = VCF2SeqHelper( args.vcffile, args.chr,
                    'NoIndel,LowQual,MissingThreshold=0.05,HetThreshold=0.25,' + args.opts )
    vcf2seqhelper.parse()
    mseq = vcf2seqhelper.get_multisequence()
    cout('Report:')
    for k,v in vcf2seqhelper.chr_used.items():
        cout(' %s\t%d' % (k,v))
    cout('Writing to %s' % args.outfile)
    bioio.save( mseq, args.outfile )


class VCFLineParser(object):

    def __init__(self, vcffile, chroms = None, filters='', **kwargs):
        self.vcffile = vcffile
        self._hdl = open( self.vcffile )

        # parse filter
        cerr('Filters: %s' % filters)
        self.filters = {}
        for filter_item in filters.split(','):
            if not filter_item: continue
            if '=' in filter_item:
                k,v = filter_item.split('=',1)
                k = k.strip()
                self._check_keyword(k)
                self.filters[k.strip()] = float(v.strip())
            else:
                filter_item = filter_item.strip()
                self._check_keyword(filter_item)
                self.filters[filter_item] = True

        self.sample_labels = None
        if ',' in chroms:
            chroms = [ c.strip() for c in chroms.split(',') ]
        self.chroms = chroms
        self.installed_filters = [
                self._filter_MissingThreshold,
                self._filter_HetThreshold,
                self._filter_MAF,
                self._filter_MAC,
        ]
        self.dp = 0 if 'DP' not in kwargs else int(kwargs['DP'])
        self.ad = 0 if 'AD' not in kwargs else int(kwargs['AD'])
        self.init_params(**kwargs)


    def init_params(self):
        raise NotImplementedError()


    def init_samples(self, line):
        """ extend this for initializing variables """
        self.sample_labels = line.split('\t')[9:]


    def add_filters(self, *args):
        for arg in args:
            self.installed_filters.append( arg )


    def parse_samples(self, snp_info, data_items, line = None):
        """ extend this for analysing sample data for every SNP line
            return True to indicate that data_items can be processed further,
            otherwise indicate unused data_items"""

        for filter_func in self.installed_filters:
            if not filter_func(snp_info, data_items):
                return False

        return True

        raise NotImplementedError()

    def finalize_samples(self):
        raise NotImplementedError()

    def process_commentline(self, line):
        pass

    def process_sampleline(self, line):
        self.init_samples( line )


    def parse(self):

        for line in self._hdl:

            line = line.strip()

            if not line:
                continue

            if line.startswith('##'):
                self.process_commentline(line)
                continue

            if line.startswith('#CHROM'):
                self.process_sampleline( line )
                continue

            tokens = line.split('\t')
            chrom, pos, posid, ref, alt, qual, filters, info, format = tokens[:9]

            if self.chroms and chrom not in self.chroms:
                continue

            # filtering based on the filter column of each samples, aka filter by
            # VCF keyword

            if self.filters:
                if filters != '.':
                    filtered = False
                    for filterkey in filters.split(','):
                        if filterkey in self.filters:
                            filtered = True
                            break

                    if filtered:
                        continue

                if 'NoIndel' in self.filters:

                    if len(ref) > 1:
                        continue

                    if len(alt) > 1:
                        alleles = alt.split(',')
                        indels_flag = False
                        for allele in alleles:
                            if len(allele) > 1:
                                indels_flag = True
                                break
                        if indels_flag:
                            continue

            ref = ref.encode('ASCII')
            if len(alt) > 1:
                alt = ''.join( alt.split(',') ).encode('ASCII')
            else:
                alt = alt.encode('ASCII')


            if posid == '.':
                posid = '%s:%s' % (chrom, pos)


            # patching data_items so the genotype (gt) is redetermined by us

            data_items = list(enumerate([ x.split(':') + [x] for x in tokens[9:] ]))
            for (idx, data_item) in data_items:
                data_item[0] = self._get_genotype(data_item)

            self.parse_samples( (chrom, pos, posid, ref, alt, qual, filters, info, format),
                                data_items, line )

        self.finalize_samples()


    def _get_genotype(self, data_item):

        if self.dp > 0 and int(data_item[2]) < dp:
            return './.'
        if self.ad > 0:
            items = [ int(x) for x in data_item[1].split(',') ]
            gt = []
            for i in range(len(items)):
                if i >= ad:
                    gt.append(i)
            if len(gt) > 1:
                return '/'.join(gt)
            else:
                return '%d/%d' % (gt[0], gt[0])
        return data_item[0]


    def _check_keyword(self, keyword):
        if keyword not in Filters.keywords:
            cerr('Filter keyword: %s is not recognized' % keyword)
            sys.exit(1)

    def _filter_MissingThreshold(self, snp_info, data_items):
        """ This filters the proportion of samples with missing SNP at particular
            SNP position """

        if 'MissingThreshold' in self.filters:

            # count missing haplotype
            missing = 0
            for (idx, data_item) in data_items:
                gt = data_item[0]
                if gt in ['./.', '.']:
                    missing += 1

            if missing/len(data_items) >= self.filters['MissingThreshold']:
                cout( 'SNP ID: %s did not pass missing threshold.' % snp_info[2] )
                return False

        return True


    def _filter_HetThreshold(self, snp_info, data_items):
        """ This filters the proportion of samples with heterozygote SNP at particular
            SNP position """

        if 'HetThreshold' in self.filters:

            # count heterozygosity
            hets = 0
            for (idx, data_item) in data_items:
                gt = data_item[0]
                if gt not in ['0/0', '1/1', '2/2', '3/3', '0', '1', '2', '3']:
                    hets += 1

            if hets/len(data_items) >= self.filters['HetThreshold']:
                cout( 'SNP ID: %s did not pass heterozygosity threshold.' % snp_info[2] )
                return False
        return True


    def _count_alleles(self, data_items):

        alleles = [0,0,0,0]
        missing = 0
        for (idx, data_item) in data_items:
            gt = data_item[0]
            if gt in ['.', './.']:
                missing += 1
                continue
            for a in [ int(x) for x in gt.split('/') ]:
                alleles[a] += 1
        return (alleles, missing)


    def _filter_MAC(self, snp_info, data_items, asutil=False):
        """ This filter the minimum MAC (Minor Allele Count) at particular SNP position """
        # set asutil=True if you want to get the Allele Count, instead of using this
        # as a filter

        if 'MAC' in self.filters or asutil:

            (alleles, missing) = self._count_alleles(data_items)
            mac = sum(alleles) - max(alleles)
            if mac < self.filters['MAC']:
                cerr('SNP IdL %s did not pass MAC threshold.' % snp_info[2])
                return False

        return True

        if False:

            # count MAC
            alleles = [0,0,0,0,0,0,0]
            for (idx, data_item) in data_items:
                gt = data_item[0]
                if gt in ['0', '0/0']:
                    alleles[0] += 1
                elif gt in ['1', '1/1']:
                    alleles[1] += 1
                elif gt == '0/1':
                    alleles[0] += 1
                    alleles[1] += 1
                elif gt in ['2', '2/2']:
                    alleles[2] += 1
                elif gt == '0/2':
                    alleles[0] += 1
                    alleles[2] += 1
                elif gt == '1/2':
                    alleles[1] += 1
                    alleles[2] += 1

            total_count = sum(alleles)
            mac = total_count - max(alleles)
            if asutil: return
            if mac < self.filters['MAC']:
                cerr( 'SNP ID: %s did not pass MAC threshold.' % snp_info[2] )
                return False


    def _filter_MAF(self, snp_info, data_items):
        """ This filters the cumulative minimum MAF threshold at particular SNP position """

        if 'MAF' in self.filters:

            # count MAF
            (alleles, missing) = self._count_alleles( data_items )
            maf = 1 - max(alleles) / (sum(alleles))
            cerr('SNP ID: %s >> ' % snp_info[2] + str(alleles) + ' ' + str(missing))
            if maf < self.filters['MAF']:
                cerr('SNP ID: %s with MAF: %1.3f did not pass threshold.' % (snp_info[2], maf))
                return False

        return True



class VCF2Filter(VCFLineParser):

    def init_params(self):
        self.mseq = bioio.multisequence()
        self.chr_used = defaultdict(int)


    def init_samples(self, line):
        super().init_samples(line)
        # create multisequence and populate with similar samples
        for label in self.sample_labels:
            self.mseq.append( bioio.biosequence(label=label) )


    def parse_samples(self, snp_info, data_items):

        if not super().parse_samples( snp_info, data_items ):
            return

        (chrom, pos, posid, ref, alt, qual, filters, info, format) = snp_info

        if 'MissingThreshold' in self.filters:

            # count missing haplotype
            missing = 0
            for (idx, data_item) in data_items:
                gt = data_item[0]
                if gt == './.':
                    missing += 1

            if missing/len(data_items) >= 0.05:
                cout( 'SNP ID: %s did not pass missing threshold.' % posid )
                return


        if 'HetThreshold' in self.filters:

            # count heterozygosity
            hets = 0
            for (idx, data_item) in data_items:
                gt = data_item[0]
                if gt not in ['0/0', '1/1', '2/2']:
                    hets += 1

            if hets/len(data_items) >= 0.33:
                cout( 'SNP ID: %s did not pass heterozygosity threshold.' % posid )
                return

        if 'MAF' in self.filters:

            # count MAF
            refs = 0
            for (idx, data_item) in data_items:
                gt = data_item[0]
                if gt == '0/0':
                    refs += 1
            maf = refs/len(data_items)
            if maf > 0.5:
                maf = 1 - maf
            if maf < self.filters['MAF']:
                print( 'SNP ID: %s did not pass MAF threshold.' % posid )
                return


        for (idx, data_item) in data_items:
            gt = data_item[0]
            if gt == '0/0':
                self.mseq[idx].append( ref[0] )
            elif gt == '1/1':
                self.mseq[idx].append( alt[0] )
            elif gt == '2/2':
                self.mseq[idx].append( alt[1] )
            else:
                self.mseq[idx].append( ord('N') )

        # reporting purposes
        self.chr_used[chrom] += 1


    def finalize_samples(self):
        pass


    def get_multisequence(self):
        return self.mseq
