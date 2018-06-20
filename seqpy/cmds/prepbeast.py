# read a fasta file and tab-delimited genbank record file (from gb2table)
# and provide a new fasta file that is proper for BEAST

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

import re

re_date = re.compile('(\d{4})')

def init_argparser():

    p = arg_parser('generate new fasta file suitable for BEAST')
    p.add_argument('-o', '--outfile')
    p.add_argument('--tabfile')
    p.add_argument('infile')

    return p


def main( args ):

    # read tables

    tables = {}
    tabfile = open(args.tabfile)
    next(tabfile)
    for line in tabfile:
        items = line.strip().split('\t')
        tables[items[0]] = items


    mseq = bioio.load( args.infile )
    for s in mseq:
        rec = tables[s.label]
        mo = re_date.search( rec[2] )
        if mo:
            year = mo.group()
        else:
            year = '-'
        #print('%s/%s/%s' % (s.label, rec[3], year))
        s.label = '%s/%s' % (s.label, year)

    bioio.save( mseq, args.outfile )
