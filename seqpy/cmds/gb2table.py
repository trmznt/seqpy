# read genbank records and provide label-translated fasta and tab-delimited file
# containing index, acc no, date of collection, country and label

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser('generate feature table from Genbank records')
    p.add_argument('-o', '--outfile')
    p.add_argument('--tabfile')
    p.add_argument("files", nargs='+')

    return p

def main( args ):

    tables = []
    container = bioio.multisequence()

    n = 1
    for infile in args.files:

        mseqs = bioio.load( infile )

        mseqs.sort( lambda x: x.label)

        for s in mseqs:
            tables.append( (n,
                            s.label,
                            s.attr.get('collection_date',''),
                            s.attr.get('country',''),
                            s.attr.get('isolate',''),
                            s.definition,
                ) )
            container.append( bioio.biosequence('%04d' % n, s.seq.upper()) )
            n += 1

    # write to output file
    tabfile = open(args.tabfile, 'w')
    tabfile.write('LABEL\tACCNO\tDATE\tCOUNTRY\tISOLATE\tDEFINITION\n')
    tables.sort()
    for r in tables:
        tabfile.write('%04d\t%s\t%s\t%s\t%s\t%s\n' % r)
    tabfile.close()

    bioio.save( container, args.outfile)

