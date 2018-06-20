

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("read from various sequence format")
    p.add_argument('-o', '--outfile')
    p.add_argument('--collect', type=int, default=0)
    p.add_argument('--reffile')
    p.add_argument('--dbfile')
    p.add_argument('treefile')

    return p


def main( args ):

    import dendropy

    tree = dendropy.Tree.get(path=args.treefile, schema="newick")

    pdc = tree.phylogenetic_distance_matrix()

    cerr('Reading: %d taxa' % len(tree.taxon_namespace))

    if args.collect > 0:

        ref_seqs = bioio.load( args.reffile )

        ref_taxa = []
        for taxon in tree.taxon_namespace:
            if ref_seqs.get_by_label(taxon.label) != None:
                print('appended')
                ref_taxa.append( taxon )

        cerr('Referenced: %d taxa' % len(ref_taxa))

        collected_taxa = set()
        for t1 in ref_taxa:
            d = []
            for t2 in tree.taxon_namespace[:-1]:
                d.append( (pdc(t1,t2),t2) )
            d.sort()

            for i in range(args.collect):
                collected_taxa.add( d[i][1] )
            collected_taxa.add( t1 )

        cerr('Collected: %d taxa' % len(collected_taxa))

        db_seqs = bioio.load( args.dbfile )
        mseq = bioio.multisequence()
        for taxon in collected_taxa:
            mseq.append( db_seqs.get_by_label(taxon.label) )

        bioio.save( mseq, args.outfile )
