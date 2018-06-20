
from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("relabel sequences with numerics")
    p.add_argument('-o', '--outfile')
    p.add_argument('-t', '--tabfile')
    p.add_argument("files", nargs='+')
    return p

def main( args ):

    container = None

    for infile in args.files:

        obj = bioio.load( infile, options = args.io_opts or [] )
        cout('reading %d sequences from %s' % (len(obj), infile))
        if container is None:
            container = obj
        else:
            container += obj

    indexes = []
    counter = 0
    for s in container:
        counter += 1
        new_label = '%04d' % counter
        indexes.append( ( new_label, s.label ))
        s.label = new_label

    if args.outfile:
        bioio.save( container, args.outfile, options = args.io_opts or [] )

    if args.tabfile:
        with open(args.tabfile, 'w') as f:
            for i in indexes:
                f.write('%s\t%s\n' % i)
