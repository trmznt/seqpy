

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("read from various sequence format")
    p.add_argument('-o', '--outfile')
    p.add_argument('--src-isolate', action='store_true')
    p.add_argument('--degap', action='store_true', default=False)
    p.add_argument('--src', action='append')
    p.add_argument('--definition')
    p.add_argument('--summary', action='store_true')
    p.add_argument('--translatename', action='append_const', dest='io_opts',
            const='translatename', help = 'translate name in NEXUS format')
    p.add_argument('--accno', action='store_true', default=False)
    p.add_argument('--minlen', type=int, default=-1)
    p.add_argument('--maxlen', type=int, default=-1)
    p.add_argument('--maxN', type=float, default=0.0)
    p.add_argument('--sort', default=None, choices=['length', 'label'])
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

    append_attributes( container, args.src, args.src_isolate, args.definition )

    if args.accno:
        set_label_to_accno( container )

    if args.degap:
        container.degap()

    if args.minlen > 0 or args.maxlen > 0 or args.maxN > 0:
        new_container = bioio.multisequence()
        for s in container:
            if args.minlen > 0 and len(s) < args.minlen:
                continue
            if args.maxlen > 0 and len(s) > args.maxlen:
                continue
            if args.maxN > 0 and s.seq.count(b'N')/len(s) > args.maxN:
                continue
            new_container.append( s )

        container = new_container

    if args.sort:
        if args.sort.startswith('len'):
            container.sort(lambda x: len(x), reverse=True)

    if args.summary:
        for s in container:
            seq = s.seq.upper()
            print( ">%s\nA:%d\tC:%d\tG:%d\tT:%d\t-:%d" % (s.label.decode('ASCII'),
                seq.count(b'A'), seq.count(b'C'), seq.count(b'G'), seq.count(b'T'),
                seq.count(b'-')) )

    if args.outfile:
        bioio.save( container, args.outfile, options = args.io_opts or [] )


def append_attributes( multiseq, attributes, isolatelabel=False, definition=None ):

    for s in multiseq:
        if s.attr is None:
            s.attr = {}
        if attributes:
            for attr in attributes:
                attr = attr.encode('ASCII')
                k,v = attr.split(b'=',1)
                s.attr[k] = v
        if isolatelabel:
            s.attr[b'isolate'] = s.label
        if definition:
            s.definition = definition.encode('ASCII')
            if b'%isolate' in s.definition and b'isolate' in s.attr:
                s.definition = s.definition.replace(b'%isolate', s.attr[b'isolate'])

def set_label_to_accno( multiseq ):

    for s in multiseq:
        if s.label.startswith('gi|'):
            label_components = s.label.split('|')
            if label_components[2] in ['gb', 'dbj', 'ref', 'emb']:
                s.label = label_components[3]
