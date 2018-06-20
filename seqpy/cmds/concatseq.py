

from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

def init_argparser():

    p = arg_parser("concatenate from various sequence alignments")
    p.add_argument('-o', '--outfile')
    p.add_argument('--src-isolate', action='store_true')
    p.add_argument('--src', action='append')
    p.add_argument('--definition')
    p.add_argument('--summary', action='store_true')
    p.add_argument('--translatename', action='append_const', dest='io_opts',
            const='translatename', help = 'translate name in NEXUS format')
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

    if args.summary:
        for s in container:
            seq = s.seq.upper()
            print( ">%s\nA:%d\tC:%d\tG:%d\tT:%d\t-:%d" % (s.label.decode('ASCII'),
                seq.count(b'A'), seq.count(b'C'), seq.count(b'G'), seq.count(b'T'),
                seq.count(b'-')) )

    if args.outfile:
        bioio.save( container, args.outfile, options = args.io_opts or [] )


def concat_sequencs( multiseqs ):

    seqnames = {}
    new_mseq = bioio.multisequence()
    for seq in multiseqs[0]:
        s = seq.clone()
        s.set_sequence( s.get_sequence() )
        new_mseq.append( s )
        seqnames[s.get_label()] = True

    for multiseq in multiseqs[1:]:
        label = s.get_label()
        if 
        
        

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

