__copyright__ = '''
seqpy - sequence processing library for python

(c) 2006 - 2013 Hidayat Trimarsanto <anto@eijkman.go.id> / <trimarsanto@gmail.com>

All right reserved.
This software is licensed under GPL v3 or later version.
Please read the README.txt of this software.
'''

## IMPLEMENTATION
##
## file/stream is opened in binary mode (sequence is considered as binary data)
## sequence label is Unicode string to accomodate special characters
##


from .multisequence import multisequence
from .biosequence import biosequence
import re


def d(bytedata):
    """ decode to unicode UTF-8 string """
    return bytedata.decode('UTF-8')

def e(text):
    """ encode from unicode UTF-8 to bytes """
    return text.encode('UTF-8')


def parse_definition_line(line):
    mo = re.match(r'^([^\[\]]+)(\[.*\])*([^\[\]]+)*$', line)
    seqlabel = mo.group(1)
    attrline = mo.group(2)
    definition = mo.group(3)
    if attrline:
        attr_groups = re.findall(r'\[(.*?)\]', attrline)
        attrs = []
        for attr in attr_groups:
            if ';' in attr:
                splits = attr.split(';')
                attrs.extend([ s.strip() for s in splits ])
            else:
                attrs.append( attr.strip() )

    else:
        attrs = None
    return seqlabel.strip(), attrs, definition


def read_fasta(stream_file, multiseq = None, options = {} ):
    """ reading from fasta-formatted stream

            >seq_name1 [organism = abc; group=samples] [ flags=go ]
            ATCGAGTCTCGAGGCTCGATCGAGGAGATCGATAGGCCTCGGAT
            ATCGCGGATTCGAGGATCGAGTCAGTCAGGAGAGATCGAGAGTC

        definition line is UTF-8 unicode string
        attributes are inside square brackets []
        sequence line is bytearray
    """
    istream = generic_open(stream_file, 'rb')

    if multiseq is None:
        multiseq = multisequence()

    current_seq = None
    master_seq = None

    for line in istream:
        line = line.strip()
        if not line:
            continue
        if line.startswith(b'>'):
            # definition line
            def_line = d(line[1:])
            if def_line.startswith('@@'):
                # controls
                label = def_line[2:].strip()
                current_seq = biosequence( label )
                multiseq.add_control(label, current_seq)

            elif def_line.startswith('@'):
                # external data
                label = def_line[1:].strip()
                if not master_seq:
                    master_seq = current_seq
                current_seq = biosequence( label )
                master_seq._extdata[ label ] = current_seq

            else:
                seqlabel, seqattrs, definition = parse_definition_line( def_line )

                current_seq = biosequence(seqlabel)
                multiseq.append( current_seq )
                if seqattrs:
                    current_seq.attr = {}
                    for attr in seqattrs:
                        if '=' in attr:
                            key, value = attr.split('=',1)
                            current_seq.attr[key.strip()] = value.strip()
                if definition:
                    current_seq.definition = definition

        else:
            current_seq.extend( line.replace(b' ',b'') )

    if hasattr(istream, 'name'):
        multiseq.filename = istream.name
    else:
        multiseq.filename = ''

    return multiseq


def build_attrline(attributes, definition):

    lines = []
    for k,v in attributes.items():
        if v is True:
            lines.append(' [%s]' % k)
        else:
            lines.append(' [%s=%s]' % (k,v))
    if definition:
        lines.append(' ')
        lines.append(definition)
    return ''.join(lines)


def write_fasta( stream_file, mseqs, options={} ):

    ostream = generic_open( stream_file, 'wb' )

    for s in mseqs:

        if 'noattr' not in options and s.attr:
            defline = '>' + s.label + build_attrline(s.attr, s.definition) + '\n'
        else:
            defline = '>' + s.label + '\n'
        ostream.write( defline.encode('UTF-8') )
        ostream.write(s.seq)
        ostream.write(b'\n')

        if s._extdata:
            for (k,v) in s._extdata.items():
                ostream.write('>@ %s\n' % k)
                ostream.write(v.seq)
                ostream.write('\n')

    if mseqs._controls:
        for (k,v) in mseqs._controls.items():
            ostream.write(b'>@@ ')
            ostream.write(k.encode('ASCII'))
            ostream.write(b'\n')
            ostream.write(v.seq)
            ostream.write(b'\n')


def read_phylip():
    pass

def write_phylip( stream_file, mseqs, options={} ):

    ostream = generic_open( stream_file, 'wb' )

    show_len = 60
    seq_no = len(mseqs)
    seq_length = len(mseqs[0])
    ostream.write( e("%s %s\n" % (seq_no, seq_length)) )

    next_pos = 0
    for seq in mseqs:
        ostream.write( e("%-13s %s\n" % (seq.label, d(seq.seq[:show_len]))) )
    next_pos += show_len
    last_pos = next_pos + show_len
    ostream.write(b"\n")
    while next_pos < seq_length:
        for seq in mseqs:
            ostream.write(b"              ")
            ostream.write(seq.seq[next_pos:last_pos])
            ostream.write(b"\n")
        next_pos = last_pos
        last_pos += show_len
        ostream.write(b"\n")

    return ostream



def read_nexus():
    pass


def write_nexus( stream_file, mseq, options = {} ):

    ostream = generic_open( stream_file, 'wb' )

    name_list = None
    if 'translatename' in options:
        name_list = []
        i = 1
        for s in mseq:
            label = str(i)
            name_list.append( (label, s.label) )
            s.set_label( label )
            i += 1

    # write NEXUS header
    ostream.write( e('#NEXUS\n') )

    # write NEXUS data block
    max_label_len = max( [ len(s.label) for s in mseq ] )
    fmt = "   %%-%ds " % max_label_len
    ostream.write( b"begin data;\n" )
    ostream.write( e("   dimensions ntax=%d nchar=%d;\n" % (len(mseq), max( [ len(s) for s in mseq ] ))) )
    ostream.write( b"   format datatype=dna interleave=no gap=- missing=?;\n" )
    ostream.write( b"   matrix\n" )
    for seq in mseq:
        ostream.write( e(fmt % seq.label) )
        ostream.write( seq.seq )
        ostream.write( b"\n" )
    ostream.write( b"   ;\n" )
    ostream.write( b"end;\n" )

    if name_list:
        ostream.write( b"begin seqpy;\n" )
        ostream.write( b"   translate\n" )
        for idx, label in name_list:
            ostream.write( e("\t\t%s\t%s\n" % (idx, label)) )
        ostream.write( b"\t\t;\n" )
        ostream.write( b"end;\n" )



class GenbankRecord(object):

    def __init__(self):
        self.locus = None
        self.accno = None
        self.version = None
        self.dblink = None
        self.keywords = None
        self.source = None
        self.organism = None
        self.features = {}
        self.ranges = []
        self.origin = bytearray()


class GenbankFeature(object):

    def __init__(self, tag, ranges):
        self.tag = tag
        self.ranges = None
        self.fields = {}


class FeatureParser(object):

    def __init__(self, gbparser):
        self.gbparser = gbparser

    def parse(self, lines, bioseq):

        for line_sections in self.line_parser(lines):

            splits = line_sections[0].split()
            tag = splits[0]
            ranges = b''.join(splits[1:])

            feature = GenbankFeature(tag, ranges)

            for line in line_sections[1:]:
                if b'=' in line:
                    key, val = d(line).lstrip()[1:].split('=', 1)
                else:
                    key, val = d(line).strip()[1:], True

                try:
                    feature.fields[key].append( val )
                except KeyError:
                    feature.fields[key] = [ val ]

            if tag.startswith(b'source'):
                bioseq.features['source'] = feature

            else:
                try:
                    bioseq.features[tag].append(feature)
                except KeyError:
                    bioseq.features[tag] = [ feature ]


    def line_parser(self, lines):

        line_buffer = None
        for line in lines:

            line = line[5:]

            if line[0] == 32 and line_buffer is not None:
                if line[16] == 47:
                    line_buffer.append( line )
                else:
                    line_buffer[-1] += b' ' + line

            else:
                if line_buffer is not None:
                    yield line_buffer
                line_buffer = [ line ]

        yield line_buffer



class GenbankParser(object):

    def __init__(self, istream):
        self.istream = istream
        self.featparser = FeatureParser(self)

        self.tag_readers = {
            b'ACCESSION': self.tag_ACCESSION,
            b'ORIGIN': self.tag_ORIGIN,
            b'DEFINITION': self.tag_DEFINITION,
            b'FEATURES': self.tag_FEATURES
        }


    def tag_PASS(self, lines, bioseq):
        pass


    def tag_LOCUS(self, lines, bioseq):
        line = b' '.join( l.strip() for l in lines )
        splits = line.split()
        bioseq.add_attr('locus', d( splits[1]) )

    def tag_DEFINITION(self, lines, bioseq):
        line = b' '.join( l.strip() for l in lines )
        bioseq.definition = d(line.split(b' ',1)[1].strip())


    def tag_ACCESSION(self, lines, bioseq):
        line = b' '.join(lines)
        accession = d(line.split()[1].strip())
        bioseq.set_label(accession)
        bioseq.add_attr('accession', accession)


    def tag_ORIGIN(self, lines, bioseq):

        for line in lines[1:]:
            seqs = line.strip().split()[1:]
            for seq in seqs:
                bioseq.extend(seq)


    def tag_FEATURES(self, lines, bioseq):

        self.featparser.parse(lines[1:], bioseq.rec)

        # moved some information from source feature to bioseq
        if 'source' in bioseq.rec.features:
            source = bioseq.rec.features['source']

            fields = source.fields

            if 'isolate' in fields:
                bioseq.add_attr('isolate', fields['isolate'][0])
            if 'strain' in fields:
                bioseq.add_attr('strain', fields['strain'][0])
            if 'country' in fields:
                bioseq.add_attr('country', fields['country'][0])
            if 'collection_date' in fields:
                bioseq.add_attr('collection_date', fields['collection_date'][0])
            if 'organism' in fields:
                bioseq.add_attr('organism', fields['organism'][0])
            if 'organelle' in fields:
                bioseq.add_attr('organelle', fields['organelle'][0])


    def line_parser(self):

        line_buffer = None
        for line in self.istream:

            line = line.rstrip()
            if not line:
                continue

            if (line[0] == 32 or 48 <= line[0] <= 57) and line_buffer is not None:
                line_buffer.append( line )

            else:
                if line_buffer is not None:
                    yield line_buffer
                line_buffer = [ line ]

        yield line_buffer


    def parse(self):

        bioseq = None
        for line_buffer in self.line_parser():

            firstline = line_buffer[0]

            if firstline.startswith(b'//'):
                yield bioseq
                bioseq = None
                continue

            elif firstline.startswith(b'LOCUS'):
                bioseq = biosequence('',b'')
                bioseq.rec = GenbankRecord()

            tag = firstline.split()[0]

            func = self.tag_readers.get(tag, self.tag_PASS)
            func(line_buffer, bioseq)



def read_genbank(stream_file, multiseq=None, options=None):

    istream = generic_open(stream_file)
    if multiseq is None:
        multiseq = multisequence()

    parser = GenbankParser(istream)
    for seq in parser.parse():
        multiseq.append(seq)

    return multiseq


def write_genbank():
    pass


class EMBLParser(object):

    def __init__(self, istream):
        self.istream = istream
        self.tag_readers = {
            b'ID': self.tag_ID,
            b'SQ': self.tag_SQ,
            b'FT': self.tag_FT
        }


    def tag_PASS(self, lines, bioseq):
        # just pass
        pass


    def tag_ID(self, lines, bioseq):
        line = b' '.join(lines)
        _, text = line.split(b' ', 1)
        accession = d(text.split(b';')[0].strip())
        bioseq.set_label(accession)
        bioseq.add_attr('accession', accession)


    def tag_FT(self, lines, bioseq):
        pass

    def tag_SQ(self, lines, bioseq):

        for line in lines[1:]:
            seqs = line.strip().split()[:-1]
            for seq in seqs:
                bioseq.extend(seq)

    def line_parser(self):

        line_buffer = None
        for line in self.istream:

            line = line.rstrip()
            if not line:
                continue

            if (line[0] == 32 or 48 <= line[0] <= 57) and line_buffer is not None:
                line_buffer.append( line )

            else:
                if line_buffer is not None:
                    yield line_buffer
                line_buffer = [ line ]

        yield line_buffer


    def parse(self):

        bioseq = None

        for line_buffer in self.line_parser():

            firstline = line_buffer[0]

            if firstline.startswith(b'//'):
                yield bioseq
                bioseq = None
                continue

            elif firstline.startswith(b'ID'):
                bioseq = biosequence('',b'')

            tag = firstline.split()[0]

            func = self.tag_readers.get(tag, self.tag_PASS)
            func(line_buffer, bioseq)


def read_embl(stream_file, multiseq=None, options=None):

    istream = generic_open(stream_file)
    if multiseq is None:
        multiseq = multisequence()

    parser = EMBLParser(istream)
    for seq in parser.parse():
        multiseq.append(seq)

    return multiseq


def write_embl():
    pass


def read_flatfile(stream_file, multiseq=None, options=None, delim='\t'):

    istream = generic_open(stream_file)
    multiseq = multisequence() if multiseq == None else multiseq

    for line in istream:
        label, seq = line.split(delim, 1)
        multiseq.append( bioseq(label, seq))

    return multiseq


def read_flatfile_csv(stream_file, multiseq=None, options=None, delim=','):

    return read_flatfile(stream_file, multiseq, options, delim)


class UnknownFormatErr(RuntimeError):
    pass

def generic_open(stream_file, mode='rb'):
    if type(stream_file) == str:
        if stream_file.lower().endswith('.gz'):
            import gzip
            return gzip.open( stream_file, mode )
        return open(stream_file, mode)
    return stream_file

def guess_parser(filename):
    """ guess filename format based on file extension """

    (base, ext) = filename.rsplit('.', 1)
    if ext.lower() == 'gz':
        (_, ext) = base.rsplit('.', 1)
    elif ext.lower() == 'seqdb':
        import seqdb
        return seqdb.open, seqdb.open
    ext = ext.upper()
    if ext in [ 'FA', 'FAS', 'FASTA', 'FST', 'EFAS', 'EFST', 'EFASTA' ]:
        return read_fasta, write_fasta
    elif ext in [ 'PHY', 'PHYLIP' ]:
        return read_phylip, write_phylip
    elif ext in [ 'GENBANK', 'GB' ]:
        return read_genbank, None
    elif ext in [ 'NEXUS', 'NXS', 'NEX' ]:
        return read_nexus, write_nexus
    elif ext in [ 'EMBL', 'EMB' ]:
        return read_embl, write_embl
    elif ext in [ 'ALN']:
        return read_clustalw, None
    elif ext in [ 'TAB', 'TSV' ]:
        return read_flatfile, None
    elif ext in [ 'CSV' ]:
        return read_flatfile_csv, None
    elif ext in [ 'ABI', 'AB1', 'SCF']:
        from seqpy.core import traceio
        return traceio.read_trace, None
    else:
        raise UnknownFormatErr( filename )

