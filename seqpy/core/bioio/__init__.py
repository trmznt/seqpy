__copyright__ = '''
bioio.py - part of seqpy

(c) 2006 - 2014 Hidayat Trimarsanto <anto@eijkman.go.id> / <trimarsanto@gmail.com>

All right reserved.
This software is licensed under GPL v3 or later version.
Please read the README.txt of this software.
'''

__version__ = '2012072601'

import sys
from .parser import guess_parser
from .seqtype import DNA, RNA, PROTEIN
from .biosequence import biosequence
from .multisequence import multisequence

def load(filename, options={}, fmt=None):
    reader, _ = guess_parser(filename)
    return reader(filename, options=options)

def save(obj, filename=None, options={}, fmt=None ):
    if not filename:
        filename = obj.filename
    else:
        obj.filename = filename

    _, writer = guess_parser(filename)
    writer(filename, obj, options=options)

