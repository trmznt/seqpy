__copyright__ = '''
seqpy/__init__.py - part of seqpy

(c) 2006 - 2012 Hidayat Trimarsanto <anto@eijkman.go.id> / <trimarsanto@gmail.com>

All right reserved.
This software is licensed under GPL v3 or later version.
Please read the README.txt of this software.
'''

import sys, os

# prepare seqpy/extlib path to sys.path
sys.path.append( os.path.join(os.path.split( os.path.split(__file__)[0] )[0], 'seqpy/extlib' ))

def _COUT( text ):
    print( text, file=sys.stdout )

def _CERR( text ):
    print( text, file=sys.stderr, flush = True )

_cout = _COUT
_cerr = _CERR

def cout( text ):
    _cout( text )

def cerr( text ):
    _cerr( text )

def cexit( text, err=1):
    _cerr( text )
    sys.exit(err)

def set_cout(func):
    global _cout
    _cout = func

def set_cerr(func):
    global _cerr
    _cerr = func

def gzopen(filename, options='rt'):

    if filename.endswith('.gz'):
        import gzip
        return gzip.open( filename, options )
    else:
        return open( filename, options )
