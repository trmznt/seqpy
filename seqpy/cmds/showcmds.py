
from seqpy import cout, cerr
from seqpy.cmds import arg_parser

import os

def init_argparser():

    p = arg_parser("show available commands")
    return p


def main( args ):

    cerr("cmds directory: %s" % os.path.dirname(__file__))
    cmds = []
    for filename in os.listdir( os.path.dirname( __file__ ) ):
        if filename.endswith('.py') and filename != '__init__.py':
            cmds.append( filename[:-3] )

    cmds.sort()
    cerr("available commands:")
    for cmd in cmds:
        cerr('  %s' % cmd)


