
from seqpy import cout, cerr
from seqpy.cmds import arg_parser, list_commands


def init_argparser():

    p = arg_parser("show available commands")
    return p


def main(args):

    cerr("available commands:")
    for cmd in list_commands():
        cerr('  %s' % cmd)

# EOF
