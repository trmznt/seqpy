__copyright__ = '''
seqpy/cmds/__init__.py - part of seqpy

(c) 2006 - 2022 Hidayat Trimarsanto <anto@eijkman.go.id> / <trimarsanto@gmail.com>

All right reserved.
This software is licensed under GPL v3 or later version.
Please read the README.txt of this software.
'''

import argparse
import importlib
import seqpy


def execute(args):

    command = args[0]
    M = importlib.import_module('seqpy.cmds.' + command)
    print(M)
    parser = M.init_argparser()
    if not parser:
        seqpy.cexit('Fatal ERR: init_argparser() does not return properly')
    parser.add_argument('--debug', default=False, action='store_true',
                        help='open ipython3 pdb console when exception occurs')

    args = parser.parse_args(args[1:])
    if args.debug:
        from ipdb import launch_ipdb_on_exception
        with launch_ipdb_on_exception():
            seqpy.cerr('WARN: running in debug mode')
            M.main(args)
    else:
        M.main(args)


def arg_parser(description=None):

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--noattr', action='append_const', dest='io_opts', const='noattr',
                        help='do not read or write sequence attribute within fasta format')

    return parser

# EOF
