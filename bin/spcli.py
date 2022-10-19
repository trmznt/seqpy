#!/usr/bin/env python3

__copyright__ = '''
spcli - part of seqpy

(c) 2011-2022 Hidayat Trimarsanto <trimarsanto@gmail.com>

All right reserved.
This software is licensed under GPL v3 or later version.
Please read the README.txt of this software.
'''

##
# spcli.py - seqpy command line interface
##

import sys
import os
import argparse
import platform
import pathlib

# prepare path to seqpy
sys.path.append(os.path.split(os.path.split(__file__)[0])[0])

import seqpy


def greet():
    seqpy.cerr('spcli - seqpy command line interface')
    seqpy.cerr('(C) 2011-2022 Hidayat Trimarsanto <trimarsanto@gmail.com>')
    seqpy.cerr(f'Host: {platform.uname().node}')


def usage():
    seqpy.cerr('  usage:')
    seqpy.cerr('    spcli scriptfile/CMD [ARGS]')
    seqpy.cerr('  try: spcli showcmds')
    sys.exit(0)


def main():

    if '_ARGCOMPLETE' in os.environ:
        autocomplete()

    greet()
    if len(sys.argv) == 1:
        usage()

    if sys.argv[1].endswith('.py'):
        # will execute a script file
        seqpy.cerr('Attempting to run script: %s' % sys.argv[1])
        with open(sys.argv[1]) as fh:
            code = compile(fh.read(), sys.argv[1], 'exec')
            sys.argv = sys.argv[1:]
            _l = {'__name__': '__spcli_main__'}
            exec(code, None, _l)
            if 'main' in _l:
                globals().update(_l)
                main = _l['main']
                if 'init_argparser' in _l:
                    init_argparser = _l['init_argparser']
                    p = init_argparser()
                    if not isinstance(p, argparse.ArgumentParser):
                        seqpy.cerr(
                            'ERR: init_argparser() did not return ArgumentParser instance')
                        sys.exit(1)

                    try:
                        p.add_argument('--debug', default=False, action='store_true',
                                       help='open ipython3 pdb console when exception occurs')
                    except argparse.ArgumentError:
                        pass

                    argp = p.parse_args(sys.argv[1:])
                    if argp.debug:
                        from ipdb import launch_ipdb_on_exception
                        with launch_ipdb_on_exception():
                            seqpy.cerr('WARN: running in debug mode')
                            main(argp)
                    else:
                        main(argp)
                else:
                    import inspect
                    if 'args' in inspect.signature(main).parameters:
                        main(args=sys.argv)
                    else:
                        main()

    elif sys.argv[1] == '-i':
        # interactive
        import IPython
        IPython.embed()

    else:
        from seqpy import cmds
        cmds.execute(sys.argv[1:])


def autocomplete():

    from seqpy.cmds import list_commands

    # prepare line
    line = os.environ.get('COMP_LINE', '')
    tokens = line.split()
    last_token = tokens[-1]

    if len(tokens) > 1 and (last_token.startswith(',') or last_token.startswith('~')):
        # let bash complete  use directory listing
        sys.exit(1)

    # prepare the completion lists
    cmds = list_commands()

    if len(tokens) == 1:
        completions = sorted(cmds)

    else:
        completions = sorted(
            [opt for opt in cmds if opt.startswith(last_token)]
        )

    # send results through fd 8
    ifs = os.environ.get('IFS', '\013')
    out_stream = os.fdopen(8, 'w')
    out_stream.write(ifs.join(completions))
    sys.exit(0)


if __name__ == '__main__':
    main()

# EOF
