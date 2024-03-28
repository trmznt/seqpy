#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

__copyright__ = '(c) 2024, Hidayat Trimarsanto <trimarsanto@gmail.com'
__license__ = 'MIT'
__version__ = '2024.03.03'


import sys
import pathlib

sys.path.append(pathlib.Path(__file__).parent.parent.as_posix())

from seqpy import subcommands


def main():

    cmds = subcommands.SubCommands(
        modules=['seqpy.cmds'],
        allow_any_script=True,
        allow_shell=True
    )

    cmds.main()


if __name__ == '__main__':
    main()

# EOF
