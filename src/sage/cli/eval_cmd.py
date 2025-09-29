import argparse

from sage.cli.options import CliOptions
from sage.repl.preparse import preparse
from sage.all import sage_globals
import os, sys


class EvalCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the "run" command.

        INPUT:

        - ``parsers`` -- the parsers to extend.

        OUTPUT:

        - the extended parser.
        """
        parser.add_argument(
            "-c",
            "--command",
            nargs="?",
            help="execute the given command as sage code",
        )

    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options

    def run(self) -> int:
        r"""
        Execute the given command.
        """
        # Allow importing modules from the current directory, matching python behavior
        sys.path.append(os.getcwd())

        code = preparse(self.options.command)
        eval(compile(code, "<cmdline>", "exec"), sage_globals())
        return 0
