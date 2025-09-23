import argparse

from sage.cli.options import CliOptions
from sage.repl.preparse import preparse_file_named
from sage.repl.load import load_cython
from sage.misc.temporary_file import tmp_filename
from sage.all import sage_globals


class RunFileCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the "run file" command.

        INPUT:

        - ``parsers`` -- the parsers to extend.

        OUTPUT:

        - the extended parser.
        """
        parser.add_argument(
            "file",
            nargs="*",
            help="execute the given file as sage code",
        )

    def __init__(self, options: CliOptions):
        r"""
        Initialize the command.
        """
        self.options = options
        # shift sys.argv for compatibility with the old sage bash script and python command when consuming arguments from the command line
        import sys
        del sys.argv[0]

    def run(self) -> int:
        r"""
        Execute the given command.
        """
        input_file = self.options.file[0]
        if input_file.endswith('.sage'):
            input_file = str(preparse_file_named(input_file))
        if input_file.endswith('.pyx') or input_file.endswith('.spyx'):
            s = load_cython(input_file)
            eval(compile(s, tmp_filename(), 'exec'), sage_globals())
        else:
            eval(compile(open(input_file, 'rb').read(), input_file, 'exec'), sage_globals())
        return 0
