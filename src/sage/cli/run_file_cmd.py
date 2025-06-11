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
            nargs="?",
            help="execute the given file as sage code",
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
        input_file = preparse_file_named(self.options.file) if self.options.file.endswith('.sage') else self.options.file
        try:
            if self.options.file.endswith('.pyx') or self.options.file.endswith('.spyx'):
                s = load_cython(input_file)
                eval(compile(s, tmp_filename(), 'exec'), sage_globals())
            else:
                eval(compile(open(input_file, 'rb').read(), input_file, 'exec'), sage_globals())
        except Exception as e:
            print(f"An error occurred while executing the file: {e}")
            return 1
        return 0
