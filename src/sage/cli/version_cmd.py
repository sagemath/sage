import argparse

from sage.version import version


class VersionCmd:
    @staticmethod
    def extend_parser(parser: argparse.ArgumentParser):
        r"""
        Extend the parser with the version command.

        INPUT:

        - ``parsers`` -- the parsers to extend.

        OUTPUT:

        - the extended parser.
        """
        parser.add_argument(
            "--version",
            action="version",
            version=version,
            help="print the version number and exit",
        )
