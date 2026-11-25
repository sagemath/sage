from __future__ import annotations

from sage.version import version
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import argparse


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
            "-V",
            "--version",
            action="version",
            version=version,
            help="print the version number and exit",
        )
