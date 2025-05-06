#!/usr/bin/env python3

import argparse
import logging
import sys

from sage.cli.eval_cmd import EvalCmd
from sage.cli.interactive_shell_cmd import InteractiveShellCmd
from sage.cli.notebook_cmd import JupyterNotebookCmd
from sage.cli.options import CliOptions
from sage.cli.version_cmd import VersionCmd


def main() -> int:
    input_args = sys.argv[1:]
    parser = argparse.ArgumentParser(
        prog="sage",
        description="If no command is given, starts the interactive interpreter where you can enter statements and expressions, immediately execute them and see their results.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="print additional information",
    )

    VersionCmd.extend_parser(parser)
    JupyterNotebookCmd.extend_parser(parser)
    EvalCmd.extend_parser(parser)

    if not input_args:
        return InteractiveShellCmd(CliOptions()).run()

    args = parser.parse_args(input_args)
    options = CliOptions(**vars(args))

    logging.basicConfig(level=logging.DEBUG if options.verbose else logging.INFO)

    if args.command:
        return EvalCmd(options).run()
    elif args.notebook:
        return JupyterNotebookCmd(options).run()
    return InteractiveShellCmd(options).run()
