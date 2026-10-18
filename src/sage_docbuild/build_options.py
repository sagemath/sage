r"""
Build options

This module defines options for building Sage documentation.
"""

import argparse
import os
from pathlib import Path

SPHINXOPTS: list[str] = []
PAPER = ""
OMIT = ["introspect"]  # docs/dirs to omit when listing and building 'all'

if PAPER:
    PAPEROPTS = ["-D", "latex_paper_size=" + PAPER]
else:
    PAPEROPTS: list[str] = []

# Options passed to every Sphinx invocation, as a list of arguments: they are
# handed to sphinx-build as such, so an option or a path may contain spaces.
ALLSPHINXOPTS: list[str] = SPHINXOPTS + PAPEROPTS
WEBSITESPHINXOPTS: list[str] = []

# Number of threads to use for parallel-building the documentation.
NUM_THREADS = int(os.environ.get('SAGE_NUM_THREADS', 1))

# Error out on errors
ABORT_ON_ERROR = True

# Run Sphinx in nitpicky mode (-n) to warn about unresolved links; this is
# applied per build, skipping inventory builds where cross-references are
# not yet resolvable.
WARN_LINKS = False

class BuildOptions(argparse.Namespace):
    source_dir: Path
    output_dir: Path
    #: Whether ``output_dir`` was given on the command line, as opposed to
    #: being filled in with the default. Single-file builds use their own
    #: default, so that they never write into the installed documentation.
    output_dir_given: bool
