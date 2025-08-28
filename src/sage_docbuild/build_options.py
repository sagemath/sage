r"""
Build options

This module defines options for building Sage documentation.
"""

import argparse
import os
from pathlib import Path

SPHINXOPTS = ""
PAPER = ""
OMIT = ["introspect"]  # docs/dirs to omit when listing and building 'all'

if PAPER:
    PAPEROPTS = "-D latex_paper_size=" + PAPER
else:
    PAPEROPTS = ""

# Note that this needs to have the doctrees dir
ALLSPHINXOPTS = SPHINXOPTS + " " + PAPEROPTS + " "
WEBSITESPHINXOPTS = ""

# Number of threads to use for parallel-building the documentation.
NUM_THREADS = int(os.environ.get('SAGE_NUM_THREADS', 1))

# Error out on errors
ABORT_ON_ERROR = True

class BuildOptions(argparse.Namespace):
    source_dir: Path
    output_dir: Path
