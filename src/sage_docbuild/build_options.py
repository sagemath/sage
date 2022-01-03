###############################################
# Options for building the Sage documentation #
###############################################

import os, re

from sage.env import SAGE_DOC_SRC, SAGE_DOC

LANGUAGES = [d for d in os.listdir(SAGE_DOC_SRC) if re.match('^[a-z][a-z]$', d)]
SPHINXOPTS = ""
PAPER = ""
OMIT = ["introspect"]  # docs/dirs to omit when listing and building 'all'

if PAPER:
    PAPEROPTS = "-D latex_paper_size=" + PAPER
else:
    PAPEROPTS = ""

#Note that this needs to have the doctrees dir
ALLSPHINXOPTS = SPHINXOPTS + " " + PAPEROPTS + " "
WEBSITESPHINXOPTS = ""

# Number of threads to use for parallel-building the documentation.
NUM_THREADS = int(os.environ.get('SAGE_NUM_THREADS', 1))

INCREMENTAL_BUILD = os.path.isdir(SAGE_DOC)

# Error out on errors
ABORT_ON_ERROR = True
