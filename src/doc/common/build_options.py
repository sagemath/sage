###############################################
# Options for building the Sage documentation #
###############################################

import os, re
SAGE_DOC = os.environ['SAGE_DOC']
LANGUAGES = [d for d in os.listdir(SAGE_DOC) if re.match('^[a-z][a-z]$', d)]
SPHINXOPTS = ""
PAPER = ""
OMIT = ["introspect"]  # docs/dirs to omit when listing and building 'all'

if PAPER:
    PAPEROPTS = "-D latex_paper_size=" + PAPER
else:
    PAPEROPTS = ""

#Note that this needs to have the doctrees dir
ALLSPHINXOPTS   = SPHINXOPTS + " " + PAPEROPTS + " "
WEBSITESPHINXOPTS = ""

# Number of threads to use for parallel-building the documentation.
NUM_THREADS = int(os.environ.get('SAGE_NUM_THREADS', 1))

# Minimize GAP/libGAP RAM usage in the builder, docbuild already uses too much
from sage.interfaces.gap import set_gap_memory_pool_size
set_gap_memory_pool_size(80 * 1024 * 1024)

INCREMENTAL_BUILD = os.path.exists(os.path.join(SAGE_DOC, 'output'))

# Error out on errors
ABORT_ON_ERROR = True
