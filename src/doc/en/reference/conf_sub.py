# nodoctest
# Sage documentation build configuration file for a sub-document of the
# reference manual.
#
# The conf.py of every sub-document is a symbolic link to this file, so keep
# here only what all of them share.  The values that depend on the individual
# sub-document are computed from its directory by reference_subdocument();
# see sage_docbuild.conf.

from sage_docbuild.conf import *
from sage_docbuild.conf import reference_subdocument

globals().update(reference_subdocument())
