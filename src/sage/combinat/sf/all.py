r"""
Symmetric function features that are imported by default in the interpreter namespace
"""
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

# In the long run, this will be the single entry point
# Nothing else will be exported
lazy_import('sage.combinat.sf.sf', 'SymmetricFunctions')

# Advanced stuff:

lazy_import('sage.combinat.sf.kfpoly', 'KostkaFoulkesPolynomial')

lazy_import('sage.combinat.sf.ns_macdonald', ['NonattackingFillings',
                                              'AugmentedLatticeDiagramFilling',
                                              'LatticeDiagram'])

del absolute_import
