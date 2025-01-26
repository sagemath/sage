r"""
Symmetric functions

- :class:`Introduction to Symmetric Functions <sage.combinat.sf.sf.SymmetricFunctions>`

- :ref:`sage.combinat.sf.sfa`
- :ref:`sage.combinat.sf.sf`
- :ref:`sage.combinat.sf.classical`
- :ref:`sage.combinat.sf.schur`
- :ref:`sage.combinat.sf.monomial`
- :ref:`sage.combinat.sf.multiplicative`
- :ref:`sage.combinat.sf.elementary`
- :ref:`sage.combinat.sf.homogeneous`
- :ref:`sage.combinat.sf.powersum`
- :ref:`sage.combinat.sf.character`
- :ref:`sage.combinat.sf.orthogonal`
- :ref:`sage.combinat.sf.symplectic`
- :ref:`sage.combinat.sf.dual`
- :ref:`sage.combinat.sf.orthotriang`
- :ref:`sage.combinat.sf.kfpoly`
- :ref:`sage.combinat.sf.hall_littlewood`
- :ref:`sage.combinat.sf.hecke`
- :ref:`sage.combinat.sf.jack`
- :ref:`k-Schur Functions <sage.combinat.sf.new_kschur>`
- :ref:`sage.combinat.sf.k_dual`
- :ref:`sage.combinat.sf.llt`
- :ref:`sage.combinat.sf.macdonald`
- :ref:`sage.combinat.sf.ns_macdonald`
- :ref:`sage.combinat.sf.witt`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

# In the long run, this will be the single entry point
# Nothing else will be exported
lazy_import('sage.combinat.sf.sf', 'SymmetricFunctions')

# Advanced stuff:

lazy_import('sage.combinat.sf.kfpoly', 'KostkaFoulkesPolynomial')

lazy_import('sage.combinat.sf.ns_macdonald', ['NonattackingFillings',
                                              'AugmentedLatticeDiagramFilling',
                                              'LatticeDiagram'])
