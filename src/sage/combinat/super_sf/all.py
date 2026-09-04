r"""
Supersymmetric Functions

- :ref:`sage.combinat.super_sf.super_sf`
- :ref:`sage.combinat.super_sf.super_sfa`
- :ref:`sage.combinat.super_sf.powersum`
- :ref:`sage.combinat.super_sf.hom_el`
- :ref:`sage.combinat.super_sf.schur`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc

install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

# In the long run, this will be the single entry point
# Nothing else will be exported
lazy_import('sage.combinat.super_sf.super_sf', 'SupersymmetricFunctions')
