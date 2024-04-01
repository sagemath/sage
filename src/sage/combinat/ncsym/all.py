r"""
Symmetric functions in non-commuting variables

- :class:`Introduction to Symmetric Functions in Non-Commuting Variables <sage.combinat.ncsym.ncsym.SymmetricFunctionsNonCommutingVariables>`

- :ref:`sage.combinat.ncsym.bases`
- :ref:`sage.combinat.ncsym.dual`
- :ref:`sage.combinat.ncsym.ncsym`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.ncsym.ncsym', 'SymmetricFunctionsNonCommutingVariables')
lazy_import('sage.combinat.ncsym.dual', 'SymmetricFunctionsNonCommutingVariablesDual')

del install_doc
del lazy_import
