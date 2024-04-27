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

from .ncsym import SymmetricFunctionsNonCommutingVariables
from .dual import SymmetricFunctionsNonCommutingVariablesDual
