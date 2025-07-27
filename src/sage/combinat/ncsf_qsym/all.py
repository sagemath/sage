r"""
Noncommutative symmetric functions and quasi-symmetric functions

- :ref:`sage.combinat.ncsf_qsym.tutorial`

- :ref:`Non-Commutative Symmetric Functions (NCSF) <sage.combinat.ncsf_qsym.ncsf>`
- :ref:`Quasi-Symmetric Functions (QSym) <sage.combinat.ncsf_qsym.qsym>`
- :ref:`sage.combinat.ncsf_qsym.generic_basis_code`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.ncsf_qsym.qsym', 'QuasiSymmetricFunctions')
lazy_import('sage.combinat.ncsf_qsym.ncsf', 'NonCommutativeSymmetricFunctions')

del install_doc
del lazy_import
