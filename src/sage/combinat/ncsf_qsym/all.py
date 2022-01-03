r"""
Non-commutative symmetric functions and quasi-symmetric functions

- :ref:`sage.combinat.ncsf_qsym.tutorial`

- :ref:`Non-Commutative Symmetric Functions (NCSF) <sage.combinat.ncsf_qsym.ncsf>`
- :ref:`Quasi-Symmetric Functions (QSym) <sage.combinat.ncsf_qsym.qsym>`
- :ref:`sage.combinat.ncsf_qsym.generic_basis_code`

"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from .qsym import QuasiSymmetricFunctions
from .ncsf import NonCommutativeSymmetricFunctions

