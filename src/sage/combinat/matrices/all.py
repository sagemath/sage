# sage_setup: distribution = sagemath-combinat
r"""
Combinatorics on matrices

- :ref:`sage.combinat.matrices.dancing_links`
- :ref:`sage.combinat.matrices.dlxcpp`
- :ref:`sage.combinat.matrices.hadamard_matrix`
- :ref:`sage.combinat.matrices.latin`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.combinat.matrices.all__sagemath_categories import *

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.matrices.latin',
            ['LatinSquare', 'LatinSquare_generator'])
lazy_import('sage.combinat.matrices.hadamard_matrix',
            ['hadamard_matrix', 'hadamard_matrix_www'])

del install_doc
del lazy_import
