r"""
Combinatorics on matrix features that are imported by default in the interpreter namespace
"""
from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.matrices.latin',
            ['LatinSquare', 'LatinSquare_generator'])
lazy_import('sage.combinat.matrices.dlxcpp', 'DLXCPP')
lazy_import('sage.combinat.matrices.hadamard_matrix',
            ['hadamard_matrix', 'hadamard_matrix_www'])

del absolute_import
