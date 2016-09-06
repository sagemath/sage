from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .matrix_space import MatrixSpace
from .constructor import (matrix, Matrix, column_matrix, random_matrix,
                         diagonal_matrix, identity_matrix, block_matrix,
                         block_diagonal_matrix, jordan_block, zero_matrix,
                         ones_matrix, elementary_matrix, companion_matrix)
from .berlekamp_massey import berlekamp_massey

lazy_import('sage.matrix.joint_spectral_radius', 'joint_spectral_radius')

Mat = MatrixSpace
