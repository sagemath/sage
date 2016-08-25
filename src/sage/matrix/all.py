from __future__ import absolute_import
from .matrix_space import MatrixSpace
from .constructor import (matrix, Matrix, column_matrix, random_matrix,
                         diagonal_matrix, identity_matrix, block_matrix,
                         block_diagonal_matrix, jordan_block, zero_matrix,
                         ones_matrix, elementary_matrix, companion_matrix)
from .berlekamp_massey import berlekamp_massey

Mat = MatrixSpace
