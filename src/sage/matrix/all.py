from sage.misc.lazy_import import lazy_import
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import (matrix, Matrix, column_matrix, random_matrix,
                                     diagonal_matrix, identity_matrix, block_matrix,
                                     block_diagonal_matrix, jordan_block, zero_matrix,
                                     ones_matrix, elementary_matrix, companion_matrix)
Mat = MatrixSpace
del lazy_import
