from sage.misc.lazy_import import lazy_import
from .matrix_space import MatrixSpace
from .constructor import (Matrix, column_matrix, random_matrix,
                         diagonal_matrix, identity_matrix, block_matrix,
                         block_diagonal_matrix, jordan_block, zero_matrix,
                         ones_matrix, elementary_matrix, companion_matrix)

# Importing matrix from constructor module directly
from .constructor import matrix as original_matrix

def matrix(*args, **kwargs):
    if 'base_ring' in kwargs:
        return original_matrix(*args, base_ring=kwargs['base_ring'])
    else:
        return original_matrix(*args, **kwargs)

Mat = MatrixSpace
