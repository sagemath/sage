from typing import assert_type

from sage.matrix.special import (
    column_matrix,
    diagonal_matrix,
    identity_matrix,
    random_matrix,
    zero_matrix,
)
from sage.rings.integer_ring import ZZ
from sage.structure.element import Matrix


assert_type(column_matrix([1, 2, 3]), Matrix)
assert_type(diagonal_matrix([1, 2, 3]), Matrix)
assert_type(identity_matrix(ZZ, 3), Matrix)
assert_type(zero_matrix(ZZ, 2, 3), Matrix)
assert_type(random_matrix(ZZ, 2, 3), Matrix)
