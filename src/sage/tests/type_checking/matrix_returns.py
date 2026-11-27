from typing import Any, assert_type

from sage.modules.free_module_element import FreeModuleElement
from sage.rings.polynomial.polynomial_element import Polynomial

from sage.matrix.constructor import matrix
from sage.matrix.matrix0 import Matrix as Matrix0
from sage.matrix.matrix1 import Matrix as Matrix1
from sage.matrix.matrix2 import Matrix
from sage.modules.free_module import FreeModule_generic
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ

assert_type(matrix(QQ, 2, 2), Matrix)

m: Matrix = matrix(QQ, [[2, 1], [1, 1]])

# matrix0-defined accessors
assert_type(m.change_ring(QQ), Matrix0)
assert_type(m.dimensions(), tuple[int, int])
assert_type(m.dict(), dict[tuple[int, int], Any])
assert_type(m.ncols(), int)
assert_type(m.rank(), Integer)
assert_type(m.nullity(), Integer)
assert_type(m.kernel(), FreeModule_generic)
assert_type(m.pivots(), tuple[int, ...])
assert_type(m.nonzero_positions(), list[tuple[int, int]])
assert_type(m.str(), str)

# matrix1-defined accessors
assert_type(m.augment(m), Matrix1)
assert_type(m.column(0), FreeModuleElement)
assert_type(m.delete_columns([0]), Matrix1)
assert_type(m.stack(m), Matrix1)
assert_type(m.submatrix(0, 0, 1, 1), Matrix1)

# matrix2-defined accessors
assert_type(m.charpoly(), Polynomial)
assert_type(m.minpoly(), Polynomial)
assert_type(m.denominator(), Integer)
assert_type(m.echelon_form(), Matrix)
assert_type(m.eigenmatrix_right(), tuple[Matrix, Matrix])
assert_type(m.inverse(), Matrix)
assert_type(m.is_diagonalizable(), bool)
assert_type(m.is_one(), bool)
assert_type(m.left_kernel(), FreeModule_generic)
assert_type(m.left_nullity(), int)
assert_type(m.rref(), Matrix)
assert_type(m.solve_right(m), FreeModuleElement | Matrix)
