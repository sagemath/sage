from __future__ import annotations

from typing import TYPE_CHECKING, Any, assert_type, cast

from sage.matrix.constructor import matrix
from sage.matrix.matrix_dense import Matrix_dense
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ, crt_basis
from sage.rings.rational_field import QQ
from sage.rings.ring import Ring
from sage.structure.parent import Parent

if TYPE_CHECKING:
    from sage.structure.category_object import CategoryObject

R: CategoryObject = QQ['x']
P = cast("Parent", QQ)

assert_type(R.base_ring(), Ring)
assert_type(R.variable_name(), str)
assert_type(R.variable_names(), tuple[str, ...])
assert_type(R.latex_name(), str)
assert_type(R.latex_variable_names(), list[str])
assert_type(R.categories(), list[Any])
assert_type(R.gens_dict(), dict[str, Any])
assert_type(R.gens_dict_recursive(), dict[str, Any])
assert_type(R.inject_variables(), None)
assert_type(R.objgen(), tuple[Any, Any])
assert_type(R.objgens(), tuple[Any, Any])

# register_action/coercion/conversion/embedding are annotated -> None in the
# stub but not callable here: registration is only legal on fresh parents.
assert_type(P, Parent)


assert_type(crt_basis([ZZ(3), ZZ(5)]), list[Integer])

mi = cast("Matrix_integer_dense", matrix(ZZ, [[1, 2], [3, 4]]))
assert_type(mi.LLL(), Matrix_integer_dense)
md = cast("Matrix_dense", mi)
assert_type(md.transpose(), Matrix_dense)
