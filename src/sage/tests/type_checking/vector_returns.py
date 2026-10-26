from collections.abc import Iterator
from typing import Any, assert_type

from sage.matrix.matrix0 import Matrix

from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.rings.ring import Ring

v: FreeModuleElement = vector(QQ, [1, 2])

assert_type(vector(QQ, [1, 2]), FreeModuleElement)
assert_type(v.additive_order(), Integer)
assert_type(v.change_ring(QQ), FreeModuleElement)
assert_type(v.column(), Matrix)
assert_type(v.conjugate(), FreeModuleElement)
assert_type(v.coordinate_ring(), Ring)
assert_type(v.degree(), int)
assert_type(v.denominator(), Integer)
assert_type(v.dense_vector(), FreeModuleElement)
assert_type(v.dict(), dict[int, Any])
assert_type(v.items(), Iterator[tuple[int, Any]])
assert_type(v.list(), list[Any])
assert_type(v.nonzero_positions(), list[int])
assert_type(v.outer_product(v), Matrix)
assert_type(v.row(), Matrix)
assert_type(v.sparse_vector(), FreeModuleElement)
assert_type(v.support(), list[int])
assert_type(v.Mod(5), FreeModuleElement)
assert_type(v.apply_map(lambda c: c), FreeModuleElement)
assert_type(v.concatenate(v), FreeModuleElement)
assert_type(v.pairwise_product(v), FreeModuleElement)
