from typing import assert_type

from sage.matrix.matrix_space import MatrixSpace
from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.sets.primes import Primes

P = Primes()

assert_type(P.first(), Integer)
assert_type(P.next(7), Integer)
assert_type(P.unrank(4), Integer)
assert_type(P.cardinality(), PlusInfinity)

M = MatrixSpace(QQ, 2, 3)
assert_type(M.matrix_space(), MatrixSpace)
assert_type(M.matrix_space(4, 5), MatrixSpace)
