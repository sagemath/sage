from typing import assert_type

from sage.misc.functional import multiplicative_order
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.infinity import PlusInfinity
from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.rings.rational_field import QQ


assert_type(multiplicative_order(GF(11)(3)), Integer | PlusInfinity)
assert_type(multiplicative_order(QQ(2 / 3)), Integer | PlusInfinity)

R = PolynomialRing(QQ, 'x')
assert_type(R.is_exact(), bool)
