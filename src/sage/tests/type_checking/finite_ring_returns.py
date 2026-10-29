from __future__ import annotations

from typing import TYPE_CHECKING, assert_type

from sage.rings.polynomial.polynomial_element import Polynomial

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer

if TYPE_CHECKING:
    from sage.rings.finite_rings.element_base import FinitePolyExtElement
    from sage.rings.finite_rings.integer_mod import IntegerMod_abstract

F = GF(3**3, 'a')

assert_type(F.cardinality(), Integer)
assert_type(F.characteristic(), Integer)
assert_type(F.degree(), Integer)
assert_type(F.order(), Integer)
assert_type(F.unit_group_exponent(), Integer)
assert_type(F.zeta_order(), Integer)
assert_type(F.ngens(), int)
assert_type(F.is_conway(), bool)
assert_type(F.is_finite(), bool)
assert_type(F.is_prime_field(), bool)
assert_type(F.modulus(), Polynomial)
assert_type(F.polynomial(), Polynomial)
assert_type(F.dual_basis(), list)
assert_type(F.embeddings(F), list)
assert_type(F.some_elements(), list)
assert_type(F.construction(), tuple)

a: FinitePolyExtElement = F.gen()
assert_type(a.additive_order(), Integer)
assert_type(a.multiplicative_order(), Integer)
assert_type(a.minimal_polynomial(), Polynomial)
assert_type(a.is_square(), bool)

Z11 = GF(11)
e: IntegerMod_abstract = Z11(3)
assert_type(e.additive_order(), Integer)
assert_type(e.is_square(), bool)
assert_type(e.lift_centered(), Integer)
assert_type(e.log(Z11(2)), Integer)
assert_type(e.multiplicative_order(), Integer)
assert_type(e.polynomial(), Polynomial)
