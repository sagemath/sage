from typing import Any, assert_type

from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.rational_field import QQ
from sage.structure.factorization import Factorization

R = QQ['x']
x = R.gen()
p: Polynomial = x**2 + 1

assert_type(p.degree(), Integer)
assert_type(p.denominator(), Polynomial)
assert_type(p.dispersion(), Integer)
assert_type(p.euclidean_degree(), Integer)
assert_type(p.exponents(), list[int])
assert_type(p.factor(), Factorization)
assert_type(p.gradient(), list)
assert_type(p.homogenize(), Polynomial)
assert_type(p.integral(), Polynomial)
assert_type(p.is_cyclotomic_product(), bool)
assert_type(p.is_gen(), bool)
assert_type(p.is_irreducible(), bool)
assert_type(p.is_monic(), bool)
assert_type(p.is_squarefree(), bool)
assert_type(p.is_unit(), bool)
assert_type(p.is_weil_polynomial(), bool)
assert_type(p.lcm(x + 1), Polynomial)
assert_type(p.lt(), Polynomial)
assert_type(p.monic(), Polynomial)
assert_type(p.monomials(), list)
assert_type(p.numerator(), Polynomial)
assert_type(p.number_of_real_roots(), int)
assert_type(p.padded_list(), list)
assert_type(p.roots(), list)
assert_type(p.square(), Polynomial)
assert_type(p.subresultants(x + 1), list)
assert_type(p.variable_name(), str)
assert_type(p.variables(), tuple)
assert_type(p.quo_rem(x + 1), tuple[Any, Any])
assert_type(p.xgcd(x + 1), tuple[Any, Any, Any])
