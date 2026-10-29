from typing import assert_type

from sage.rings.integer import Integer
from sage.rings.rational import Rational

q = Rational(22, 7)

assert_type(q.absolute_norm(), Rational)
assert_type(q.ceil(), Integer)
assert_type(q.denominator(), Integer)
assert_type(q.floor(), Integer)
assert_type(q.height(), Integer)
assert_type(q.imag(), Rational)
assert_type(q.is_integral(), bool)
assert_type(q.is_one(), bool)
assert_type(q.is_perfect_power(), bool)
assert_type(q.is_square(), bool)
assert_type(q.list(), list[Rational])
assert_type(q.numerator(), Integer)
assert_type(q.period(), Integer)
assert_type(q.real(), Rational)
assert_type(q.round(), Integer)
assert_type(q.sign(), Integer)
assert_type(q.squarefree_part(), Integer)
assert_type(q.str(), str)
assert_type(q.support(), list[Integer])
assert_type(q.trunc(), Integer)
assert_type(q.conjugate(), Rational)
assert_type(q.norm(), Rational)
assert_type(Rational(4).nth_root(2), Rational)
assert_type(q.as_integer_ratio(), tuple[int, int])
assert_type(q.mod_ui(5), int)
assert_type(q.is_S_unit(), bool)
