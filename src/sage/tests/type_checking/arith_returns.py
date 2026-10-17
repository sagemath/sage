from typing import Any, assert_type

from sage.arith.misc import (
    CRT_list,
    crt,
    factor,
    inverse_mod,
    is_square,
    power_mod,
    xgcd,
    is_pseudoprime,
    jacobi_symbol,
    kronecker_symbol,
    legendre_symbol,
    next_prime,
    previous_prime,
    prime_divisors,
    primitive_root,
    random_prime,
)
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.factorization import Factorization

assert_type(is_pseudoprime(91), bool)
assert_type(next_prime(10), Integer)
assert_type(previous_prime(10), Integer)
assert_type(random_prime(30), Integer)
assert_type(factor(12), Factorization)
assert_type(prime_divisors(12), list)
assert_type(kronecker_symbol(2, 15), Integer)
assert_type(legendre_symbol(2, 7), int)
assert_type(jacobi_symbol(2, 15), Integer)
assert_type(primitive_root(17), Integer)


# overload paths for the input-polymorphic helpers (Sage 10.9 runtime-verified)
R = PolynomialRing(QQ, 'x')
x = R.gen()

assert_type(xgcd(12, 8), tuple[Integer, Integer, Integer])
assert_type(xgcd(Integer(12), 8), tuple[Integer, Integer, Integer])
assert_type(inverse_mod(3, 10), Integer)
assert_type(power_mod(3, 4, 10), int)
assert_type(power_mod(Integer(3), 4, 10), Integer)
assert_type(crt(2, 1, 3, 5), Integer)
assert_type(CRT_list([2, 1], [3, 5]), Integer | int)
assert_type(is_square(16), bool)
assert_type(is_square(16, root=True), tuple[bool, Integer])
assert_type(xgcd(x ** 2 + 1, x + 1), tuple[Any, Any, Any])
assert_type(inverse_mod(x, x ** 2 + 1), Any)
