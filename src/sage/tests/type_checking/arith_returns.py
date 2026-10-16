from typing import assert_type

from sage.arith.misc import (
    factor,
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
