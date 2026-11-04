from typing import assert_type

from sage.functions.other import binomial, factorial
from sage.rings.integer import Integer

assert_type(binomial(5, 2), int)
assert_type(binomial(Integer(5), 2), Integer)
assert_type(factorial(5), int)
assert_type(factorial(Integer(5)), Integer)
