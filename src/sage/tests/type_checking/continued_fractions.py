from typing import assert_type

from sage.rings.continued_fraction import (
    ContinuedFraction_infinite,
    ContinuedFraction_periodic,
    ContinuedFraction_real,
    continued_fraction,
    convergents,
)
from sage.rings.rational import Rational


assert_type(continued_fraction(22 / 7), ContinuedFraction_periodic | ContinuedFraction_real | ContinuedFraction_infinite)
assert_type(convergents(22 / 7), list[Rational])

cf = continued_fraction([1, 2, 3])
assert_type(cf.convergents(), list[Rational])
