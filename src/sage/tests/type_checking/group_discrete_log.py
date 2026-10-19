from typing import assert_type

from sage.groups.generic import bsgs, discrete_log_rho
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer


F = GF(11)
g = F(2)

assert_type(bsgs(g, g ** 7, (1, 10)), Integer)
assert_type(discrete_log_rho(g ** 2, g, ord=5), Integer)
