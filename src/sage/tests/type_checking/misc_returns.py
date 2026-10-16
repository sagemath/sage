from typing import assert_type

from sage.misc.functional import multiplicative_order
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.infinity import PlusInfinity


assert_type(multiplicative_order(GF(11)(3)), Integer | PlusInfinity)
assert_type(multiplicative_order(QQ(2 / 3)), Integer | PlusInfinity)

# is_exact is annotated in source but @cached_method-wrapped; ty cannot see
# through the decorator yet, so it is not assert_type-checked here.
