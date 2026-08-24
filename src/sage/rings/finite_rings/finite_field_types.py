from typing import assert_type

from sage.rings.finite_rings.element_base import FiniteRingElement
from sage.rings.finite_rings.finite_field_constructor import GF


for F in (GF(29), GF(2**8, 'a')):
    assert_type(F.from_integer(1), FiniteRingElement)
    assert_type(F.gen(), FiniteRingElement)
    assert_type(F.multiplicative_generator(), FiniteRingElement)
    assert_type(F.random_element(), FiniteRingElement)
