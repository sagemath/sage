from typing import assert_type, cast

from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.ring import Ring
from sage.structure.element import (
    CommutativeRingElement,
    Element,
    PrincipalIdealDomainElement,
)

e = cast(Element, Integer(0))
assert_type(e.is_zero(), bool)

r = cast(CommutativeRingElement, QQ(3))
assert_type(r.inverse_mod(QQ(7)), CommutativeRingElement)

p = cast(PrincipalIdealDomainElement, Integer(12))
assert_type(p.gcd(18), Element)
assert_type(p.lcm(18), Element)

Z = cast(Ring, ZZ)
assert_type(Z.order(), Integer | PlusInfinity)
