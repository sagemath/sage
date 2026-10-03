from typing import Any

from sage.rings.finite_rings.element_base import FiniteRingElement
from sage.rings.ring import Field


class FiniteField(Field):
    def from_integer(self, n: Any, reverse: bool = False) -> FiniteRingElement: ...

    def gen(self) -> FiniteRingElement: ...

    def multiplicative_generator(self) -> FiniteRingElement: ...

    def random_element(self, *args: Any, **kwds: Any) -> FiniteRingElement: ...
