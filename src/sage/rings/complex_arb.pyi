from typing import Any

from sage.libs.flint.acb import acb_t
from sage.rings.complex_interval import ComplexIntervalFieldElement
from sage.rings.real_arb import RealBall
from sage.structure.element import RingElement
from sage.rings.ring import Field

def ComplexIntervalFieldElement_to_acb(target: acb_t, source: ComplexIntervalFieldElement) -> None:
    ...

def acb_to_ComplexIntervalFieldElement(target: ComplexIntervalFieldElement, source: acb_t) -> int:
    ...

class ComplexBall(RingElement):
    value: acb_t

    def _new(self) -> 'ComplexBall':
        ...

    def _add_(self, other: Any) -> Any:
        ...

    def _mul_(self, other: Any) -> Any:
        ...

    def _complex_mpfi_(self, parent: Any) -> ComplexIntervalFieldElement:
        ...

    def real(self) -> RealBall:
        ...

    def imag(self) -> RealBall:
        ...

    def pow(self, expo: Any, analytic: Any = None) -> Any:
        ...
