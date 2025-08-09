from typing import Any

from sage.libs.mpfr.types import mpfr_t, mpfr_prec_t
from sage.structure.element import FieldElement
from sage.rings.real_mpfr import RealNumber

class ComplexNumber(FieldElement):
    __re: mpfr_t
    __im: mpfr_t
    _prec: mpfr_prec_t
    _multiplicative_order: Any

    def _add_(self, other: Any) -> Any:
        ...

    def _mul_(self, other: Any) -> Any:
        ...

    def abs_c(self) -> RealNumber:
        ...

    def norm_c(self) -> RealNumber:
        ...

    def _new(self) -> 'ComplexNumber':
        ...
