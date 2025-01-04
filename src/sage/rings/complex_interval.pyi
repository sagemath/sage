from typing import Any

from sage.libs.mpfi.types import mpfi_t
from sage.libs.mpfr.types import mpfr_prec_t
from sage.rings.real_mpfi import RealIntervalFieldElement, RealIntervalField_class
from sage.structure.element import FieldElement

class ComplexIntervalFieldElement(FieldElement):
    __re: mpfi_t
    __im: mpfi_t
    _prec: mpfr_prec_t
    _multiplicative_order: Any

    def _new(self) -> 'ComplexIntervalFieldElement':
        ...

    def _new_real(self) -> RealIntervalFieldElement:
        ...
