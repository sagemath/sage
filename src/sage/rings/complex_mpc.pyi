from typing import Any

from sage.libs.mpc.types import mpc_t, mpc_rnd_t
from sage.structure.element import FieldElement
from sage.rings.ring import Field

class MPComplexNumber(FieldElement):
    value: mpc_t
    init: bool

    def _new(self) -> 'MPComplexNumber':
        ...

    def _add_(self, other: Any) -> Any:
        ...

    def _mul_(self, other: Any) -> Any:
        ...

class MPComplexField_class(Field):
    _prec: int
    __rnd: mpc_rnd_t
    __rnd_str: Any
    __real_field: Any
    __imag_field: Any

    def _new(self) -> MPComplexNumber:
        ...

    def _an_element_(self) -> MPComplexNumber:
        ...
