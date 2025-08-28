from typing import Any

from sage.libs.gsl.types import gsl_complex
from sage.structure.element import FieldElement

class ComplexDoubleElement(FieldElement):
    _complex: gsl_complex

    def _new_c(self, x: gsl_complex) -> 'ComplexDoubleElement':
        ...

    def _add_(self, other: Any) -> Any:
        ...

    def _mul_(self, other: Any) -> Any:
        ...

    def _pow_(self, other: Any) -> Any:
        ...

def new_ComplexDoubleElement() -> ComplexDoubleElement:
    ...
