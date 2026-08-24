from typing import Any

from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteFieldBase


class FiniteFieldFactory:
    def __call__(self, *args: Any, **kwds: Any) -> FiniteFieldBase: ...


GF: FiniteFieldFactory
FiniteField: FiniteFieldFactory
