from typing import Any
from sage.structure.richcmp import richcmp
from sage.structure.element import FieldElement
from sage.rings.function_field.element import FunctionFieldElement


class FunctionFieldElement_polymod(FunctionFieldElement):
    def __init__(self, parent: Any, x: Any, reduce: bool = True) -> None:
        ...

    def element(self) -> Any:
        ...

    def _repr_(self) -> str:
        ...

    def __bool__(self) -> bool:
        ...

    def __hash__(self) -> int:
        ...

    def _richcmp_(self, other: Any, op: int) -> Any:
        ...

    def _add_(self, right: Any) -> Any:
        ...

    def _sub_(self, right: Any) -> Any:
        ...

    def _mul_(self, right: Any) -> Any:
        ...

    def _div_(self, right: Any) -> Any:
        ...

    def __invert__(self) -> Any:
        ...

    def list(self) -> list:
        ...

    def nth_root(self, n: int) -> FunctionFieldElement:
        ...

    def is_nth_power(self, n: int) -> bool:
        ...

    def _pth_root(self) -> FunctionFieldElement:
        ...
