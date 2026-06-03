from __future__ import annotations

from collections.abc import Sequence
from types import TracebackType

from sage.structure.category_object import NameSpec
from sage.structure.parent import Parent

class localvars:
    _obj: Parent
    _names: Sequence[str]
    _latex_names: Sequence[str] | None
    _orig: tuple[Sequence[str] | None, Sequence[str] | None] | None

    def __init__(
        self,
        obj: ParentWithGens,
        names: NameSpec,
        latex_names: Sequence[str] | None = None,
        normalize: bool = True,
    ) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None: ...
