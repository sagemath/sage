from __future__ import annotations

from types import TracebackType
from typing import Any, Mapping, Sequence

from sage.structure.category_object import NameSpec
from sage.structure.parent_base import ParentWithBase

class ParentWithGens(ParentWithBase):
    _base: Any
    _gens: tuple[Any, ...] | None
    _names: tuple[str, ...] | None
    _latex_names: tuple[str, ...] | None
    _list: Any

    def __init__(
        self,
        base: Any,
        names: NameSpec = ...,
        normalize: bool = ...,
        category: Any | None = ...,
    ) -> None: ...
    def ngens(self) -> int: ...
    def gen(self, i: int = ...) -> Any: ...
    def gens(self) -> tuple[Any, ...]: ...
    def _assign_names(self, names: NameSpec = ..., normalize: bool = ...) -> None: ...
    def __getstate__(self) -> dict[str, Any]: ...
    def __setstate__(self, d: Mapping[str, Any]) -> None: ...
    def hom(
        self,
        im_gens: Any,
        codomain: Any | None = ...,
        base_map: Any | None = ...,
        category: Any | None = ...,
        check: bool = ...,
    ) -> Any: ...

class localvars:
    _obj: ParentWithGens
    _names: Sequence[str]
    _latex_names: Sequence[str] | None
    _orig: tuple[Sequence[str] | None, Sequence[str] | None] | None

    def __init__(
        self,
        obj: ParentWithGens,
        names: NameSpec,
        latex_names: Sequence[str] | None = ...,
        normalize: bool = ...,
    ) -> None: ...
    def __enter__(self) -> None: ...
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None: ...
