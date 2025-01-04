from collections.abc import Iterable
from typing import Any

def normalize_index(arg: Any, nops: int, err_msg: str) -> int:
    ...

def normalize_index_for_doctests(arg: int, nops: int) -> int:
    ...

class OperandsWrapper:
    def __getitem__(self, arg: int | slice | list[int] | tuple[int, ...]) -> Any:
        ...

    def _repr_(self) -> str:
        ...

    def _latex_(self) -> str:
        ...

    def __reduce__(self) -> tuple:
        ...

def restore_op_wrapper(expr: Any) -> OperandsWrapper:
    ...
