from collections.abc import Callable
from typing import Union

def abc(f: Callable | None = None, optional: bool = False) -> Callable:
    ...

class ABC:
    def __init__(self, f: Callable, optional: bool = False) -> None:
        ...

    def __repr__(self) -> str:
        ...

    def _sage_src_lines_(self) -> Union[str, int]:
        ...

    def __get__(self, instance: object, cls: type) -> Union[Callable, NotImplementedError]:
        ...

    def is_optional(self) -> bool:
        ...

def abstract_methods_of_class(cls: type) -> dict[str, list[str]]:
    ...
