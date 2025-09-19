from collections.abc import Callable
from typing import Any
from sage.symbolic.expression import Expression

class E(Expression):
    def __init__(self) -> None:
        ...

    def __pow__(self, left: Any, right: Any, dummy: Any) -> Any:
        ...
