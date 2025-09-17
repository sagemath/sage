from collections.abc import Callable
from typing import Any

def add_vararg(first: Any, *rest: Any) -> Any:
    ...

def mul_vararg(first: Any, *rest: Any) -> Any:
    ...

arithmetic_operators: dict[Callable, str] = {
    add_vararg: '+',
    mul_vararg: '*',
    operator.add: '+',
    operator.sub: '-',
    operator.mul: '*',
    operator.truediv: '/',
    operator.floordiv: '//',
    operator.pow: '^'
}

relation_operators: dict[Callable, str] = {
    operator.eq: '==',
    operator.lt: '<',
    operator.gt: '>',
    operator.ne: '!=',
    operator.le: '<=',
    operator.ge: '>='
}

class FDerivativeOperator:
    def __init__(self, function: Callable, parameter_set: list[int]) -> None:
        ...

    def __call__(self, *args: Any) -> Any:
        ...

    def __repr__(self) -> str:
        ...

    def function(self) -> Callable:
        ...

    def change_function(self, new: Callable) -> FDerivativeOperator:
        ...

    def parameter_set(self) -> list[int]:
        ...

class DerivativeOperator:
    class DerivativeOperatorWithParameters:
        def __init__(self, parameter_set: list[int]) -> None:
            ...

        def __call__(self, function: Callable) -> FDerivativeOperator:
            ...

        def __repr__(self) -> str:
            ...

    def __getitem__(self, args: int | tuple[int, ...]) -> DerivativeOperatorWithParameters:
        ...
