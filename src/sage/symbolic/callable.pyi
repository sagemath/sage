from collections.abc import Callable
from typing import Any, Union
from sage.rings.abc import CallableSymbolicExpressionRing as CallableSymbolicExpressionRingABC
from sage.symbolic.ring import SymbolicRing, SR
from sage.categories.pushout import ConstructionFunctor
from sage.structure.factory import UniqueFactory

class CallableSymbolicExpressionFunctor(ConstructionFunctor):
    def __init__(self, arguments: tuple):
        ...

    def __repr__(self) -> str:
        ...

    def merge(self, other: 'CallableSymbolicExpressionFunctor') -> 'CallableSymbolicExpressionFunctor':
        ...

    def __call__(self, R: SymbolicRing) -> 'CallableSymbolicExpressionRing':
        ...

    def arguments(self) -> tuple:
        ...

    def unify_arguments(self, x: 'CallableSymbolicExpressionFunctor') -> tuple:
        ...

class CallableSymbolicExpressionRing_class(SymbolicRing, CallableSymbolicExpressionRingABC):
    def __init__(self, arguments: tuple):
        ...

    def _coerce_map_from_(self, R: Any) -> bool:
        ...

    def construction(self) -> tuple:
        ...

    def _element_constructor_(self, x: Any) -> Any:
        ...

    def _repr_(self) -> str:
        ...

    def arguments(self) -> tuple:
        ...

    args = arguments

    def _repr_element_(self, x: Any) -> str:
        ...

    def _latex_element_(self, x: Any) -> str:
        ...

    def _call_element_(self, _the_element: Any, *args: Any, **kwds: Any) -> Any:
        ...

    __reduce__ = object.__reduce__

class CallableSymbolicExpressionRingFactory(UniqueFactory):
    def create_key(self, args: Any, check: bool = True) -> tuple:
        ...

    def create_object(self, version: int, key: tuple, **extra_args: Any) -> 'CallableSymbolicExpressionRing_class':
        ...

CallableSymbolicExpressionRing = CallableSymbolicExpressionRingFactory('sage.symbolic.callable.CallableSymbolicExpressionRing')
