from typing import Any, Dict, Tuple, Union

def is_FpTElement(x: Any) -> bool:
    ...

class FpTElement:
    _numer: Any
    _denom: Any
    initialized: bool
    p: int

    def __init__(self, parent: Any, numer: Any, denom: Any = 1, coerce: bool = True, reduce: bool = True) -> None:
        ...

    def _new_c(self) -> 'FpTElement':
        ...

    def _copy_c(self) -> 'FpTElement':
        ...

    def numer(self) -> Any:
        ...

    def numerator(self) -> Any:
        ...

    def denom(self) -> Any:
        ...

    def denominator(self) -> Any:
        ...

    def __call__(self, *args: Any, **kwds: Any) -> 'FpTElement':
        ...

    def subs(self, in_dict: Dict = None, *args: Any, **kwds: Any) -> 'FpTElement':
        ...

    def valuation(self, v: Any) -> int:
        ...

    def factor(self) -> Any:
        ...

    def _repr_(self) -> str:
        ...

    def _latex_(self) -> str:
        ...

    def _richcmp_(self, other: 'FpTElement', op: int) -> bool:
        ...

    def __hash__(self) -> int:
        ...

    def __neg__(self) -> 'FpTElement':
        ...

    def __invert__(self) -> 'FpTElement':
        ...

    def _add_(self, other: 'FpTElement') -> 'FpTElement':
        ...

    def _sub_(self, other: 'FpTElement') -> 'FpTElement':
        ...

    def _mul_(self, other: 'FpTElement') -> 'FpTElement':
        ...

    def _div_(self, other: 'FpTElement') -> 'FpTElement':
        ...

    def next(self) -> 'FpTElement':
        ...

    def _sqrt_or_None(self) -> 'FpTElement | None':
        ...

    def is_square(self) -> bool:
        ...

    def sqrt(self, extend: bool = True, all: bool = False) -> 'FpTElement | list[FpTElement]':
        ...

    def __pow__(self, e: int, dummy: Any) -> 'FpTElement':
        ...

class FpT:
    INTEGER_LIMIT: int

    def __init__(self, R: Any, names: Any = None) -> None:
        ...

    def __iter__(self) -> 'FpT_iter':
        ...

    def iter(self, bound: Any = None, start: Any = None) -> 'FpT_iter':
        ...

class FpT_iter:
    parent: Any
    degree: int
    cur: FpTElement
    g: Any

    def __init__(self, parent: Any, degree: int = None, start: FpTElement = None) -> None:
        ...

    def __cinit__(self, parent: Any, *args: Any, **kwds: Any) -> None:
        ...

    def __dealloc__(self) -> None:
        ...

    def __iter__(self) -> 'FpT_iter':
        ...

    def __next__(self) -> FpTElement:
        ...

class Polyring_FpT_coerce:
    p: int

    def __init__(self, R: Any) -> None:
        ...

    def _extra_slots(self) -> Dict:
        ...

    def _update_slots(self, _slots: Dict) -> None:
        ...

    def _call_(self, _x: Any) -> FpTElement:
        ...

    def _call_with_args(self, _x: Any, args: Tuple = (), kwds: Dict = {}) -> FpTElement:
        ...

    def section(self) -> 'FpT_Polyring_section':
        ...

class FpT_Polyring_section:
    p: int

    def __init__(self, f: Polyring_FpT_coerce) -> None:
        ...

    def _extra_slots(self) -> Dict:
        ...

    def _update_slots(self, _slots: Dict) -> None:
        ...

    def _call_(self, _x: Any) -> Any:
        ...

class Fp_FpT_coerce:
    p: int

    def __init__(self, R: Any) -> None:
        ...

    def _extra_slots(self) -> Dict:
        ...

    def _update_slots(self, _slots: Dict) -> None:
        ...

    def _call_(self, _x: Any) -> FpTElement:
        ...

    def _call_with_args(self, _x: Any, args: Tuple = (), kwds: Dict = {}) -> FpTElement:
        ...

    def section(self) -> 'FpT_Fp_section':
        ...

class FpT_Fp_section:
    p: int

    def __init__(self, f: Fp_FpT_coerce) -> None:
        ...

    def _extra_slots(self) -> Dict:
        ...

    def _update_slots(self, _slots: Dict) -> None:
        ...

    def _call_(self, _x: Any) -> Any:
        ...

class ZZ_FpT_coerce:
    p: int

    def __init__(self, R: Any) -> None:
        ...

    def _extra_slots(self) -> Dict:
        ...

    def _update_slots(self, _slots: Dict) -> None:
        ...

    def _call_(self, _x: Any) -> FpTElement:
        ...

    def _call_with_args(self, _x: Any, args: Tuple = (), kwds: Dict = {}) -> FpTElement:
        ...

    def section(self) -> Any:
        ...

def unpickle_FpT_element(K: Any, numer: Any, denom: Any) -> FpTElement:
    ...
