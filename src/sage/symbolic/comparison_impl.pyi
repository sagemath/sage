from collections.abc import Iterable
from typing import Any

def print_order(lhs: Any, rhs: Any) -> int:
    ...

class _print_key:
    def __init__(self, ex: Any):
        ...

    def __lt__(self, other: '_print_key') -> bool:
        ...

def print_sorted(expressions: Iterable) -> list:
    ...

class _math_key:
    def __init__(self, ex: Any):
        ...

    def __lt__(self, other: '_math_key') -> bool:
        ...

def math_sorted(expressions: Iterable) -> list:
    ...

def mixed_order(lhs: Any, rhs: Any) -> int:
    ...

class _mixed_key:
    def __init__(self, ex: Any):
        ...

    def __lt__(self, other: '_mixed_key') -> bool:
        ...

def mixed_sorted(expressions: Iterable) -> list:
    ...
