from __future__ import annotations
from collections.abc import Callable
from typing import Any

from sage.symbolic.function import SymbolicFunction


def function_factory(
    name: str,
    nargs: int = 0,
    latex_name: str | None = None,
    conversions: dict[str, Any] | None = None,
    evalf_params_first: bool = True,
    eval_func: Callable | None = None,
    evalf_func: Callable | None = None,
    conjugate_func: Callable | None = None,
    real_part_func: Callable | None = None,
    imag_part_func: Callable | None = None,
    derivative_func: Callable | None = None,
    tderivative_func: Callable | None = None,
    power_func: Callable | None = None,
    series_func: Callable | None = None,
    print_func: Callable | None = None,
    print_latex_func: Callable | None = None,
) -> SymbolicFunction:
    ...


def unpickle_function(
    name: str,
    nargs: int,
    latex_name: str | None,
    conversions: dict[str, Any] | None,
    evalf_params_first: bool,
    pickled_funcs: list[Any],
) -> SymbolicFunction:
    ...


def function(
    s: str,
    **kwds: Any
) -> SymbolicFunction | list[SymbolicFunction]:
    ...
