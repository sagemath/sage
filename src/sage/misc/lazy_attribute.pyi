# This type-stub file helps pyright understand the decorator @lazy_attribute.

from collections.abc import Callable
from typing import Any

# Adapted from https://github.com/python/typeshed/blob/b9640005eb586afdbe0a57bac2b88a7a12465069/stdlib/builtins.pyi#L1237-L1254
class lazy_attribute:
    def __init__(
        self,
        f: Callable[[Any], Any] | None = ...
    ) -> None: ...
    def __get__(self, a: Any, cls: type) -> Any: ...
