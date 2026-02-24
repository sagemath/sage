from collections.abc import Callable
from typing import Any

def dict_key(o: Any) -> Any:
    ...

def cache_key(o: Any) -> Any:
    ...

def cached_method(f, name: str | None = None, key=None, do_pickle: bool = False) -> CachedMethod:
    ...

def cached_function(f, name: str | None = None, key=None, do_pickle: bool = False) -> CachedFunction:
    ...

class CachedFunction:
    def __init__(self, f: Callable, classmethod: bool = False,
                 name: str | None = None, key: Callable | None = None,
                 do_pickle: bool = False) -> None:
        ...

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        ...

    def cached(self, *args: Any, **kwds: Any) -> Any:
        ...

    def is_in_cache(self, *args: Any, **kwds: Any) -> bool:
        ...

    def set_cache(self, value: Any, *args: Any, **kwds: Any) -> None:
        ...

    def get_key(self, *args: Any, **kwds: Any) -> Any:
        ...

    def __repr__(self) -> str:
        ...

    def clear_cache(self) -> None:
        ...

    def precompute(self, arglist: Any, num_processes: int = 1) -> None:
        ...

class CachedMethod:
    def __init__(self, f: Callable, name: str | None = None,
                 key: Callable | None = None,
                 do_pickle: bool = False) -> None:
        ...

    def __call__(self, inst: Any, *args: Any, **kwds: Any) -> Any:
        ...

    def _get_instance_cache(self, inst: Any) -> dict:
        ...

    def __get__(self, inst: Any, cls: Any) -> Any:
        ...

class CacheDict(dict):
    pass

class CachedInParentMethod(CachedMethod):
    def __init__(self, f: Callable, name: str | None = None,
                 key: Callable | None = None,
                 do_pickle: bool = False) -> None:
        ...

    def _get_instance_cache(self, inst: Any) -> dict:
        ...

    def __get__(self, inst: Any, cls: Any) -> Any:
        ...

class CachedMethodCaller(CachedFunction):
    def __init__(self, cachedmethod: CachedMethod, inst: Any,
                 cache: dict | None = None, name: str | None = None,
                 key: Callable | None = None,
                 do_pickle: bool = False) -> None:
        ...

    def _instance_call(self, *args: Any, **kwds: Any) -> Any:
        ...

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        ...

    def cached(self, *args: Any, **kwds: Any) -> Any:
        ...

    def __get__(self, inst: Any, cls: Any) -> Any:
        ...

    def precompute(self, arglist: Any, num_processes: int = 1) -> None:
        ...

class CachedMethodCallerNoArgs(CachedFunction):
    def __init__(self, inst: Any, f: Callable, cache: Any = None,
                 name: str | None = None,
                 do_pickle: bool = False) -> None:
        ...

    def _instance_call(self) -> Any:
        ...

    def __call__(self) -> Any:
        ...

    def set_cache(self, value: Any) -> None:
        ...

    def clear_cache(self) -> None:
        ...

    def is_in_cache(self) -> bool:
        ...

    def __get__(self, inst: Any, cls: Any) -> Any:
        ...

class GloballyCachedMethodCaller(CachedMethodCaller):
    def get_key_args_kwds(self, args: tuple, kwds: dict) -> Any:
        ...
