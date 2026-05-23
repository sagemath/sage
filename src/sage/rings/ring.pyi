from __future__ import annotations

from collections.abc import Iterator
from typing import Any, NoReturn

from sage.categories.category import Category
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.categories.rings import Rings
from sage.structure.category_object import NameSpec
from sage.structure.parent import Parent
from sage.structure.parent_gens import ParentWithGens

_Rings: Rings
_CommutativeRings: CommutativeRings
_Fields: Fields

class Ring(ParentWithGens):
    def __init__(
        self,
        base: Parent[Any] | object,
        names: NameSpec = None,
        normalize: bool = True,
        category: Category | None = None,
    ) -> None: ...
    def __iter__(self) -> Iterator[Any]: ...
    def __len__(self) -> int: ...
    def __xor__(self, n: object) -> NoReturn: ...
    def base_extend(self, X: Parent[Any]) -> Parent[Any]: ...
    def category(self) -> Category: ...
    def __mul__(self, x: object) -> object: ...
    def zero(self) -> Any: ...
    def one(self) -> Any: ...
    def order(self) -> int: ...

class CommutativeRing(Ring):
    def __init__(self, *args, **kwds) -> None: ...

class IntegralDomain(Ring):
    def __init__(self, *args, **kwds) -> None: ...

class NoetherianRing(Ring):
    def __init__(self, *args, **kwds) -> None: ...

class DedekindDomain(Ring):
    def __init__(self, *args, **kwds) -> None: ...

class PrincipalIdealDomain(Ring):
    def __init__(self, *args, **kwds) -> None: ...

def _is_Field(x: object) -> bool: ...

class Field(Ring):
    pass

class Algebra(Ring):
    def __init__(
        self, base_ring: Parent[Any] | object, *args: object, **kwds: object
    ) -> None: ...

class CommutativeAlgebra(Ring):
    def __init__(
        self, base_ring: Parent[Any] | object, *args: object, **kwds: object
    ) -> None: ...
