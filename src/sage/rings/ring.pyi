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
    _default_category: Category = _CommutativeRings

    def __init__(
        self,
        base_ring: Parent[Any] | object,
        names: NameSpec = None,
        normalize: bool = True,
        category: Category | None = None,
    ) -> None: ...
    def fraction_field(self) -> Ring: ...
    def extension(
        self,
        poly: object,
        name: str | None = None,
        names: NameSpec = None,
        **kwds: object,
    ) -> Ring: ...

class IntegralDomain(CommutativeRing):
    _default_category: Category
    def __init__(self, *args, **kwds) -> None: ...

class NoetherianRing(CommutativeRing):
    _default_category: Category
    def __init__(self, *args, **kwds) -> None: ...

class DedekindDomain(CommutativeRing):
    _default_category: Category
    def __init__(self, *args, **kwds) -> None: ...

class PrincipalIdealDomain(CommutativeRing):
    _default_category: Category
    def __init__(self, *args, **kwds) -> None: ...

def _is_Field(x: object) -> bool: ...

class Field(CommutativeRing):
    _default_category: Category = _Fields

class Algebra(Ring):
    def __init__(
        self, base_ring: Parent[Any] | object, *args: object, **kwds: object
    ) -> None: ...

class CommutativeAlgebra(CommutativeRing):
    def __init__(
        self, base_ring: Parent[Any] | object, *args: object, **kwds: object
    ) -> None: ...

def is_Ring(x: object) -> bool: ...
