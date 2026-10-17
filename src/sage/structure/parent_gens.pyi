from typing import Any

from sage.structure.category_object import CategoryObject

class ParentWithGens(CategoryObject):
    def gens(self) -> tuple[Any, ...]: ...

    def ngens(self) -> int: ...
