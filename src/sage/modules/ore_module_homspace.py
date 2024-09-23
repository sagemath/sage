from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.homset import HomsetWithBase
from sage.matrix.matrix_space import MatrixSpace

from sage.modules.ore_module import OreModule
from sage.modules.ore_module_morphism import OreModule_morphism


class OreModule_homspace(UniqueRepresentation, HomsetWithBase):
    Element = OreModule_morphism

    def __init__(self, domain, codomain, category=None):
        if not isinstance(domain, OreModule):
            raise ValueError("domain must be a Ore module")
        if not isinstance(codomain, OreModule):
            raise ValueError("codomain must be a Ore module")
        if domain.ore_ring() is not codomain.ore_ring():
            raise ValueError("domain and codomain must be defined over the same ring with same twisting maps")
        super().__init__(domain, codomain, category)
        base = domain.base_ring()
        self._matrix_space = MatrixSpace(base, domain.dimension(), codomain.dimension())

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def matrix_space(self):
        return self._matrix_space
