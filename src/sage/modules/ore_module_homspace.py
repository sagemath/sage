r"""
Space of morphisms between Ore modules

AUTHOR:

- Xavier Caruso (2024-10)
"""

# ***************************************************************************
#    Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.homset import HomsetWithBase
from sage.matrix.matrix_space import MatrixSpace

from sage.modules.ore_module import OreModule
from sage.modules.ore_module_morphism import OreModuleMorphism

class OreModule_homspace(UniqueRepresentation, HomsetWithBase):
    Element = OreModuleMorphism

    def __init__(self, domain, codomain, category=None):
        if not isinstance(domain, OreModule):
            raise ValueError("domain must be a Ore module")
        if not isinstance(codomain, OreModule):
            raise ValueError("codomain must be a Ore module")
        if domain.ore_ring(action=False) is not codomain.ore_ring(action=False):
            raise ValueError("domain and codomain must be defined over the same ring with same twisting maps")
        super().__init__(domain, codomain, category)
        base = domain.base_ring()
        self._matrix_space = MatrixSpace(base, domain.dimension(), codomain.dimension())

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def matrix_space(self):
        return self._matrix_space
