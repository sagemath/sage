r"""
Morphisms between Ore modules

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

from sage.misc.latex import latex
from sage.matrix.constructor import matrix
from sage.categories.morphism import Morphism
from sage.modules.ore_module import OreModule, OreSubmodule, OreQuotientModule

class OreModule_morphism(Morphism):
    def __init__(self, parent, matrix):
        Morphism.__init__(self, parent)
        self._matrix = parent.matrix_space()(matrix)
        for x in parent.domain().basis():
            if self._call_(x.image()) != self._call_(x).image():
                raise ValueError("does not define a morphism of Ore modules")

    def _repr_type(self):
        return "Ore module"

    def _latex_(self):
        s = "\\begin{array}{l}\n"
        s += "\\text{\\texttt{%s morphism:}} \\\\\n" % self._repr_type()
        s += "\\text{\\texttt{{ }{ }From:}}\\hspace{1ex} %s \\\\\n" % latex(self.domain())
        s += "\\text{\\texttt{{ }{ }To:}}\\hspace{3ex} %s \n" % latex(self.codomain())
        s += "\\end{array}"
        return s

    def matrix(self):
        return self._matrix.__copy__()

    def is_zero(self):
        return self._matrix.is_zero()

    def is_identity(self):
        return self.domain() is self.codmain() and self._matrix.is_one()

    def determinant(self):
        if self.domain() is not self.codomain():
            raise ValueError("determinants are only defined for endomorphisms")
        return self._matrix.determinant()

    det = determinant

    def characteristic_polynomial(self, var='x'):
        if self.domain() is not self.codomain():
            raise ValueError("characteristic polynomials are only defined for endomorphisms")
        return self._matrix.charpoly(var)

    charpoly = characteristic_polynomial

    def is_injective(self):
        return self._matrix.rank() == self.domain().rank()

    def is_surjective(self):
        return self._matrix.rank() == self.codomain().rank()

    def is_bijective(self):
        return self.is_injective() and self.is_surjective()

    def is_isomorphism(self):
        return self.is_bijective()

    def _call_(self, x):
        return self.codomain()(x * self._matrix)

    def _composition_(self, other, homset):
        if not isinstance(other, OreModule_morphism):
            raise ValueError
        return homset(other._matrix * self._matrix)

    def kernel(self, names=None):
        ker = self._matrix.left_kernel_matrix()
        return OreSubmodule(self.domain(), ker, names)

    def image(self, names=None):
        return OreSubmodule(self.codomain(), self._matrix, names)

    def cokernel(self, names=None):
        return OreQuotientModule(self.codomain(), self._matrix, names)

    def coimage(self, names=None):
        ker = self._matrix.left_kernel_matrix()
        return OreQuotientModule(self.domain(), ker, names)
