#*****************************************************************************
#       Copyright (C) 2022 Antoine Leudière <antoine.leudiere@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# from sage.categories.morphism import Morphism
from sage.structure.element import Element
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.categories.drinfeld_modules import DrinfeldModules


class DrinfeldModuleMorphism(Element):

    def __init__(self, parent, x):
        super().__init__(parent)
        domain = parent.domain()
        codomain = parent.codomain()
        if x.parent() is parent:
            ore_polynomial = x.ore_polynomial()
        else:
            ore_polynomial = domain.ore_polring()(x)
        # Test well-definition of the morphism
        if domain.gen() * ore_polynomial != ore_polynomial * codomain.gen():
            raise ValueError('the Ore polynomial does not define a morphism')
        # Instantiation
        self._domain = domain
        self._codomain = codomain
        self._ore_polynomial = ore_polynomial

    def _repr_(self):
        return f'Drinfeld Module morphism:\n' \
                f'  From: {self._domain}\n'  \
                f'  To:   {self._codomain}\n' \
                f'  Defn: {self._ore_polynomial}'

    def codomain(self):
        return self._codomain

    def domain(self):
        return self._domain

    def ore_polynomial(self):
        return self._ore_polynomial

    def is_zero(self):
        return self._ore_polynomial.is_zero()

    def is_isogeny(self):
        return not self.is_zero()
