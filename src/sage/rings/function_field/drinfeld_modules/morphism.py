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
        # NOTE: it used to be x.parent() is parent, but this was False.
        # DrinfeldModuleHomset inherits Homset, which does NOT inherit
        # UniqueRepresentation
        if x.parent() == parent:  # x is a DrinfeldModuleMorphism
            ore_pol = x.ore_polynomial()
        else:  # x is an Ore polynomial
            ore_pol = domain.ore_polring()(x)
        if ore_pol * domain.gen() != codomain.gen() * ore_pol:
            raise ValueError('the Ore polynomial does not define a morphism')
        self._domain = domain
        self._codomain = codomain
        self._ore_polynomial = ore_pol

    # NOTE: Should I inherit UniqueRepresentation to avoid this?
    def __eq__(self, other):
        try:
            if self.parent() == other.parent():
                return self.ore_polynomial() == other.ore_polynomial()
        except AttributeError:
                return False
        return False

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

    def is_isomorphism(self):
        return self.is_isogeny() and self._ore_polynomial.degree() == 0
