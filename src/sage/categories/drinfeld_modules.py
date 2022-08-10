r"""
Drinfeld modules
"""
#*****************************************************************************
#  Copyright (C) 2022      Xavier Caruso <xavier.caruso@normalesup.org>
#                          Antoine Leudière <antoine.leudiere@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import CategoryWithParameters
from sage.misc.functional import log
from sage.rings.morphism import RingHomomorphism
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

from sage.categories.homsets import Homsets

# from sage.misc.cachefunc import cached_method
# from sage.categories.basic import Fields

class DrinfeldModules(CategoryWithParameters):

    def __init__(self, morphism, name='t'):
        gamma = morphism
        # Check input is a ring Morphism
        if not isinstance(gamma, RingHomomorphism):
            raise TypeError('input must be a Ring morphism')
        self._morphism = morphism
        self._domain = gamma.domain()
        # Check domain is Fq[X]
        function_ring = self._domain
        if not isinstance(function_ring, PolynomialRing_general):
            raise NotImplementedError('domain must be a polynomial ring')
        function_ring_base = function_ring.base_ring()
        if not function_ring_base.is_field() or not function_ring_base.is_finite() :
            raise TypeError('the base ring of the domain must be a finite field')
        Fq = function_ring_base
        FqX = function_ring
        X = FqX.gen()
        # Check domain is field
        K = gamma.codomain()
        if not K.is_field():
            raise TypeError('the codomain must be a field')
        # Build K{t}
        d = log(Fq.cardinality(), Fq.characteristic())
        tau = K.frobenius_endomorphism(d)
        self._codomain = OrePolynomialRing(K, tau, names=name)
        # Create characteristic
        self._characteristic = None
        if K.is_finite():
            f = gamma * FqX.coerce_map_from(Fq)  # Fq -> K
            E = K.over(f)
            self._characteristic = FqX(E(gamma(X)).minpoly())

    def base(self):
        return self.codomain().base_ring()

    def characteristic(self):
        if self._characteristic is None:
            raise NotImplementedError
        return self._characteristic 

    def _call_(self, gen):
        # Avoid circular import
        from sage.rings.function_field.drinfeld_module import DrinfeldModule
        # If gen is not in the codomain, an exception is raised
        gen = self._codomain(gen)
        if self.characteristic()(gen[0]) != 0:
            raise ValueError('incorrect characteristic')
        return DrinfeldModule(self._domain, gen)

    def super_categories(self):
        return []

    def _repr_(self):
        return f'Category of Drinfeld modules defined by {self._morphism}'

    def _make_named_class_key(self, name):
        return self._domain.category()

    def Homsets(self):
        return Homsets()

    def Endsets(self):
        return Homsets()

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def morphism(self):
        return self._morphism

    class ParentMethods:
        def characteristic(self):
            return self.category().characteristic()

    class ElementMethods:
        pass
