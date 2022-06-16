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
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

# from sage.misc.cachefunc import cached_method
# from sage.categories.basic import Fields

class DrinfeldModules(CategoryWithParameters):
    r"""
    The category of Drinfeld modules.

    EXAMPLES:

    We create the category of function fields::

        sage: C = FunctionFields()
        sage: C
        Category of function fields

    TESTS::

        sage: TestSuite(FunctionFields()).run()
    """

    def __init__(self, gamma, name='t'):
        r"""
        """
        self._gamma = gamma
        self._domain = FqX = gamma.domain()
        K = gamma.codomain()
        if not isinstance(FqX, PolynomialRing_general):
            raise NotImplementedError('domain must be a polynomial ring')
        Fq = FqX.base_ring()
        if not Fq.is_field() or not Fq.is_finite() :
            raise TypeError('the base ring of the domain must be a finite field')
        d = log(Fq.cardinality(), Fq.characteristic())
        tau = K.frobenius_endomorphism(d)
        self._codomain = OrePolynomialRing(K, tau, names=name)
        # Create characteristic
        self._characteristic = None
        if K.is_finite():
            f = gamma * FqX.coerce_map_from(Fq)
            E = K.over(f)
            self._characteristic = E(gamma(FqX.gen())).minpoly()

    def characteristic(self):
        if self._characteristic is None:
            raise NotImplementedError
        return self._characteristic 

    def _call_(self, gen):
        r"""
        Constructs an object in this category from the data in ``x``,
        or throws a TypeError.
        """
        from sage.rings.function_field.finite_drinfeld_module import FiniteDrinfeldModule
        gen = self._codomain(gen)
        if self.characteristic()(gen[0]) != 0:
            raise ValueError('incorrect characteristic')
        return FiniteDrinfeldModule(self._domain, gen)

    def super_categories(self):
        return []

    def _repr_(self):
        return f'Category of Drinfeld modules defined by {self._gamma}'

    def _make_named_class_key(self, name):
        return self._domain.category()

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def gamma(self):
        return self._gamma

    class ParentMethods:
        def characteristic(self):
            return self.category().characteristic()

    class ElementMethods:
        pass
