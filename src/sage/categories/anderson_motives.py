r"""
Anderson motives

AUTHOR:

- Xavier Caruso, Antoine Leudière (2025-11): initial version
"""

# *****************************************************************************
#        Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************


from sage.misc.latex import latex

from sage.categories.modules import Modules
from sage.categories.ore_modules import OreModules
from sage.categories.homsets import Homsets
from sage.categories.drinfeld_modules import DrinfeldModules

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing


class AndersonMotives(OreModules):
    r"""
    The category of Anderson motives.

    .. SEEALSO::

        :class:`sage.category.drinfeld_modules.DrinfeldModules`,
        :mod:`sage.rings.function_field.drinfeld_modules.anderson_motive`
    """
    @staticmethod
    def __classcall_private__(cls, category):
        r"""
        Normalize the input and construct the category.

        INPUT:

        - ``category`` -- the corresponding category of Drinfeld
          modules, or an object for constructing it

        TESTS::

            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: C1 = AndersonMotives(phi.category())
            sage: C1
            Category of Anderson motives over Univariate Polynomial Ring in T over Finite Field in z of size 3^3

            sage: C2 = AndersonMotives(phi.base_morphism())
            sage: C2
            Category of Anderson motives over Univariate Polynomial Ring in T over Finite Field in z of size 3^3

            sage: C1 is C2
            True
        """
        if isinstance(category, AndersonMotives):
            return category
        if isinstance(category, DrinfeldModules):
            category = DrinfeldModules(category.base_morphism())
        else:
            category = DrinfeldModules(category)
        return AndersonMotives.__classcall__(cls, category)

    def __init__(self, category):
        r"""
        Initialize this category.

        INPUT:

        - ``category`` -- the corresponding category of Drinfeld modules

        TESTS::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: TestSuite(C).run()
        """
        self._drinfeld_category = category
        self._base_morphism = category.base_morphism()
        self._A_field = category.base()
        self._function_ring = A = category.function_ring()
        self._base_over_constants_field = category.base_over_constants_field()
        self._ore_variable_name = category._ore_variable_name
        self._characteristic = category._characteristic
        K = self._base_morphism.codomain()
        self._base_combined = AK = PolynomialRing(K, A.variable_name())  # TODO: find a better name
        self._constant_coefficient = category.constant_coefficient()
        self._divisor = AK.gen() - self._constant_coefficient
        twisting_morphism = category.ore_polring().twisting_morphism()
        twisting_morphism = AK.hom([AK.gen()], base_map=twisting_morphism)
        self._ore_polring = OrePolynomialRing(AK, twisting_morphism, names=self._ore_variable_name)
        super().__init__(self._ore_polring)

    def _repr_(self):
        r"""
        Return a string representation of this category.

        TESTS::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C  # indirect doctest
            Category of Anderson motives over Univariate Polynomial Ring in T over Finite Field in z of size 3^3
        """
        return f'Category of Anderson motives over {self.base()}'

    def _latex_(self):
        r"""
        Return a LaTeX representation of this category.

        TESTS::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: latex(C)  # indirect doctest
            \text{Category{ }of{ }Anderson{ }motives{ }over{ }\Bold{F}_{3^{3}}[T]
        """
        return f'\\text{{Category{{ }}of{{ }}Anderson{{ }}motives{{ }}' \
               f'over{{ }}{latex(self.base())}'

    def __reduce__(self):
        r"""
        Return the necessary arguments to construct this object,
        as per the pickle protocol.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: loads(dumps(C)) is C
            True
        """
        return AndersonMotives, (self._drinfeld_category,)

    def Homsets(self):
        r"""
        Return the category of homsets.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())

            sage: from sage.categories.homsets import Homsets
            sage: C.Homsets() is Homsets()
            True
        """
        return Homsets()

    def Endsets(self):
        r"""
        Return the category of endsets.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())

            sage: from sage.categories.homsets import Homsets
            sage: C.Endsets() is Homsets().Endsets()
            True
        """
        return Homsets().Endsets()

    def A_field(self):
        r"""
        Return the underlying `A`-field of this category.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C.A_field()
            Finite Field in z of size 3^3 over its base
        """
        return self._A_field

    def base(self):
        r"""
        Return the base over which the Anderson motives in this
        category are defined.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C.base()
            Univariate Polynomial Ring in T over Finite Field in z of size 3^3
        """
        return self._base_combined

    def divisor(self):
        r"""
        Return the polynomial `T - z` if `T` denotes the generator of
        the function ring `A` and `z` is the image of `T` in the `A`-field.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C.divisor()
            T + 2*z
        """
        return self._divisor

    def characteristic(self):
        r"""
        Return the characteristic of the underlying `A`-field.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C.divisor()
            T + 2*z
        """
        if self._characteristic is None:
            raise NotImplementedError('function ring characteristic not '
                                      'implemented in this case')
        return self._characteristic

    def function_ring(self):
        r"""
        Return the underlying function ring of this category.

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C.function_ring()
            Univariate Polynomial Ring in T over Finite Field of size 3
        """
        return self._function_ring

    def object(self, tau=None):
        r"""
        Return the object in this category with `\tau`-action
        given by the matrix ``tau``.

        INPUT:

        - ``tau`` -- a matrix or ``None`` (default: ``None``);
          if ``None``, return the trivial Anderson module in this
          category

        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())

            sage: M = C.object()
            sage: M
            Anderson motive of rank 1 over Univariate Polynomial Ring in T over Finite Field in z of size 3^3
            sage: M.matrix()
            [1]

        ::

            sage: tau = matrix(2, 2, [[T, 1], [z, 1]])
            sage: N = C.object(tau)
            sage: N.matrix()
            [T 1]
            [z 1]
        """
        from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive
        return AndersonMotive(self, tau)

    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.anderson_motives import AndersonMotives
            sage: A.<T> = GF(3)[]
            sage: K.<z> = GF(3^3)
            sage: phi = DrinfeldModule(A, [z, z^2, z^3])
            sage: C = AndersonMotives(phi.category())
            sage: C.super_categories()
            [Category of Ore modules over Univariate Polynomial Ring in T over Finite Field in z of size 3^3 twisted by T |--> T, with map of base ring]
        """
        AKtau = self._ore_polring
        return [OreModules(AKtau.base(), AKtau)]

    class ParentMethods:

        def function_ring(self):
            r"""
            Return the underlying function ring of this Anderson motive.

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z> = Fq.extension(6)
                sage: phi = DrinfeldModule(A, [z, z^3, z^5])
                sage: M = phi.anderson_motive()
                sage: M.function_ring()
                Univariate Polynomial Ring in T over Finite Field in z2 of size 5^2
            """
            return self._category.function_ring()

        def A_field(self):
            r"""
            Return the underlying `A`-field of this Anderson motive.

            This is an instance of the class
            :class:`sage.rings.ring_extension.RingExtension`.

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z> = Fq.extension(6)
                sage: phi = DrinfeldModule(A, [z, z^3, z^5])
                sage: M = phi.anderson_motive()
                sage: M.A_field()
                Finite Field in z of size 5^12 over its base
            """
            return self._category.A_field()

        def base(self):
            r"""
            Return the base ring over which this Anderson motive
            is defined.

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z> = Fq.extension(6)
                sage: phi = DrinfeldModule(A, [z, z^3, z^5])
                sage: M = phi.anderson_motive()
                sage: M.base()
                Univariate Polynomial Ring in T over Finite Field in z of size 5^12
            """
            return self._category.base()

        def characteristic(self):
            r"""
            Return the characteristic of the underlying `A`-field.

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z> = Fq.extension(6)
                sage: phi = DrinfeldModule(A, [z, z^3, z^5])
                sage: M = phi.anderson_motive()
                sage: M.characteristic()
                T^6 + (4*z2 + 3)*T^5 + T^4 + (3*z2 + 1)*T^3 + T^2 + (4*z2 + 1)*T + z2
            """
            return self._category.characteristic()

        def ore_polring(self):
            r"""
            Return the Ore polynomial ring over which this Anderson
            motive is defined.

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z> = Fq.extension(6)
                sage: phi = DrinfeldModule(A, [z, z^3, z^5])
                sage: M = phi.anderson_motive()
                sage: M.ore_polring()
                Ore Polynomial Ring in τ over Univariate Polynomial Ring in T
                over Finite Field in z of size 5^12 twisted by T |--> T, with map of base ring
            """
            return self._category._ore_polring

        def ore_variable(self):
            r"""
            Return the generator of the Ore polynomial ring over which
            this Anderson motive is defined.

            EXAMPLES::

                sage: Fq = GF(25)
                sage: A.<T> = Fq[]
                sage: K.<z> = Fq.extension(6)
                sage: phi = DrinfeldModule(A, [z, z^3, z^5])
                sage: M = phi.anderson_motive()
                sage: tau = M.ore_variable()
                sage: tau
                τ
                sage: tau.parent()
                Ore Polynomial Ring in τ over Univariate Polynomial Ring in T
                over Finite Field in z of size 5^12 twisted by T |--> T, with map of base ring
            """
            return self._category._ore_polring.gen()
