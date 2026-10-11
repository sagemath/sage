r"""
Congruence subgroups `\Gamma_0(N)` of `\mathrm{GL}_{2}(\mathbb{F}_{q}[T])`.

By definition, it is the subgroup of matrices whose bottom left coefficient
is divislble by `N`.

AUTHORS:

- Cécile Armana, Xavier Caruso (2026-02): initial version
"""

# ****************************************************************************
#       Copyright (C) 2026 Cécile Armana <armana@math.cnrs.fr>
#                          Xavier Caruso <xavier@caruso.ovh>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.finite_fields import FiniteFields
from sage.categories.homset import Hom
from sage.categories.map import Map
from sage.categories.pushout import pushout
from sage.groups.group import Group
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation


class Gamma0Element(MultiplicativeGroupElement):
    r"""
    An element in a congruence subgroup.
    """
    def __init__(self, parent, elt):
        r"""
        Initialize this element of the congruence subgroup.

        TESTS:

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: g = G.an_element()
            sage: TestSuite(g).run()

        ::

            sage: G([1, 0, T^4 + 2*T + 3, 1])
            [            1             0]
            [T^4 + 2*T + 3             1]

        ::

            sage: G([1, T, T^4 + 2*T + 3, 1])
            Traceback (most recent call last):
            ...
            ValueError: not in Gamma0(T^4 + 2*T + 3)

        ::

            sage: G([1, 0, T, 1])
            Traceback (most recent call last):
            ...
            ValueError: not in Gamma0(T^4 + 2*T + 3)
        """
        MultiplicativeGroupElement.__init__(self, parent)
        MS = parent.matrix_space()
        level = parent.level()
        if isinstance(elt, Gamma0Element):
            mat = MS(elt._mat)
        else:
            mat = MS(elt)
        if not mat.is_invertible():
            raise ValueError("not in Gamma0(%s)" % level)
        if mat[1,0] % level != 0:
            raise ValueError("not in Gamma0(%s)" % level)
        self._mat = mat

    def __repr__(self):
        r"""
        Return the string representation of the element this element.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: g = G.an_element()
            sage: g.__repr__()
            '[4*T^4 + 3*T + 3 4*T^4 + 3*T + 2]\n[  T^4 + 2*T + 3   T^4 + 2*T + 4]'
        """
        return self._mat.__repr__()

    def _richcmp_(self, other, op):
        r"""
        Compare this element with ``other`` for the operator ``op``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: N = T^4 + 2*T + 3
            sage: G = Gamma0(N)
            sage: g = G.an_element()
            sage: h = G(matrix([[1,0],[N,1]]))
            sage: g == h
            False
        """
        return richcmp(self._mat, other._mat, op)

    def __reduce__(self):
        r"""
        Return data defining this element (for pickling).

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: g = G.an_element()
            sage: loads(dumps(g)) == g
            True
        """
        return Gamma0Element, (self.parent(), self._mat)

    def is_one(self):
        r"""
        Return whether this matrix is the identity matrix.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: g = G.an_element()
            sage: g
            [4*T^4 + 3*T + 3 4*T^4 + 3*T + 2]
            [  T^4 + 2*T + 3   T^4 + 2*T + 4]
            sage: g.is_one()
            False
            sage: h = g * g^(-1)
            sage: h.is_one()
            True
        """
        return self._mat.is_one()

    def _mul_(self, other):
        r"""
        Return the product of this element with ``other``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]; N = T^4+2*T+3
            sage: G = Gamma0(N)
            sage: g = G.an_element()
            sage: h = G(matrix([[1,0],[N,1]]))
            sage: g._mul_(h)
            [4*T^8 + T^5 + 3*T^4 + T^2 + T + 4                   4*T^4 + 3*T + 2]
            [  T^8 + 4*T^5 + 3*T^4 + 4*T^2 + T                     T^4 + 2*T + 4]
        """
        return self.parent()(self._mat * other._mat)

    def __invert__(self):
        r"""
        Return the inverse of this element.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4+2*T+3)
            sage: g = G.an_element()
            sage: g.__invert__()
            [  T^4 + 2*T + 4   T^4 + 2*T + 3]
            [4*T^4 + 3*T + 2 4*T^4 + 3*T + 3]
        """
        return self.parent()(self._mat.inverse())

    def level(self):
        r"""
        Return the level of the congruence subgroup in which this
        element lives.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: g = G.an_element()
            sage: g.level()
            T^4 + 2*T + 3
        """
        return self.parent().level()


class Gamma0_drinfeld(Group, UniqueRepresentation):
    r"""
    A congruence subgroup.
    """
    Element = Gamma0Element

    @staticmethod
    def __classcall__(cls, level):
        r"""
        Normalize and check the input.

        INPUT:

        - ``level`` -- a polynomial

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T + 1)
            sage: H = Gamma0(2*T + 2)
            sage: G is H
            True

        ::

            sage: B.<T> = GF(5^2)[]
            sage: G2 = Gamma0(T + 1)
            sage: G is G2
            False

        ::

            sage: Gamma0(A(0))
            Traceback (most recent call last):
            ...
            ValueError: level must be nonzero
        """
        base = level.parent()
        if not (isinstance(level, Polynomial) and base.base_ring() in FiniteFields()):
            raise TypeError("base ring must be a polynomial ring over a finite field")
        if level == 0:
            raise ValueError("level must be nonzero")
        level = level.monic()
        return super().__classcall__(cls, base, level)

    def __init__(self, base, level):
        r"""
        The congruence subgroup `\Gamma_0(N)`.

        INPUT:

        - ``base`` -- a polynomial ring

        - ``level`` -- a monic polynomial over ``base``

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(A(1))
            Congruence Subgroup Gamma0(1)
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: G
            Congruence Subgroup Gamma0(T^4 + 2*T + 3)
            sage: TestSuite(G).run()
        """
        Group.__init__(self)
        self._base = base
        self._level = level
        self._matrix_space = MS = MatrixSpace(base, 2)
        self._q = base.base_ring().cardinality()
        coerce = InclusionIntoMatrixSpace(Hom(self, MS))
        MS.register_coercion(coerce)

    def __reduce__(self):
        r"""
        Return data defining this parent (for pickling).

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 + 2*T + 3)
            sage: g = G.an_element()
            sage: loads(dumps(G)) is G
            True
        """
        return Gamma0_drinfeld, (self._level,)

    def _repr_(self):
        r"""
        Return the string representation of this congruence subgroup.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 + 2*T + 3)._repr_()
            'Congruence Subgroup Gamma0(T^4 + 2*T + 3)'
        """
        return "Congruence Subgroup Gamma0(%s)" % self._level

    def _coerce_map_from_(self, other):
        r"""
        Return a coerce map from ``other`` to `self``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^2 + 3*T + 2)
            sage: H = Gamma0(T + 1)
            sage: G.has_coerce_map_from(H)  # indirect doctest
            False
            sage: H.has_coerce_map_from(G)  # indirect doctest
            True
        """
        if isinstance(other, Gamma0_drinfeld):
            return other.level() % self.level() == 0

    def _pushout_(self, other):
        r"""
        Return a coerce map from ``other`` to `self``.

        EXAMPLES::

            sage: from sage.categories.pushout import pushout
            sage: A.<T> = GF(5^2)[]
            sage: G = Gamma0(T^2 + 3*T + 2)
            sage: B.<T> = GF(5^3)[]
            sage: H = Gamma0(T^2 + 4*T + 4)
            sage: P = pushout(G, H)
            sage: P
            Congruence Subgroup Gamma0(T + 2)
            sage: P.base_ring()
            Univariate Polynomial Ring in T over Finite Field in z6 of size 5^6
        """
        if isinstance(other, Gamma0_drinfeld):
            base = pushout(self.base_ring(), other.base_ring())
            if base is not None:
                N = base(self.level()).gcd(base(other.level()))
                return Gamma0_drinfeld(N)

    def base_ring(self):
        r"""
        Return the base ring of this congruence subgroup.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 + 2*T + 3).base_ring()
            Univariate Polynomial Ring in T over Finite Field of size 5
        """
        return self._base

    def level(self):
        r"""
        Return the level of this congruence subgroup.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 + 2*T + 3).level()
            T^4 + 2*T + 3
        """
        return self._level

    @lazy_attribute
    def _level_factorized(self):
        r"""
        Return the factorization of the level of this congruence subgroup.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 + 2*T + 3)._level_factorized
            (T + 2)^2 * (T^2 + T + 2)
        """
        return self._level.factor()

    def _an_element_(self):
        r"""
        Return an element of this congruence subgroup.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 + 2*T + 3).an_element()
            [4*T^4 + 3*T + 3 4*T^4 + 3*T + 2]
            [  T^4 + 2*T + 3   T^4 + 2*T + 4]
        """
        N = self.level()
        return self([1-N, -N, N, 1+N])

    def matrix_space(self):
        r"""
        Return the matrix space of this congruence subgroup.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 + 2*T + 3).matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial
            Ring in T over Finite Field of size 5
        """
        return self._matrix_space

    def genus(self):
        r""""
        Return the genus of this congruence subgroup, i.e. the genus of
        the attached Drinfeld modular curve.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T - 1).genus()
            0
            sage: Gamma0(T^3 + 1).genus()
            5
            sage: Gamma0(T^5 - 3*T^4 + 1).genus()
            155

        If the level `N` is irreducible, the formula simplifies as follows::

            sage: def genus_irr(N):
            ....:     q = N.base_ring().cardinality()
            ....:     d = N.degree()
            ....:     if is_even(d):
            ....:          g = (q^d - q^2) / (q^2 - 1)
            ....:     else:
            ....:          g = (q^d - q) / (q^2 - 1)
            ....:     return g

            sage: N = T^6 + T^4 + 4*T^3 + T^2 + 2
            sage: N.is_irreducible()
            True
            sage: Gamma0(N).genus() == genus_irr(N)
            True

            sage: N = T^7 + 3*T + 3
            sage: N.is_irreducible()
            True
            sage: Gamma0(N).genus() == genus_irr(N)
            True

        TESTS::

            sage: Gamma0(A(1)).genus()
            0

        REFERENCE:

        [Gek2001]_, Theorem 8.1
        """
        q = self._q
        L = self._level_factorized
        s = len(L)
        epsilon = kappa = r = 1
        for P, mult in L:
            d = P.degree()
            epsilon *= q**(d*(mult-1)) * (1 + q**d)
            kappa *= q**(d * (mult // 2)) + q**(d * ((mult-1) // 2))
            if d % 2 == 1:
                r = 0
        g = 1 + (epsilon - (q+1)*kappa - 2**(s-1) * (r*q*(q-1) + (q+1)*(q-2))) // (q**2 - 1)
        return ZZ(g)

    def ncusps(self):
        r"""
        Return the number of cusps of this congruence subgroup, i.e. the
        number of cusps of the attached Drinfeld modular curve.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T - 1).ncusps()
            2
            sage: Gamma0(T^3 + 1).ncusps()
            4
            sage: Gamma0(T^7 - 3*T^4 + 1).ncusps()
            8

        If the level `N` is irreducible, the number of cusps is `2`.
        We check it below::

            sage: N = T^6 + T^4 + 4*T^3 + T^2 + 2
            sage: N.is_irreducible()
            True
            sage: Gamma0(N).ncusps() == 2
            True

        ::

            sage: N = T^7 + 3*T + 3
            sage: N.is_irreducible()
            True
            sage: Gamma0(N).ncusps() == 2
            True

        REFERENCE:

        [Gek2001]_, Proposition 6.7
        """
        q = self._q
        L = self._level_factorized
        s = len(L)
        kappa = 1
        for P, mult in L:
            d = P.degree()
            kappa *= q**(d * (mult // 2)) + q**(d * ((mult-1) // 2))
        h = 2**s + (kappa - 2**s) // (q - 1)
        return ZZ(h)

    def index(self):
        r"""
        Return the index of this congruence subgroup as a subgroup of
        `\mathrm{GL}_{2}(A)` (where `A` is the base ring).

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: Gamma0(T^4 - 3*T^2 + 1).index()
            900

        If the level `N` is irreducible, the index is `1 + q^{\mathrm{deg}(N)}`.
        We check it on an example below::

            sage: N = T^4 + 4*T^2 + 4*T + 2
            sage: N.is_irreducible()
            True
            sage: q = A.base_ring().cardinality()
            sage: Gamma0(N).index() == 1 + q^4
            True
        """
        q = self._q
        ind = q ** self._level.degree()
        ind *= prod(1 + q**(-P.degree()) for P, _ in self._level_factorized)
        return ZZ(ind)


class InclusionIntoMatrixSpace(Map):
    r"""
    Inclusion of a congruence subgroup into the corresponding matrix space.
    """
    def __init__(self, parent):
        r"""
        Initialize this morphism.

        TESTS::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 - 3*T^2 + 1)
            sage: M = G.matrix_space()
            sage: f = M.coerce_map_from(G)
            sage: f
            Generic map:
              From: Congruence Subgroup Gamma0(T^4 + 2*T^2 + 1)
              To:   Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in T over Finite Field of size 5
            sage: type(f)
            <class 'sage.modular.drinfeld_modform.congroup_gamma0.InclusionIntoMatrixSpace'>
            sage: TestSuite(f).run(skip='_test_category')
        """
        Map.__init__(self, parent)

    def _call_(self, g):
        r"""
        Return `g` viewed as an element of a matrix space.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 - 3*T^2 + 1)
            sage: M = G.matrix_space()
            sage: f = M.coerce_map_from(G)
            sage: g = G.an_element()
            sage: h = f(g)
            sage: h.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in T over Finite Field of size 5
        """
        return g._mat

    def _richcmp_(self, other, op):
        r"""
        Compare this map with ``other`` for the operator ``op``.

        EXAMPLES::

            sage: A.<T> = GF(5)[]
            sage: G = Gamma0(T^4 - 3*T^2 + 1)
            sage: M = G.matrix_space()
            sage: g1 = M.coerce_map_from(G)
            sage: g2 = M.coerce_map_from(G)
            sage: g1 == g2
            True

        ::

            sage: H = Gamma0(T^2 - 3*T + 1)
            sage: h = M.coerce_map_from(H)
            sage: g1 == h
            False
        """
        return richcmp(self.parent(), other.parent(), op)
