r"""
Skew Univariate Polynomial Rings

This module provides the
:class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general``,
which constructs a general dense univariate skew polynomials over commutative
base rings with automorphisms over the base rings. This is usual accessed only
indirectly through the constructor
:func:`sage.rings.polynomial.skew_polynomial_constructor.SkewPolynomialRing`.

See :class:`SkewPolynomialRing_general` for a definition of a univariate skew
polynomial ring.

AUTHOR:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests
  and refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and
  doctest errors

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from __future__ import print_function, absolute_import, division

from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.rings.ring import Algebra
from sage.categories.rings import Rings
from sage.rings.integer import Integer
from sage.structure.category_object import normalize_names
from sage.categories.morphism import Morphism
from sage.categories.morphism import IdentityMorphism
from sage.rings.polynomial.skew_polynomial_element import (SkewPolynomial,
                                                           SkewPolynomialBaseringInjection)

#########################################################################################

class SkewPolynomialRing_general(Algebra, UniqueRepresentation):
    r"""
    A general implementation of univariate skew polynomial ring over a commutative ring.

    Let `R` be a commutative ring, and let `\sigma` be an automorphism of
    `R`.  The ring of skew polynomials `R[X, \sigma]` is the polynomial
    ring `R[X]`, where the addition is the usual polynomial addition, but
    the multiplication operation is defined by the modified rule

    .. MATH::

        X*a = \sigma(a) X.

    This means that `R[X, \sigma]` is a non-commutative ring. Skew polynomials
    were first introduced by Ore [Ore33].

    EXAMPLES::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = SkewPolynomialRing(R,sigma); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    One can also use a shorter syntax::

        sage: S.<x> = R['x',sigma]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    If we omit the diamond notation, the variable holding the indeterminate is
    not assigned::

        sage: Sy = R['y',sigma]
        sage: y
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined
        sage: Sy.gen()
        y

    Note however that contrary to usual polynomial rings, we cannot omit the
    variable name on the RHS, since this collides with the notation for creating polynomial rings::

        sage: Sz.<z> = R[sigma]
        Traceback (most recent call last):
        ...
        ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n
        Defn: t |--> t + 1' is not alphanumeric

    As for polynomials, skew polynomial rings with different variable names
    are not equal::

        sage: R['x',sigma] == R['y',sigma]
        False

    Of course, skew polynomial rings with different twist maps are not
    equal either::

        sage: R['x',sigma] == R['x',sigma^2]
        False

    Saving and loading of polynomial rings works::

        sage: loads(dumps(R['x',sigma])) == R['x',sigma]
        True

    There is a coercion map from the base ring of the skew polynomial rings::

        sage: S.has_coerce_map_from(R)
        True
        sage: x.parent()
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1
        sage: t.parent()
        Univariate Polynomial Ring in t over Integer Ring
        sage: y = x+t; y
        x + t
        sage: y.parent() is S
        True

    .. SEEALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :mod:`sage.rings.polynomial.skew_polynomial_element`

    REFERENCES:

    .. [Ore33] Oystein Ore.
       *Theory of Non-Commutative Polynomials*
       Annals of Mathematics, Second Series, Volume 34,
       Issue 3 (Jul., 1933), 480-508.
    """
    @staticmethod
    def __classcall__(cls, base_ring, twist_map=None, name=None, sparse=False,
                      element_class=None):
        r"""
        Set the default values for ``name``, ``sparse`` and ``element_class``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R, sigma)
            sage: S.__class__(R, sigma, x)
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
             twisted by t |--> t + 1
        """
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_element
                element_class = skew_polynomial_element.SkewPolynomial_generic_dense
        if twist_map is None:
            twist_map = IdentityMorphism(base_ring)
        else:
            if not isinstance(twist_map, Morphism):
                raise TypeError("given map is not a ring homomorphism")
            if twist_map.domain() != base_ring or twist_map.codomain() != base_ring:
                raise TypeError("given map is not an automorphism of %s" % base_ring)
        return super(SkewPolynomialRing_general,cls).__classcall__(cls,
                         base_ring, twist_map, name, sparse, element_class)

    def __init__(self, base_ring, twist_map, name, sparse, element_class):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``twist_map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to
          be used in ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S([1]) + S([-1])
            0
            sage: TestSuite(S).run()

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x', Frob]; T
            Skew Polynomial Ring in x over Finite Field in t of size 5^3
             twisted by t |--> t^5

        We skip the pickling tests currently because ``Frob`` does not
        pickle correctly (see note on :trac:`13215`)::

            sage: TestSuite(T).run(skip=["_test_pickling", "_test_elements"])
        """
        self.__is_sparse = sparse
        self._polynomial_class = element_class
        self._map = twist_map
        self._maps = {0: IdentityMorphism(base_ring), 1: self._map}
        self._no_generic_basering_coercion = True
        Algebra.__init__(self, base_ring, names=name,
                         normalize=True, category=Rings())
        base_inject = SkewPolynomialBaseringInjection(base_ring, self)
        self._populate_coercion_lists_(
                coerce_list = [base_inject],
                convert_list = [list, base_inject])

    def _element_constructor_(self, a=None, check=True, construct=False, **kwds):
        r"""
        Convert a base ring element ``a`` into a constant of this univariate
        skew polynomial ring, possibly non-canonically.

        INPUT:

        - ``a`` -- (default: ``None``) an element of the base ring
          of ``self`` or a ring that has a coerce map from ``self``

        - ``check`` -- boolean (default: ``True``)

        - ``construct`` -- boolean (default: ``False``)

        OUTPUT:

        An zero-degree skew polynomial in ``self``, equal to ``a``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S(1 + x + x^2 + x^3)
            x^3 + x^2 + x + 1
            sage: S(1 + t)
            t + 1
            sage: S(1 + t).degree()
            0
            sage: S(0).list()
            []
        """
        C = self._polynomial_class
        if isinstance(a, list):
            return C(self, a, check=check, construct=construct)
        if isinstance(a, Element):
            P = a.parent()
            def build(check):
                if a.is_zero():
                    return P.zero()
                else:
                    return C(self, [a], check=check, construct=construct)
            if P is self:
                return a
            elif P is self.base_ring():
                build(False)
            elif P == self.base_ring() or self.base_ring().has_coerce_map_from(P):
                build(True)
        try:
            return a._polynomial_(self)
        except AttributeError:
            pass
        if isinstance(a, str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(a))
            except NameError:
                raise TypeError("unable to coerce string")
        return C(self, a, check, construct=construct, **kwds)

    def _coerce_map_from_(self, P):
        r"""
        Check whether ``self`` has a coerce map from ``P``.

        The rings that canonically coerce into this ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring

        - skew polynomial rings in the same variable and automorphism over
          any base ring that canonically coerces to the base ring of this ring

        INPUT:

        - ``P`` -- a ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.has_coerce_map_from(S)
            True
            sage: S.has_coerce_map_from(R)
            True
            sage: S.has_coerce_map_from(ZZ)
            True
            sage: S.has_coerce_map_from(GF(5^3))
            False

            sage: S.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
              Defn:   Polynomial base injection morphism:
                      From: Integer Ring
                      To:   Univariate Polynomial Ring in t over Integer Ring
                    then
                      Skew Polynomial base injection morphism:
                      From: Univariate Polynomial Ring in t over Integer Ring
                      To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S.coerce_map_from(S)
            Identity endomorphism of Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        try:
            connecting = self.base_ring().coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(self.base_ring()) * connecting
        except TypeError:
            pass
        try:
            if isinstance(P, SkewPolynomialRing_general):
                if self.__is_sparse and not P.is_sparse():
                    return False
                if P.variable_name() == self.variable_name():
                    if (P.base_ring() is self.base_ring()
                            and self.base_ring() is ZZ_sage):
                       if self._implementation_names == ('NTL',):
                            return False
                    return self.base_ring().has_coerce_map_from(P.base_ring())
        except AttributeError:
            pass

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        s = "Skew Polynomial Ring in %s over %s twisted by %s" % (self.variable_name(),
                                                                  self.base_ring(),
                                                                  self._map._repr_short())
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: latex(S)
            \Bold{Z}[t][x,\begin{array}{l}
            \text{\texttt{Ring{ }endomorphism...}}
            \end{array}]
        """
        from sage.misc.latex import latex
        return "%s[%s,%s]" % (latex(self.base_ring()), self.latex_variable_names()[0],
                              latex(self._map))

    def change_var(self, var):
        r"""
        Return the skew polynomial ring in variable ``var`` with the same base
        ring and twist map as ``self``.

        INPUT:

        - ``var`` -- a string representing the name of the new variable.

        OUTPUT:

        ``self`` with variable name name changed to ``var``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: R.<x> = SkewPolynomialRing(k,Frob); R
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ry = R.change_var('y'); Ry
            Skew Polynomial Ring in y over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ry is R.change_var('y')
            True
        """
        from sage.rings.polynomial.skew_polynomial_ring_constructor import SkewPolynomialRing
        return SkewPolynomialRing(self.base_ring(), self._map, names=var,
                                  sparse=self.__is_sparse)

    def characteristic(self):
        r"""
        Return the characteristic of the base ring of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: R['x',sigma].characteristic()
            0

            sage: k.<u> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: k['y',Frob].characteristic()
            5
        """
        return self.base_ring().characteristic()

    @cached_method
    def twist_map(self, n=1):
        r"""
        Return the twist map, the automorphism of the base ring of
        ``self``, iterated ``n`` times.

        INPUT:

        -  ``n`` - an integer (default: 1)

        OUTPUT:

        ``n``-th iterative of the twist map of this skew polynomial ring.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.twist_map()
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 1
            sage: S.twist_map() == sigma
            True
            sage: S.twist_map(10)
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 10

        If ``n`` in negative, Sage tries to compute the inverse of the
        twist map::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: T.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3

        Sometimes it fails, even if the twist map is actually invertible::

            sage: S.twist_map(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                  Defn: t |--> t + 1
        """
        try:
            return self._map ** n
        except TypeError as e:
            if n < 0:
                raise NotImplementedError("inversion of the twist map %s" % self._map)
            else:
                raise ValueError("Unexpected error in iterating the twist map: %s", e)

    @cached_method
    def gen(self, n=0):
        r"""
        Return the indeterminate generator of this skew polynomial ring.

        INPUT:

        - ``n`` -- index of generator to return (default: 0). Exists for
          compatibility with other polynomial rings.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
            sage: y = S.gen(); y
            x
            sage: y == x
            True
            sage: y is x
            True
            sage: S.gen(0)
            x

        This is also known as the parameter::

            sage: S.parameter() is S.gen()
            True
        """
        if n != 0:
            raise IndexError("generator %s not defined" % n)
        return self._polynomial_class(self, [0,1])

    parameter = gen

    def gens_dict(self):
        r"""
        Return a {name: variable} dictionary of the generators of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def is_finite(self):
        r"""
        Return ``False`` since skew polynomial rings are not finite
        (unless the base ring is `0`.)

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: k.is_finite()
            True
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_finite()
            False
        """
        R = self.base_ring()
        return R.is_finite() and R.order() == 1

    def is_exact(self):
        r"""
        Return ``True`` if elements of this skew polynomial ring are exact.
        This happens if and only if elements of the base ring are exact.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_exact()
            True
            sage: S.base_ring().is_exact()
            True

            sage: R.<u> = k[[]]
            sage: sigma = R.hom([u+u^2])
            sage: T.<y> = R['y',sigma]
            sage: T.is_exact()
            False
            sage: T.base_ring().is_exact()
            False
        """
        return self.base_ring().is_exact()

    def is_sparse(self):
        r"""
        Return ``True`` if the elements of this polynomial ring are sparsely
        represented.

        .. WARNING::

            Since sparse skew polynomials are not yet implemented, this
            function always returns ``False``.

        EXAMPLES:

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.is_sparse()
            False
        """
        return self.__is_sparse

    def ngens(self):
        r"""
        Return the number of generators of this skew polynomial ring,
        which is 1.

        EXAMPLES::

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.ngens()
            1
        """
        return 1

    def random_element(self, degree=2, monic=False, *args, **kwds):
        r"""
        Return a random skew polynomial in ``self``.

        INPUT:

        - ``degree`` -- (default: 2) integer with degree
          or a tuple of integers with minimum and maximum degrees

        - ``monic`` -- (default: ``False``) if ``True``, return a monic
          skew polynomial

        - ``*args, **kwds`` -- passed on to the ``random_element`` method
          for the base ring

        OUTPUT:

        Skew polynomial such that the coefficients of `x^i`, for `i` up
        to ``degree``, are random elements from the base ring, randomized
        subject to the arguments ``*args`` and ``**kwds``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: S.random_element()  # random
            (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: S.random_element(monic=True)  # random
            x^2 + (2*t^2 + t + 1)*x + 3*t^2 + 3*t + 2

        Use ``degree`` to obtain polynomials of higher degree

            sage: p = S.random_element(degree=5)   # random
            (t^2 + 3*t)*x^4 + (4*t + 4)*x^3 + (4*t^2 + 4*t)*x^2 + (2*t^2 + 1)*x + 3

        When ``monic`` is ``False``, the returned skew polynomial may have
        a degree less than ``degree`` (it happens when the random leading
        coefficient is zero). However, if ``monic`` is ``True``, this can't
        happen::

            sage: p = S.random_element(degree=4, monic=True)
            sage: p.leading_coefficient() == S.base_ring().one()
            True
            sage: p.degree() == 4
            True

        If a tuple of two integers is given for the degree argument, a random
        integer will be chosen between the first and second element of the
        tuple as the degree, both inclusive::

            sage: S.random_element(degree=(2,7))  # random
            (3*t^2 + 1)*x^4 + (4*t + 2)*x^3 + (4*t + 1)*x^2
             + (t^2 + 3*t + 3)*x + 3*t^2 + 2*t + 2

        If the first tuple element is greater than the second, a a
        ``ValueError`` is raised::

            sage: S.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: first degree argument must be less or equal to the second
        """
        R = self.base_ring()
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("first degree argument must be less or equal to the second")
            degree = randint(*degree)
        if monic:
            return self([R.random_element(*args, **kwds) for _ in range(degree)] + [R.one()])
        else:
            return self([R.random_element(*args, **kwds) for _ in range(degree+1)])

    def is_commutative(self):
        r"""
        Return ``True`` if this skew polynomial ring is commutative, i.e. if the
        twist map is the identity.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_commutative()
            False

            sage: T.<y> = k['y',Frob^3]
            sage: T.is_commutative()
            True
        """
        return self.twist_map().is_identity()
