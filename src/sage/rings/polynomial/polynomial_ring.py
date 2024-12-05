"""
Univariate Polynomial Rings

Sage implements sparse and dense polynomials over commutative and
non-commutative rings.  In the non-commutative case, the polynomial
variable commutes with the elements of the base ring.

AUTHOR:

- William Stein

- Kiran Kedlaya (2006-02-13): added macaulay2 option

- Martin Albrecht (2006-08-25): removed it again as it isn't needed anymore

- Simon King (2011-05): Dense and sparse polynomial rings must not be equal.

- Simon King (2011-10): Choice of categories for polynomial rings.

EXAMPLES::

    sage: z = QQ['z'].0
    sage: (z^3 + z - 1)^3
    z^9 + 3*z^7 - 3*z^6 + 3*z^5 - 6*z^4 + 4*z^3 - 3*z^2 + 3*z - 1

Saving and loading of polynomial rings works::

    sage: loads(dumps(QQ['x'])) == QQ['x']
    True
    sage: k = PolynomialRing(QQ['x'],'y'); loads(dumps(k))==k
    True
    sage: k = PolynomialRing(ZZ,'y'); loads(dumps(k)) == k
    True
    sage: k = PolynomialRing(ZZ,'y', sparse=True); loads(dumps(k))
    Sparse Univariate Polynomial Ring in y over Integer Ring

Rings with different variable names are not equal; in fact,
by :issue:`9944`, polynomial rings are equal if and only
if they are identical (which should be the  case for all parent
structures in Sage)::

    sage: QQ['y'] != QQ['x']
    True
    sage: QQ['y'] != QQ['z']
    True

We create a polynomial ring over a quaternion algebra::

    sage: # needs sage.combinat sage.modules
    sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
    sage: R.<w> = PolynomialRing(A, sparse=True)
    sage: f = w^3 + (i+j)*w + 1
    sage: f
    w^3 + (i + j)*w + 1
    sage: f^2
    w^6 + (2*i + 2*j)*w^4 + 2*w^3 - 2*w^2 + (2*i + 2*j)*w + 1
    sage: f = w + i ; g = w + j
    sage: f * g
    w^2 + (i + j)*w + k
    sage: g * f
    w^2 + (i + j)*w - k

:issue:`9944` introduced some changes related with
coercion. Previously, a dense and a sparse polynomial ring with the
same variable name over the same base ring evaluated equal, but of
course they were not identical. Coercion maps are cached - but if a
coercion to a dense ring is requested and a coercion to a sparse ring
is returned instead (since the cache keys are equal!), all hell breaks
loose.

Therefore, the coercion between rings of sparse and dense polynomials
works as follows::

    sage: R.<x> = PolynomialRing(QQ, sparse=True)
    sage: S.<x> = QQ[]
    sage: S == R
    False
    sage: S.has_coerce_map_from(R)
    True
    sage: R.has_coerce_map_from(S)
    False
    sage: (R.0 + S.0).parent()
    Univariate Polynomial Ring in x over Rational Field
    sage: (S.0 + R.0).parent()
    Univariate Polynomial Ring in x over Rational Field

It may be that one has rings of dense or sparse polynomials over
different base rings. In that situation, coercion works by means of
the :func:`~sage.categories.pushout.pushout` formalism::

    sage: R.<x> = PolynomialRing(GF(5), sparse=True)
    sage: S.<x> = PolynomialRing(ZZ)
    sage: R.has_coerce_map_from(S)
    False
    sage: S.has_coerce_map_from(R)
    False
    sage: S.0 + R.0
    2*x
    sage: (S.0 + R.0).parent()
    Univariate Polynomial Ring in x over Finite Field of size 5
    sage: (S.0 + R.0).parent().is_sparse()
    False

Similarly, there is a coercion from the (non-default) NTL
implementation for univariate polynomials over the integers
to the default FLINT implementation, but not vice versa::

    sage: R.<x> = PolynomialRing(ZZ, implementation='NTL')                              # needs sage.libs.ntl
    sage: S.<x> = PolynomialRing(ZZ, implementation='FLINT')
    sage: (S.0 + R.0).parent() is S                                                     # needs sage.libs.flint sage.libs.ntl
    True
    sage: (R.0 + S.0).parent() is S                                                     # needs sage.libs.flint sage.libs.ntl
    True

TESTS::

    sage: K.<x> = FractionField(QQ['x'])
    sage: V.<z> = K[]
    sage: x+z
    z + x

Check that :issue:`5562` has been fixed::

    sage: R.<u> = PolynomialRing(RDF, 1)
    sage: v1 = vector([u])                                                              # needs sage.modules
    sage: v2 = vector([CDF(2)])                                                         # needs sage.modules
    sage: v1 * v2                                                                       # needs sage.modules
    2.0*u
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import sys

from sage.structure.element import Element
from sage.structure.category_object import check_default_category

import sage.categories as categories
from sage.categories.morphism import IdentityMorphism
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.rings import Rings

from sage.rings.ring import (Ring, IntegralDomain)
from sage.structure.element import RingElement
import sage.rings.rational_field as rational_field
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.number_field.number_field_base import NumberField

try:
    from sage.libs.pari.all import pari_gen
except ImportError:
    pari_gen = ()

from sage.rings.polynomial.polynomial_ring_constructor import polynomial_default_category

import sage.misc.latex as latex
from sage.misc.prandom import randint
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

import sage.rings.abc
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.finite_rings.element_base import FiniteRingElement
from sage.rings.polynomial.polynomial_singular_interface import PolynomialRing_singular_repr
from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular
from sage.rings.power_series_ring_element import PowerSeries

_CommutativeRings = categories.commutative_rings.CommutativeRings()

import sage.interfaces.abc


def is_PolynomialRing(x):
    """
    Return ``True`` if ``x`` is a *univariate* polynomial ring (and not a
    sparse multivariate polynomial ring in one variable).

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        sage: from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
        sage: is_PolynomialRing(2)
        doctest:warning...
        DeprecationWarning: The function is_PolynomialRing is deprecated;
        use 'isinstance(..., PolynomialRing_general)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        False

    This polynomial ring is not univariate.

    ::

        sage: is_PolynomialRing(ZZ['x,y,z'])
        False
        sage: is_MPolynomialRing(ZZ['x,y,z'])
        doctest:warning...
        DeprecationWarning: The function is_MPolynomialRing is deprecated;
        use 'isinstance(..., MPolynomialRing_base)' instead.
        See https://github.com/sagemath/sage/issues/38266 for details.
        True

    ::

        sage: is_PolynomialRing(ZZ['w'])
        True

    Univariate means not only in one variable, but is a specific data
    type. There is a multivariate (sparse) polynomial ring data type,
    which supports a single variable as a special case.

    ::

        sage: # needs sage.libs.singular
        sage: R.<w> = PolynomialRing(ZZ, implementation='singular'); R
        Multivariate Polynomial Ring in w over Integer Ring
        sage: is_PolynomialRing(R)
        False
        sage: type(R)
        <class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
    """
    from sage.misc.superseded import deprecation
    deprecation(38266,
                "The function is_PolynomialRing is deprecated; "
                "use 'isinstance(..., PolynomialRing_general)' instead.")
    return isinstance(x, PolynomialRing_general)


#########################################################################################

class PolynomialRing_general(Ring):
    """
    Univariate polynomial ring over a ring.
    """

    def __init__(self, base_ring, name=None, sparse=False, implementation=None,
                 element_class=None, category=None):
        """
        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: R(-1) + R(1)
            0
            sage: (x - 2/3)*(x^2 - 8*x + 16)
            x^3 - 26/3*x^2 + 64/3*x - 32/3

            sage: category(ZZ['x'])
            Join of Category of unique factorization domains
             and Category of algebras with basis over
              (Dedekind domains and euclidean domains
               and noetherian rings and infinite enumerated sets
               and metric spaces)
             and Category of commutative algebras over
              (Dedekind domains and euclidean domains
               and noetherian rings and infinite enumerated sets
               and metric spaces)
             and Category of infinite sets

            sage: category(GF(7)['x'])
            Join of Category of euclidean domains
             and Category of algebras with basis over
              (finite enumerated fields and subquotients of monoids
               and quotients of semigroups)
            and Category of commutative algebras over
              (finite enumerated fields and subquotients of monoids
               and quotients of semigroups)
            and Category of infinite sets

        TESTS:

        Verify that :issue:`15232` has been resolved::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: TestSuite(R).run()

        Check that category for zero ring::

            sage: PolynomialRing(Zmod(1), 'x').category()
            Category of finite commutative rings

        Check `is_finite` inherited from category (:issue:`24432`)::

            sage: Zmod(1)['x'].is_finite()
            True

            sage: GF(7)['x'].is_finite()
            False

            sage: Zmod(1)['x']['y'].is_finite()
            True

            sage: GF(7)['x']['y'].is_finite()
            False
        """
        # We trust that, if category is given, it is useful and does not need to be joined
        # with the default category
        if base_ring.is_zero():
            category = categories.rings.Rings().Commutative().Finite()
        else:
            defaultcat = polynomial_default_category(base_ring.category(), 1)
            category = check_default_category(defaultcat, category)
        self.__is_sparse = sparse
        if element_class:
            self._polynomial_class = element_class
        else:
            if sparse:
                from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_sparse
                self._polynomial_class = Polynomial_generic_sparse
            else:
                from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense
                self._polynomial_class = Polynomial_generic_dense
        self.Element = self._polynomial_class
        self.__cyclopoly_cache = {}
        self._has_singular = False
        Ring.__init__(self, base_ring, names=name, normalize=True, category=category)
        from sage.rings.semirings.non_negative_integer_semiring import NonNegativeIntegerSemiring
        self._indices = NonNegativeIntegerSemiring()
        self._populate_coercion_lists_(convert_method_name='_polynomial_')

    def __reduce__(self):
        from sage.rings.polynomial.polynomial_ring_constructor import unpickle_PolynomialRing
        args = (self.base_ring(), self.variable_names(), None, self.is_sparse())
        return unpickle_PolynomialRing, args

    def _element_constructor_(self, x=None, check=True, is_gen=False,
                              construct=False, **kwds):
        r"""
        Convert ``x`` into this univariate polynomial ring,
        possibly non-canonically.

        Conversion from power series::

            sage: R.<x> = QQ[]
            sage: R(1 + x + x^2 + O(x^3))
            x^2 + x + 1

        Stacked polynomial rings coerce into constants if possible. First,
        the univariate case::

            sage: R.<x> = QQ[]
            sage: S.<u> = R[]
            sage: S(u + 2)
            u + 2
            sage: S(x + 3)
            x + 3
            sage: S(x + 3).degree()
            0

        Second, the multivariate case::

            sage: R.<x,y> = QQ[]
            sage: S.<u> = R[]
            sage: S(x + 2*y)
            x + 2*y
            sage: S(x + 2*y).degree()
            0
            sage: S(u + 2*x)
            u + 2*x
            sage: S(u + 2*x).degree()
            1

        Foreign polynomial rings coerce into the highest ring; the point
        here is that an element of T could coerce to an element of R or an
        element of S; it is anticipated that an element of T is more likely
        to be "the right thing" and is historically consistent.

        ::

            sage: R.<x> = QQ[]
            sage: S.<u> = R[]
            sage: T.<a> = QQ[]
            sage: S(a)
            u

        Coercing in pari elements::

            sage: QQ['x'](pari('[1,2,3/5]'))                                            # needs sage.libs.pari
            3/5*x^2 + 2*x + 1
            sage: QQ['x'](pari('(-1/3)*x^10 + (2/3)*x - 1/5'))                          # needs sage.libs.pari
            -1/3*x^10 + 2/3*x - 1/5

        Coercing strings::

            sage: QQ['y']('-y')
            -y

        TESTS:

        Python 3 range is allowed::

            sage: R = PolynomialRing(ZZ,'x')
            sage: R(range(4))
            3*x^3 + 2*x^2 + x

        This shows that the issue at :issue:`4106` is fixed::

            sage: # needs sage.symbolic
            sage: x = var('x')
            sage: R = IntegerModRing(4)
            sage: S = R['x']
            sage: S(x)
            x

        Throw a :exc:`TypeError` if any of the coefficients cannot be coerced
        into the base ring (:issue:`6777`)::

            sage: RealField(300)['x']( [ 1, ComplexField(300).gen(), 0 ])               # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: unable to convert '1.00...00*I' to a real number

        Check that the bug in :issue:`11239` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^2, prefix='z')
            sage: L.<b> = GF(5^4, prefix='z')
            sage: f = K['x'].gen() + a
            sage: L['x'](f)
            x + b^3 + b^2 + b + 3

        A test from :issue:`14485` ::

            sage: x = SR.var('x')                                                       # needs sage.symbolic
            sage: QQbar[x](x^6 + x^5 + x^4 - x^3 + x^2 - x + 2/5)                       # needs sage.rings.number_field sage.symbolic
            x^6 + x^5 + x^4 - x^3 + x^2 - x + 2/5

        Check support for unicode characters (:issue:`29280`)::

            sage: QQ['λ']('λ^2')
            λ^2
        """
        C = self.element_class
        if isinstance(x, (list, tuple)):
            return C(self, x, check=check, is_gen=False, construct=construct)
        if isinstance(x, range):
            return C(self, list(x), check=check, is_gen=False,
                     construct=construct)
        if isinstance(x, Element):
            P = x.parent()
            if P is self:
                return x
            elif P is self.base_ring():
                # It *is* the base ring, hence, we should not need to check.
                # Moreover, if x is equal to zero then we usually need to
                # provide [] to the polynomial class, not [x], if we don't want
                # to check (normally, polynomials like to strip trailing zeroes).
                # However, in the padic case, we WANT that trailing
                # zeroes are not stripped, because O(5)==0, but still it must
                # not be forgotten. It should be the job of the __init__ method
                # to decide whether to strip or not to strip.
                return C(self, [x], check=False, is_gen=False,
                         construct=construct)
            elif P == self.base_ring():
                return C(self, [x], check=True, is_gen=False,
                         construct=construct)
        if isinstance(x, sage.interfaces.abc.SingularElement) and self._has_singular:
            self._singular_().set_ring()
            try:
                return x.sage_poly(self)
            except Exception:
                raise TypeError("Unable to coerce singular object")
        elif isinstance(x, str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(x))
            except NameError:
                raise TypeError("Unable to coerce string")
        elif isinstance(x, FractionFieldElement):
            if x.denominator().is_unit():
                x = x.numerator() * x.denominator().inverse_of_unit()
            else:
                raise TypeError("denominator must be a unit")
        elif isinstance(x, pari_gen):
            if x.type() == 't_RFRAC':
                raise TypeError("denominator must be a unit")
            if x.type() != 't_POL':
                x = x.Polrev()
        elif isinstance(x, FiniteRingElement):
            try:
                return self(x.polynomial())
            except AttributeError:
                pass
        elif isinstance(x, PowerSeries):
            x = x.truncate()
        return C(self, x, check, is_gen, construct=construct, **kwds)

    @classmethod
    def _implementation_names(cls, implementation, base_ring, sparse=False):
        """
        Check whether this class can handle the implementation
        ``implementation`` over the given base ring and sparseness.

        This is a simple wrapper around :meth:`_implementation_names_impl`
        which does the real work.

        .. NOTE::

            It is assumed that the ``base_ring`` argument is a base ring
            which the class can handle.

        INPUT:

        - ``implementation`` -- either a string denoting a specific
          implementation or ``None`` for the default

        - ``base_ring`` -- the base ring for the polynomial ring

        - ``sparse`` -- boolean; whether the implementation is sparse

        OUTPUT:

        - if the implementation is supported, the output is a list of
          all names (possibly including ``None``) which refer to the
          given implementation. The first element of the list is the
          canonical name. If the ``__init__`` method does not take an
          ``implementation`` keyword, then the first element must be
          ``None``.

        - if the implementation is not supported, raise a
          :exc:`ValueError`.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
            sage: PolynomialRing_general._implementation_names(None, ZZ, True)
            [None, 'generic']
            sage: PolynomialRing_general._implementation_names("generic", ZZ, True)
            [None, 'generic']
            sage: PolynomialRing_general._implementation_names("xyzzy", ZZ, True)
            Traceback (most recent call last):
            ...
            ValueError: unknown implementation 'xyzzy' for sparse polynomial rings over Integer Ring
        """
        names = cls._implementation_names_impl(implementation, base_ring, sparse)
        if names is NotImplemented:
            raise ValueError("unknown implementation %r for %s polynomial rings over %r" %
                    (implementation, "sparse" if sparse else "dense", base_ring))
        assert isinstance(names, list)
        assert implementation in names
        return names

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        The underlying implementation of :meth:`_implementation_names`.

        The behaviour is exactly the same, except that an unknown
        implementation returns ``NotImplemented`` instead of raising
        :exc:`ValueError`.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
            sage: PolynomialRing_general._implementation_names_impl("xyzzy", ZZ, True)
            NotImplemented
        """
        if implementation is None or implementation == "generic":
            return [None, "generic"]
        return NotImplemented

    def is_integral_domain(self, proof=True):
        """
        EXAMPLES::

            sage: ZZ['x'].is_integral_domain()
            True
            sage: Integers(8)['x'].is_integral_domain()
            False
        """
        return self.base_ring().is_integral_domain(proof)

    def is_unique_factorization_domain(self, proof=True):
        """
        EXAMPLES::

            sage: ZZ['x'].is_unique_factorization_domain()
            True
            sage: Integers(8)['x'].is_unique_factorization_domain()
            False
        """
        return self.base_ring().is_unique_factorization_domain(proof)

    def is_noetherian(self):
        return self.base_ring().is_noetherian()

    def some_elements(self):
        r"""
        Return a list of polynomials.

        This is typically used for running generic tests.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.some_elements()
            [x, 0, 1, 1/2, x^2 + 2*x + 1, x^3, x^2 - 1, x^2 + 1, 2*x^2 + 2]
        """
        # the comments in the following lines describe the motivation for
        # adding these elements, they are not accurate over all rings and in
        # all contexts
        R = self.base_ring()
        # Doing things this way is a little robust against rings where
        #    2 might not convert in
        one = R.one()
        return [self.gen(),
            self.zero(), self(one), self(R.an_element()), # elements of the base ring
            self([one,2*one,one]), # a square
            self([0,0,0,one]), # a power but not a square
            self([-one,0,one]), # a reducible element
            self([one,0,one]), # an irreducible element
            self([2*one,0,2*one]), # an element with non-trivial content
        ]

    @cached_method
    def flattening_morphism(self):
        r"""
        Return the flattening morphism of this polynomial ring.

        EXAMPLES::

            sage: QQ['a','b']['x'].flattening_morphism()
            Flattening morphism:
              From: Univariate Polynomial Ring in x over
                    Multivariate Polynomial Ring in a, b over Rational Field
              To:   Multivariate Polynomial Ring in a, b, x over Rational Field

            sage: QQ['x'].flattening_morphism()
            Identity endomorphism of Univariate Polynomial Ring in x over Rational Field
        """
        from .multi_polynomial_ring import MPolynomialRing_base
        base = self.base_ring()
        if isinstance(base, PolynomialRing_general) or isinstance(base, MPolynomialRing_base):
            from .flatten import FlatteningMorphism
            return FlatteningMorphism(self)
        else:
            return IdentityMorphism(self)

    def construction(self):
        """
        Return the construction functor.
        """
        return categories.pushout.PolynomialFunctor(self.variable_name(), sparse=self.__is_sparse), self.base_ring()

    def completion(self, p=None, prec=20, extras=None):
        r"""
        Return the completion of ``self`` with respect to the irreducible
        polynomial ``p``.

        Currently only implemented for ``p=self.gen()`` (the default), i.e. you
        can only complete `R[x]` with respect to `x`, the result being a ring
        of power series in `x`. The ``prec`` variable controls the precision
        used in the power series ring. If ``prec`` is `\infty`, then this
        returns a :class:`LazyPowerSeriesRing`.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: P
            Univariate Polynomial Ring in x over Rational Field
            sage: PP = P.completion(x)
            sage: PP
            Power Series Ring in x over Rational Field
            sage: f = 1 - x
            sage: PP(f)
            1 - x
            sage: 1 / f
            -1/(x - 1)
            sage: g = 1 / PP(f); g
            1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10 + x^11
             + x^12 + x^13 + x^14 + x^15 + x^16 + x^17 + x^18 + x^19 + O(x^20)
            sage: 1 / g
            1 - x + O(x^20)

            sage: # needs sage.combinat
            sage: PP = P.completion(x, prec=oo); PP
            Lazy Taylor Series Ring in x over Rational Field
            sage: g = 1 / PP(f); g
            1 + x + x^2 + O(x^3)
            sage: 1 / g == f
            True
        """
        if p is None or str(p) == self._names[0]:
            if prec == float('inf'):
                from sage.rings.lazy_series_ring import LazyPowerSeriesRing
                return LazyPowerSeriesRing(self.base_ring(), names=(self._names[0],),
                                           sparse=self.is_sparse())
            from sage.rings.power_series_ring import PowerSeriesRing
            return PowerSeriesRing(self.base_ring(), name=self._names[0],
                                   default_prec=prec, sparse=self.is_sparse())

        raise NotImplementedError("cannot complete %s with respect to %s" % (self, p))

    def _coerce_map_from_base_ring(self):
        """
        Return a coercion map from the base ring of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.coerce_map_from(QQ)
            Polynomial base injection morphism:
              From: Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field
            sage: R.coerce_map_from(ZZ)
            Composite map:
              From: Integer Ring
              To:   Univariate Polynomial Ring in x over Rational Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Polynomial base injection morphism:
                      From: Rational Field
                      To:   Univariate Polynomial Ring in x over Rational Field
            sage: R.coerce_map_from(GF(7))
        """
        from sage.rings.polynomial.polynomial_element import PolynomialBaseringInjection

        return PolynomialBaseringInjection(self.base_ring(), self)

    def _coerce_map_from_(self, P):
        """
        The rings that canonically coerce to this polynomial ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring.

        - polynomial rings in the same variable over any base ring that
          canonically coerces to the base ring of this ring.

        - a multivariate polynomial ring P such that ``self``'s variable name
          is among the variable names of P, and the ring obtained by
          removing that variable is different from the base ring of ``self``,
          but coerces into it. (see :issue:`813` for a discussion of this)

        Caveat: There is no coercion from a dense into a sparse
        polynomial ring. So, when adding a dense and a sparse
        polynomial, the result will be dense. See :issue:`9944`.

        EXAMPLES::

            sage: R = QQ['x']
            sage: R.has_coerce_map_from(ZZ['x'])
            True
            sage: R.has_coerce_map_from(ZZ['y'])
            False

        Here we test against the change in the coercions introduced
        in :issue:`9944`::

            sage: R.<x> = PolynomialRing(QQ, sparse=True)
            sage: S.<x> = QQ[]
            sage: (R.0 + S.0).parent()
            Univariate Polynomial Ring in x over Rational Field
            sage: (S.0 + R.0).parent()
            Univariate Polynomial Ring in x over Rational Field

        Here we test a feature that was implemented in :issue:`813`::

            sage: P = QQ['x','y']
            sage: Q = Frac(QQ['x'])['y']
            sage: Q.has_coerce_map_from(P)
            True
            sage: P.0 + Q.0
            y + x

        In order to avoid bidirectional coercions (which are generally
        problematic), we only have a coercion from P to Q if the base
        ring of Q is more complicated than "P minus one variable"::

            sage: Q = QQ['x']['y']
            sage: P.has_coerce_map_from(Q)
            True
            sage: Q.has_coerce_map_from(P)
            False
            sage: Q.base_ring() is P.remove_var(Q.variable_name())
            True

        Over the integers, there is a coercion from the NTL and generic
        implementation to the default FLINT implementation::

            sage: del R, S  # clear values from doctests above
            sage: R = PolynomialRing(ZZ, 't', implementation='NTL')                     # needs sage.libs.ntl
            sage: S = PolynomialRing(ZZ, 't', implementation='FLINT')                   # needs sage.libs.flint
            sage: T = PolynomialRing(ZZ, 't', implementation='generic')
            sage: R.has_coerce_map_from(S)                                              # needs sage.libs.flint sage.libs.ntl
            False
            sage: R.has_coerce_map_from(T)                                              # needs sage.libs.ntl
            False
            sage: S.has_coerce_map_from(T)                                              # needs sage.libs.flint
            True
            sage: S.has_coerce_map_from(R)                                              # needs sage.libs.flint sage.libs.ntl
            True
            sage: T.has_coerce_map_from(R)                                              # needs sage.libs.ntl
            False
            sage: T.has_coerce_map_from(S)                                              # needs sage.libs.flint
            False
        """
        base_ring = self.base_ring()

        # workaround, useful for the zero ring
        if P == base_ring:
            return self._coerce_map_from_base_ring()

        # handle constants that canonically coerce into self.base_ring()
        # first, if possible
        try:
            connecting = base_ring.coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(base_ring) * connecting
        except TypeError:
            pass

        # polynomial rings in the same variable over a base that canonically
        # coerces into self.base_ring()
        if isinstance(P, PolynomialRing_general):
            if self.construction()[0] != P.construction()[0]:
                # Construction (including variable names) must be the
                # same to allow coercion
                return False
            self_sparse = self.is_sparse()
            P_sparse = P.is_sparse()
            if self_sparse and not P_sparse:
                # Coerce only sparse -> dense
                return False

            if P.base_ring() is base_ring:
                # Same base ring but different implementations.
                # Ideally, we should avoid cyclic coercions (a coercion
                # from A to B and also from B to A), but this is
                # currently hard to do:
                # see https://github.com/sagemath/sage/issues/24319
                if not self_sparse and P_sparse:
                    # Always allow coercion sparse -> dense
                    pass
                elif base_ring is ZZ:
                    # Over ZZ, only allow coercion from any ZZ['x']
                    # implementation to the default FLINT implementation
                    try:
                        from .polynomial_integer_dense_flint import Polynomial_integer_dense_flint
                    except ImportError:
                        return None
                    if self.element_class is not Polynomial_integer_dense_flint:
                        return None
                # Other rings: always allow coercion
                # To be fixed in Issue #24319
            f = base_ring.coerce_map_from(P.base_ring())
            if f is None:
                return None
            from sage.rings.homset import RingHomset
            from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
            return PolynomialRingHomomorphism_from_base(RingHomset(P, self), f)

        # Last, we consider multivariate polynomial rings:
        from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
        if isinstance(P, MPolynomialRing_base) and self.variable_name() in P.variable_names():
            P_ = P.remove_var(self.variable_name())
            return self.base_ring() != P_ and self.base_ring().has_coerce_map_from(P_)

    def _magma_init_(self, magma):
        """
        Used in converting this ring to the corresponding ring in MAGMA.

        EXAMPLES::

            sage: # optional - magma
            sage: R = QQ['y']
            sage: R._magma_init_(magma)
            'SageCreateWithNames(PolynomialRing(_sage_ref...),["y"])'
            sage: S = magma(R)
            sage: S
            Univariate Polynomial Ring in y over Rational Field
            sage: S.1
            y
            sage: magma(PolynomialRing(GF(7), 'x'))                                     # needs sage.rings.finite_rings
            Univariate Polynomial Ring in x over GF(7)
            sage: magma(PolynomialRing(GF(49,'a'), 'x'))                                # needs sage.rings.finite_rings
            Univariate Polynomial Ring in x over GF(7^2)
            sage: magma(PolynomialRing(PolynomialRing(ZZ,'w'), 'x'))
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in w over Integer Ring

        Watch out, Magma has different semantics from Sage, i.e., in Magma
        there is a unique univariate polynomial ring, and the variable name
        has no intrinsic meaning (it only impacts printing), so can't be
        reliably set because of caching.

        ::

            sage: # optional - magma
            sage: m = Magma()
            sage: m(QQ['w'])
            Univariate Polynomial Ring in w over Rational Field
            sage: m(QQ['x'])
            Univariate Polynomial Ring in x over Rational Field
            sage: m(QQ['w'])
            Univariate Polynomial Ring in x over Rational Field

        A nested example over a Givaro finite field::

            sage: k.<a> = GF(9)                                                         # needs sage.rings.finite_rings
            sage: R.<x> = k[]                                                           # needs sage.rings.finite_rings
            sage: magma(a^2*x^3 + (a+1)*x + a)  # optional - magma                      # needs sage.rings.finite_rings
            a^2*x^3 + a^2*x + a
        """
        B = magma(self.base_ring())
        Bref = B._ref()
        s = 'PolynomialRing(%s)' % (Bref)
        return magma._with_names(s, self.variable_names())

    def _gap_init_(self, gap=None):
        """
        String for representing this polynomial ring in GAP.

        INPUT:

        - ``gap`` -- (optional GAP instance) used for representing the base ring

        EXAMPLES::

            sage: R.<z> = ZZ[]
            sage: gap(R)                                                                # needs sage.libs.gap
            PolynomialRing( Integers, ["z"] )
            sage: gap(R) is gap(R)                                                      # needs sage.libs.gap
            True
            sage: gap(z^2 + z)                                                          # needs sage.libs.gap
            z^2+z

        A univariate polynomial ring over a multivariate polynomial
        ring over a number field::

            sage: # needs sage.rings.number_field
            sage: Q.<t> = QQ[]
            sage: K.<tau> = NumberField(t^2 + t + 1)
            sage: P.<x,y> = K[]
            sage: S.<z> = P[]
            sage: gap(S)                                                                # needs sage.libs.gap
            PolynomialRing( PolynomialRing( <algebraic extension over the Rationals of degree 2>, ["x", "y"] ), ["z"] )
            sage: gap(S) is gap(S)                                                      # needs sage.libs.gap
            True
        """
        if gap is not None:
            base_ring = gap(self.base_ring()).name()
        else:
            base_ring = self.base_ring()._gap_init_()
        return 'PolynomialRing(%s, ["%s"])' % (base_ring, self.variable_name())

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: sage_input(GF(5)['x']['y'], verify=True)
            # Verified
            GF(5)['x']['y']
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: ZZ['z']._sage_input_(SageInputBuilder(), False)
            {constr_parent: {subscr: {atomic:ZZ}[{atomic:'z'}]} with gens: ('z',)}
        """
        base = sib(self.base_ring())
        sie = base[self.variable_name()]
        gens_syntax = sib.empty_subscript(base)
        return sib.parent_with_gens(self, sie, self.variable_names(), 'R',
                                    gens_syntax=gens_syntax)

    def _macaulay2_init_(self, macaulay2=None):
        """
        EXAMPLES::

            sage: R = QQ['x']
            sage: macaulay2(R).describe()  # optional - macaulay2
            QQ[x, Degrees => {1}, Heft => {1}, MonomialOrder => {MonomialSize => 16},
                                                                {GRevLex => {1}    }
                                                                {Position => Up    }
            --------------------------------------------------------------------------------
            DegreeRank => 1]

        TESTS:

        Check that results are cached (:issue:`28074`)::

            sage: R = ZZ['t']
            sage: macaulay2(R) is macaulay2(R)  # optional - macaulay2
            True
        """
        if macaulay2 is None:
            from sage.interfaces.macaulay2 import macaulay2 as m2_default
            macaulay2 = m2_default
        return macaulay2._macaulay2_input_ring(self.base_ring(), self.gens())

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R._is_valid_homomorphism_(GF(7), [5])
            False
            sage: R._is_valid_homomorphism_(Qp(7), [5])                                 # needs sage.rings.padics
            True
        """
        # Since poly rings are free, any image of the gen
        # determines a homomorphism
        if base_map is None:
            # If no base map is given, the only requirement is that the
            # base ring coerces into the codomain
            return codomain.has_coerce_map_from(self.base_ring())
        return True

    #    Polynomial rings should be unique parents. Hence,
    #    no need for any comparison method

    def __hash__(self):
        # should be faster than just relying on the string representation
        try:
            return self._cached_hash
        except AttributeError:
            pass
        h = self._cached_hash = hash((self.base_ring(),self.variable_name()))
        return h

    def _repr_(self):
        try:
            return self._cached_repr
        except AttributeError:
            pass
        s = "Univariate Polynomial Ring in %s over %s" % (
                self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = "Sparse " + s
        self._cached_repr = s
        return s

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: S.<alpha12> = ZZ[]
            sage: latex(S)
            \Bold{Z}[\alpha_{12}]
        """
        return "%s[%s]" % (latex.latex(self.base_ring()), self.latex_variable_names()[0])

    def base_extend(self, R):
        """
        Return the base extension of this polynomial ring to `R`.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: R.<x> = RR[]; R
            Univariate Polynomial Ring in x over Real Field with 53 bits of precision
            sage: R.base_extend(CC)
            Univariate Polynomial Ring in x over Complex Field with 53 bits of precision
            sage: R.base_extend(QQ)
            Traceback (most recent call last):
            ...
            TypeError: no such base extension
            sage: R.change_ring(QQ)
            Univariate Polynomial Ring in x over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if R.has_coerce_map_from(self.base_ring()):
            return PolynomialRing(R, names=self.variable_name(), sparse=self.is_sparse())
        else:
            raise TypeError("no such base extension")

    def change_ring(self, R):
        """
        Return the polynomial ring in the same variable as ``self`` over `R`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings sage.rings.real_interval_field
            sage: R.<ZZZ> = RealIntervalField()[]; R
            Univariate Polynomial Ring in ZZZ over
             Real Interval Field with 53 bits of precision
            sage: R.change_ring(GF(19^2, 'b'))
            Univariate Polynomial Ring in ZZZ over Finite Field in b of size 19^2
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        return PolynomialRing(R, names=self.variable_name(), sparse=self.is_sparse())

    def change_var(self, var):
        r"""
        Return the polynomial ring in variable ``var`` over the same base
        ring.

        EXAMPLES::

            sage: R.<x> = ZZ[]; R
            Univariate Polynomial Ring in x over Integer Ring
            sage: R.change_var('y')
            Univariate Polynomial Ring in y over Integer Ring
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        return PolynomialRing(self.base_ring(), names=var, sparse=self.is_sparse())

    def extend_variables(self, added_names, order='degrevlex'):
        r"""
        Return a multivariate polynomial ring with the same base ring but
        with ``added_names`` as additional variables.

        EXAMPLES::

            sage: R.<x> = ZZ[]; R
            Univariate Polynomial Ring in x over Integer Ring
            sage: R.extend_variables('y, z')
            Multivariate Polynomial Ring in x, y, z over Integer Ring
            sage: R.extend_variables(('y', 'z'))
            Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        if isinstance(added_names, str):
            added_names = added_names.split(',')
        return PolynomialRing(self.base_ring(), names=self.variable_names() + tuple(added_names), order=order)

    def variable_names_recursive(self, depth=sage.rings.infinity.infinity):
        r"""
        Return the list of variable names of this ring and its base rings,
        as if it were a single multi-variate polynomial.

        INPUT:

        - ``depth`` -- integer or :mod:`Infinity <sage.rings.infinity>`

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: R = QQ['x']['y']['z']
            sage: R.variable_names_recursive()
            ('x', 'y', 'z')
            sage: R.variable_names_recursive(2)
            ('y', 'z')
        """
        if depth <= 0:
            return ()
        elif depth == 1:
            return self.variable_names()
        else:
            my_vars = self.variable_names()
            try:
                return self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
            except AttributeError:
                return my_vars

    def _mpoly_base_ring(self, variables=None):
        r"""
        Return the base ring if this is viewed as a polynomial ring over
        ``variables``. See also ``Polynomial._mpoly_dict_recursive``.
        """
        if variables is None:
            variables = self.variable_names_recursive()
        variables = list(variables)
        var = self.variable_name()
        if var not in variables:
            return self
        else:
            try:
                return self.base_ring()._mpoly_base_ring(variables[:variables.index(var)])
            except AttributeError:
                return self.base_ring()

    def characteristic(self):
        """
        Return the characteristic of this polynomial ring, which is the
        same as that of its base ring.

        EXAMPLES::

            sage: # needs sage.rings.real_interval_field
            sage: R.<ZZZ> = RealIntervalField()[]; R
            Univariate Polynomial Ring in ZZZ over Real Interval Field with 53 bits of precision
            sage: R.characteristic()
            0
            sage: S = R.change_ring(GF(19^2, 'b')); S                                   # needs sage.rings.finite_rings
            Univariate Polynomial Ring in ZZZ over Finite Field in b of size 19^2
            sage: S.characteristic()                                                    # needs sage.rings.finite_rings
            19
        """
        return self.base_ring().characteristic()

    def cyclotomic_polynomial(self, n):
        """
        Return the `n`-th cyclotomic polynomial as a polynomial in this
        polynomial ring. For details of the implementation, see the
        documentation for
        :func:`sage.rings.polynomial.cyclotomic.cyclotomic_coeffs`.

        EXAMPLES::

            sage: R = ZZ['x']
            sage: R.cyclotomic_polynomial(8)
            x^4 + 1
            sage: R.cyclotomic_polynomial(12)
            x^4 - x^2 + 1

            sage: S = PolynomialRing(FiniteField(7), 'x')
            sage: S.cyclotomic_polynomial(12)
            x^4 + 6*x^2 + 1
            sage: S.cyclotomic_polynomial(1)
            x + 6

        TESTS:

        Make sure it agrees with other systems for the trivial case::

            sage: ZZ['x'].cyclotomic_polynomial(1)
            x - 1
            sage: gp('polcyclo(1)')                                                     # needs sage.libs.pari
            x - 1
        """
        if n <= 0:
            raise ArithmeticError("n=%s must be positive" % n)
        elif n == 1:
            return self.gen() - 1
        else:
            from .cyclotomic import cyclotomic_coeffs
            return self(cyclotomic_coeffs(n), check=True)

    @cached_method
    def gen(self, n=0):
        """
        Return the indeterminate generator of this polynomial ring.

        EXAMPLES::

            sage: R.<abc> = Integers(8)[]; R
            Univariate Polynomial Ring in abc over Ring of integers modulo 8
            sage: t = R.gen(); t
            abc
            sage: t.is_gen()
            True

        An identical generator is always returned.

        ::

            sage: t is R.gen()
            True
        """
        if n != 0:
            raise IndexError("generator n not defined")
        return self.element_class(self, [0,1], is_gen=True)

    def gens_dict(self):
        """
        Return a dictionary whose entries are ``{name:variable,...}``,
        where ``name`` stands for the variable names of this
        object (as strings) and ``variable`` stands for the corresponding
        generators (as elements of this object).

        EXAMPLES::

            sage: R.<y,x,a42> = RR[]
            sage: R.gens_dict()
            {'a42': a42, 'x': x, 'y': y}
        """
        gens = self.gens()
        names = self.variable_names()
        assert len(gens) == len(names)
        return dict(zip(names, gens))

    def parameter(self):
        """
        Return the generator of this polynomial ring.

        This is the same as ``self.gen()``.
        """
        return self.gen()

    @cached_method
    def is_exact(self):
        return self.base_ring().is_exact()

    def is_field(self, proof=True):
        """
        Return ``False``, since polynomial rings are never fields.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<z> = Integers(2)[]; R
            Univariate Polynomial Ring in z over Ring of integers modulo 2 (using GF2X)
            sage: R.is_field()
            False
        """
        return False

    def is_sparse(self):
        """
        Return ``True`` if elements of this polynomial ring have a sparse
        representation.

        EXAMPLES::

            sage: R.<z> = Integers(8)[]; R
            Univariate Polynomial Ring in z over Ring of integers modulo 8
            sage: R.is_sparse()
            False
            sage: R.<W> = PolynomialRing(QQ, sparse=True); R
            Sparse Univariate Polynomial Ring in W over Rational Field
            sage: R.is_sparse()
            True
        """
        return self.__is_sparse

    def monomial(self, exponent):
        """
        Return the monomial with the ``exponent``.

        INPUT:

        - ``exponent`` -- nonnegative integer

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: R.monomial(5)
            x^5
            sage: e=(10,)
            sage: R.monomial(*e)
            x^10
            sage: m = R.monomial(100)
            sage: R.monomial(m.degree()) == m
            True
        """
        return self({exponent: self.base_ring().one()})

    def krull_dimension(self):
        """
        Return the Krull dimension of this polynomial ring, which is one
        more than the Krull dimension of the base ring.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.krull_dimension()
            1

            sage: # needs sage.rings.finite_rings
            sage: R.<z> = GF(9, 'a')[]; R
            Univariate Polynomial Ring in z over Finite Field in a of size 3^2
            sage: R.krull_dimension()
            1
            sage: S.<t> = R[]
            sage: S.krull_dimension()
            2
            sage: for n in range(10):
            ....:     S = PolynomialRing(S, 'w')
            sage: S.krull_dimension()
            12
        """
        return self.base_ring().krull_dimension() + 1

    def ngens(self):
        """
        Return the number of generators of this polynomial ring, which is 1
        since it is a univariate polynomial ring.

        EXAMPLES::

            sage: R.<z> = Integers(8)[]; R
            Univariate Polynomial Ring in z over Ring of integers modulo 8
            sage: R.ngens()
            1
        """
        return 1

    def random_element(self, degree=(-1, 2), monic=False, *args, **kwds):
        r"""
        Return a random polynomial of given degree (bounds).

        INPUT:

        - ``degree`` -- (default: ``(-1, 2)``) integer for fixing the degree or
          a tuple of minimum and maximum degrees

        - ``monic`` -- boolean (default: ``False``); indicate whether the sampled
          polynomial should be monic

        - ``*args, **kwds`` -- additional keyword parameters passed on to the
          ``random_element`` method for the base ring

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = R.random_element(10, x=5, y=10)
            sage: f.degree()
            10
            sage: f.parent() is R
            True
            sage: all(a in range(5, 10) for a in f.coefficients())
            True
            sage: R.random_element(6).degree()
            6

        If a tuple of two integers is given for the ``degree`` argument, a
        polynomial is chosen among all polynomials with degree between them. If
        the base ring can be sampled uniformly, then this method also samples
        uniformly::

            sage: R.random_element(degree=(0, 4)).degree() in range(0, 5)
            True
            sage: found = [False]*5
            sage: while not all(found):
            ....:     found[R.random_element(degree=(0, 4)).degree()] = True

        Note that the zero polynomial has degree `-1`, so if you want to
        consider it set the minimum degree to `-1`::

            sage: while R.random_element(degree=(-1,2), x=-1, y=1) != R.zero():
            ....:     pass

        Monic polynomials are chosen among all monic polynomials with degree
        between the given ``degree`` argument::

            sage: all(R.random_element(degree=(-1, 1), monic=True).is_monic() for _ in range(10^3))
            True
            sage: all(R.random_element(degree=(0, 1), monic=True).is_monic() for _ in range(10^3))
            True

        TESTS::

            sage: R.random_element(degree=[5])
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)

            sage: R.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: minimum degree must be less or equal than maximum degree

        Check that :issue:`16682` is fixed::

            sage: R = PolynomialRing(GF(2), 'z')
            sage: for _ in range(100):
            ....:     d = randint(-1, 20)
            ....:     P = R.random_element(degree=d)
            ....:     assert P.degree() == d

        In :issue:`37118`, ranges including integers below `-1` no longer raise
        an error::

            sage: R.random_element(degree=(-2, 3))  # random
            z^3 + z^2 + 1

        ::

            sage: 0 in [R.random_element(degree=(-1, 2), monic=True) for _ in range(500)]
            False

        Testing error handling::

            sage: R.random_element(degree=-5)
            Traceback (most recent call last):
            ...
            ValueError: degree (=-5) must be at least -1

            sage: R.random_element(degree=(-3, -2))
            Traceback (most recent call last):
            ...
            ValueError: maximum degree (=-2) must be at least -1

        Testing uniformity::

            sage: from collections import Counter
            sage: R = GF(3)["x"]
            sage: samples = [R.random_element(degree=(-1, 2)) for _ in range(27000)]    # long time
            sage: assert all(750 <= f <= 1250 for f in Counter(samples).values())       # long time

            sage: samples = [R.random_element(degree=(-1, 2), monic=True) for _ in range(13000)] # long time
            sage: assert all(750 <= f <= 1250 for f in Counter(samples).values())       # long time
        """
        R = self.base_ring()

        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("minimum degree must be less or equal than maximum degree")
            if degree[1] < -1:
                raise ValueError(f"maximum degree (={degree[1]}) must be at least -1")
        else:
            if degree < -1:
                raise ValueError(f"degree (={degree}) must be at least -1")
            degree = (degree, degree)

        if degree[0] <= -2:
            degree = (-1, degree[1])

        # If the coefficient range only contains 0, then
        # * if the degree range includes -1, return the zero polynomial,
        # * otherwise raise a value error
        if args == (0, 1):
            if degree[0] == -1:
                return self.zero()
            else:
                raise ValueError("No polynomial of degree >= 0 has all coefficients zero")

        if degree == (-1, -1):
            return self.zero()

        # If `monic` is set, zero should be ignored
        if degree[0] == -1 and monic:
            if degree[1] == -1:
                raise ValueError("the maximum degree of monic polynomials needs to be at least 0")
            if degree[1] == 0:
                return self.one()
            degree = (0, degree[1])

        # Pick random coefficients
        end = degree[1]
        if degree[0] == -1:
            return self([R.random_element(*args, **kwds) for _ in range(end + 1)])

        nonzero = False
        coefs = [None] * (end + 1)

        while not nonzero:
            # Pick leading coefficients, if `monic` is set it's handle here.
            if monic:
                for i in range(degree[1] - degree[0] + 1):
                    coefs[end - i] = R.random_element(*args, **kwds)
                    if not nonzero and not coefs[end - i].is_zero():
                        coefs[end - i] = R.one()
                        nonzero = True
            else:
                # Fast path
                for i in range(degree[1] - degree[0] + 1):
                    coefs[end - i] = R.random_element(*args, **kwds)
                    nonzero |= not coefs[end - i].is_zero()

        # Now we pick the remaining coefficients.
        for i in range(degree[1] - degree[0] + 1, degree[1] + 1):
            coefs[end - i] = R.random_element(*args, **kwds)

        return self(coefs)

    def _monics_degree(self, of_degree):
        """
        Refer to monics() for full documentation.
        """
        base = self.base_ring()
        for coeffs in sage.misc.mrange.xmrange_iter([[base.one()]]+[base]*of_degree):
            # Each iteration returns a *new* list!
            # safe to mutate the return
            coeffs.reverse()
            yield self(coeffs)

    def _monics_max(self, max_degree):
        """
        Refer to monics() for full documentation.
        """
        for degree in range(max_degree + 1):
            yield from self._monics_degree(degree)

    def _polys_degree(self, of_degree):
        """
        Refer to polynomials() for full documentation.
        """
        base = self.base_ring()
        base0 = base.zero()
        for leading_coeff in base:
            if leading_coeff != base0:
                for lt1 in sage.misc.mrange.xmrange_iter([base]*(of_degree)):
                    # Each iteration returns a *new* list!
                    # safe to mutate the return
                    coeffs = [leading_coeff] + lt1
                    coeffs.reverse()
                    yield self(coeffs)

    def _polys_max( self, max_degree ):
        """
        Refer to polynomials() for full documentation.
        """
        base = self.base_ring()
        for coeffs in sage.misc.mrange.xmrange_iter([base]*(max_degree+1)):
            # Each iteration returns a *new* list!
            # safe to mutate the return
            coeffs.reverse()
            yield self(coeffs)

    @lazy_attribute
    def _Karatsuba_threshold(self):
        """
        Return the default Karatsuba threshold.

        EXAMPLES::

            sage: R.<x> = QQbar[]                                                       # needs sage.rings.number_field
            sage: R._Karatsuba_threshold                                                # needs sage.rings.number_field
            8
            sage: MS = MatrixSpace(ZZ, 2, 2)                                            # needs sage.modules
            sage: R.<x> = MS[]                                                          # needs sage.modules
            sage: R._Karatsuba_threshold                                                # needs sage.modules
            0
        """
        base_ring = self.base_ring()
        if isinstance(base_ring, PolynomialRing_general):
            return 0
        try:
            from sage.matrix.matrix_space import MatrixSpace
        except ImportError:
            pass
        else:
            if isinstance(base_ring, MatrixSpace):
                return 0
        from sage.rings.fraction_field import FractionField_generic
        if isinstance(base_ring, FractionField_generic):
            return 1 << 60
        # Generic default value
        return 8

    def karatsuba_threshold(self):
        """
        Return the Karatsuba threshold used for this ring by the method
        :meth:`_mul_karatsuba` to fall back to the schoolbook algorithm.

        EXAMPLES::

            sage: K = QQ['x']
            sage: K.karatsuba_threshold()
            8
            sage: K = QQ['x']['y']
            sage: K.karatsuba_threshold()
            0
        """
        return self._Karatsuba_threshold

    def set_karatsuba_threshold(self, Karatsuba_threshold):
        """
        Changes the default threshold for this ring in the method
        :meth:`_mul_karatsuba` to fall back to the schoolbook algorithm.

        .. warning::

           This method may have a negative performance impact in polynomial
           arithmetic. So use it at your own risk.

        EXAMPLES::

            sage: K = QQ['x']
            sage: K.karatsuba_threshold()
            8
            sage: K.set_karatsuba_threshold(0)
            sage: K.karatsuba_threshold()
            0
        """
        self._Karatsuba_threshold = int(Karatsuba_threshold)

    def polynomials( self, of_degree=None, max_degree=None ):
        """
        Return an iterator over the polynomials of specified degree.

        INPUT: Pass exactly one of:

        - ``max_degree`` -- an int; the iterator will generate all polynomials
          which have degree less than or equal to ``max_degree``

        - ``of_degree`` -- an int; the iterator will generate
          all polynomials which have degree ``of_degree``

        OUTPUT: an iterator

        EXAMPLES::

            sage: P = PolynomialRing(GF(3), 'y')
            sage: for p in P.polynomials(of_degree=2): print(p)
            y^2
            y^2 + 1
            y^2 + 2
            y^2 + y
            y^2 + y + 1
            y^2 + y + 2
            y^2 + 2*y
            y^2 + 2*y + 1
            y^2 + 2*y + 2
            2*y^2
            2*y^2 + 1
            2*y^2 + 2
            2*y^2 + y
            2*y^2 + y + 1
            2*y^2 + y + 2
            2*y^2 + 2*y
            2*y^2 + 2*y + 1
            2*y^2 + 2*y + 2
            sage: for p in P.polynomials(max_degree=1): print(p)
            0
            1
            2
            y
            y + 1
            y + 2
            2*y
            2*y + 1
            2*y + 2
            sage: for p in P.polynomials(max_degree=1, of_degree=3): print(p)
            Traceback (most recent call last):
            ...
            ValueError: you should pass exactly one of of_degree and max_degree

        AUTHORS:

        - Joel B. Mohler
        """

        if self.base_ring().order() is sage.rings.infinity.infinity:
            raise NotImplementedError
        if of_degree is not None and max_degree is None:
            return self._polys_degree( of_degree )
        if max_degree is not None and of_degree is None:
            return self._polys_max( max_degree )
        raise ValueError("you should pass exactly one of of_degree and max_degree")

    def monics( self, of_degree=None, max_degree=None ):
        """
        Return an iterator over the monic polynomials of specified degree.

        INPUT: Pass exactly one of:


        - ``max_degree`` -- an int; the iterator will generate all monic
          polynomials which have degree less than or equal to ``max_degree``

        - ``of_degree`` -- an int; the iterator will generate
          all monic polynomials which have degree ``of_degree``

        OUTPUT: an iterator

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: P = PolynomialRing(GF(4, 'a'), 'y')
            sage: for p in P.monics(of_degree=2): print(p)
            y^2
            y^2 + a
            y^2 + a + 1
            y^2 + 1
            y^2 + a*y
            y^2 + a*y + a
            y^2 + a*y + a + 1
            y^2 + a*y + 1
            y^2 + (a + 1)*y
            y^2 + (a + 1)*y + a
            y^2 + (a + 1)*y + a + 1
            y^2 + (a + 1)*y + 1
            y^2 + y
            y^2 + y + a
            y^2 + y + a + 1
            y^2 + y + 1
            sage: for p in P.monics(max_degree=1): print(p)
            1
            y
            y + a
            y + a + 1
            y + 1
            sage: for p in P.monics(max_degree=1, of_degree=3): print(p)
            Traceback (most recent call last):
            ...
            ValueError: you should pass exactly one of of_degree and max_degree

        AUTHORS:

        - Joel B. Mohler
        """

        if self.base_ring().order() is sage.rings.infinity.infinity:
            raise NotImplementedError
        if of_degree is not None and max_degree is None:
            return self._monics_degree( of_degree )
        if max_degree is not None and of_degree is None:
            return self._monics_max( max_degree )
        raise ValueError("you should pass exactly one of of_degree and max_degree")


class PolynomialRing_commutative(PolynomialRing_general):
    """
    Univariate polynomial ring over a commutative ring.
    """
    def __init__(self, base_ring, name=None, sparse=False, implementation=None,
                 element_class=None, category=None):
        if base_ring not in _CommutativeRings:
            raise TypeError("Base ring %s must be a commutative ring." % repr(base_ring))
        # We trust that, if a category is given, that it is useful.
        if base_ring.is_zero():
            category = categories.algebras.Algebras(base_ring.category()).Commutative().Finite()
        else:
            defaultcat = polynomial_default_category(base_ring.category(), 1)
            category = check_default_category(defaultcat, category)
        PolynomialRing_general.__init__(self, base_ring, name=name,
                                        sparse=sparse, implementation=implementation,
                                        element_class=element_class, category=category)

    def quotient_by_principal_ideal(self, f, names=None, **kwds):
        """
        Return the quotient of this polynomial ring by the principal
        ideal (generated by) `f`.

        INPUT:

        - ``f`` -- either a polynomial in ``self``, or a principal
          ideal of ``self``
        - further named arguments that are passed to the quotient constructor

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: I = (x^2 - 1) * R
            sage: R.quotient_by_principal_ideal(I)                                      # needs sage.libs.pari
            Univariate Quotient Polynomial Ring in xbar
             over Rational Field with modulus x^2 - 1

        The same example, using the polynomial instead of the ideal,
        and customizing the variable name::

            sage: R.<x> = QQ[]
            sage: R.quotient_by_principal_ideal(x^2 - 1, names=('foo',))                # needs sage.libs.pari
            Univariate Quotient Polynomial Ring in foo
             over Rational Field with modulus x^2 - 1

        TESTS:

        Quotienting by the zero ideal returns ``self`` (:issue:`5978`)::

            sage: R = QQ['x']
            sage: R.quotient_by_principal_ideal(R.zero_ideal()) is R
            True
            sage: R.quotient_by_principal_ideal(0) is R
            True
        """
        from sage.rings.ideal import Ideal
        I = Ideal(f)
        if I.is_zero():
            return self
        f = I.gen()
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing
        return PolynomialQuotientRing(self, f, names, **kwds)

    def weyl_algebra(self):
        """
        Return the Weyl algebra generated from ``self``.

        EXAMPLES::

            sage: R = QQ['x']
            sage: W = R.weyl_algebra(); W                                               # needs sage.modules
            Differential Weyl algebra of polynomials in x over Rational Field
            sage: W.polynomial_ring() == R                                              # needs sage.modules
            True
        """
        from sage.algebras.weyl_algebra import DifferentialWeylAlgebra
        return DifferentialWeylAlgebra(self)

    def _roots_univariate_polynomial(self, p, ring=None, multiplicities=True, algorithm=None, degree_bound=None):
        """
        Return the list of roots of ``p``.

        INPUT:

        - ``p`` -- the polynomial whose roots are computed
        - ``ring`` -- the ring to find roots (default: the base ring of ``p``)
        - ``multiplicities`` -- boolean (default: ``True``); if ``True``,
          return a list of pairs ``(root, multiplicity)``; if ``False`` return
          a list of roots
        - ``algorithm`` -- ignored (TODO: remove)
        - ``degree_bound`` -- if not ``None``, return only roots of degree at
          most ``degree_bound``

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<x> = QQ[]
            sage: S.<y> = R[]
            sage: p = y^3 + (-x^2 - 3)*y^2 + (2*x^3 - x^2 + 3)*y - x^4 + 2*x^2 - 1
            sage: p.roots()
            [(x^2 - 2*x + 1, 1), (x + 1, 2)]
            sage: p.roots(multiplicities=False)
            [x^2 - 2*x + 1, x + 1]
            sage: p.roots(degree_bound=1)
            [(x + 1, 2)]

        TESTS:

        Check that :issue:`23639` is fixed::

            sage: foo = QQ['x']['y'].one()
            sage: foo.roots(QQ)
            []
        """
        if ring is not None and ring is not self:
            p = p.change_ring(ring)
            if degree_bound is None:
                return p.roots(multiplicities=multiplicities, algorithm=algorithm)
            return p.roots(multiplicities=multiplicities, algorithm=algorithm, degree_bound=degree_bound)

        roots = p._roots_from_factorization(p.factor(), multiplicities)
        if degree_bound is not None:
            if multiplicities:
                roots = [(r,m) for (r,m) in roots if r.degree() <= degree_bound]
            else:
                roots = [r for r in roots if r.degree() <= degree_bound]
        return roots


class PolynomialRing_integral_domain(PolynomialRing_commutative, PolynomialRing_singular_repr,
                                     IntegralDomain):
    def __init__(self, base_ring, name='x', sparse=False, implementation=None,
            element_class=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x'); R
            Univariate Polynomial Ring in x over Integer Ring
            sage: type(R.gen())                                                         # needs sage.libs.flint
            <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>

            sage: R = PRing(ZZ, 'x', implementation='NTL'); R                           # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
            sage: type(R.gen())                                                         # needs sage.libs.ntl
            <class 'sage.rings.polynomial.polynomial_integer_dense_ntl.Polynomial_integer_dense_ntl'>
        """
        self._implementation_repr = ''
        if element_class is None:
            given_implementation = implementation
            for implementation in self._implementation_names(implementation, base_ring, sparse):
                if base_ring is ZZ:
                    if implementation == 'NTL':
                        try:
                            from sage.rings.polynomial.polynomial_integer_dense_ntl \
                                import Polynomial_integer_dense_ntl as element_class
                        except ImportError:
                            if given_implementation:
                                raise
                            continue
                        self._implementation_repr = ' (using NTL)'
                    elif implementation == 'FLINT':
                        try:
                            from .polynomial_integer_dense_flint \
                                import Polynomial_integer_dense_flint as element_class
                        except ImportError:
                            if given_implementation:
                                raise
                            continue
                break
        PolynomialRing_commutative.__init__(self, base_ring, name=name,
                                            sparse=sparse, implementation=implementation,
                                            element_class=element_class, category=category)
        self._has_singular = can_convert_to_singular(self)

    @cached_method(key=lambda self, d, q, sign, lead: (d, q, sign, tuple([x if isinstance(x, (tuple, list)) else (x, 0) for x in lead]) if isinstance(lead, (tuple, list)) else ((lead, 0))))
    def weil_polynomials(self, d, q, sign=1, lead=1):
        r"""
        Return all integer polynomials whose complex roots all have a specified
        absolute value.

        Such polynomials `f` satisfy a functional equation

        .. MATH::

            T^d f(q/T) = s q^{d/2} f(T)

        where `d` is the degree of `f`, `s` is a sign and `q^{1/2}` is the
        absolute value of the roots of `f`.

        INPUT:

        - ``d`` -- integer; the degree of the polynomials

        - ``q`` -- integer; the square of the complex absolute value of the
          roots

        - ``sign`` -- integer (default: `1`); the sign `s` of the functional
          equation

        - ``lead`` -- integer; list of integers or list of pairs of integers
          (default: `1`), constraints on the leading few coefficients of the
          generated polynomials. If pairs `(a, b)` of integers are given, they
          are treated as a constraint of the form `\equiv a \pmod{b}`; the
          moduli must be in decreasing order by divisibility, and the modulus
          of the leading coefficient must be 0.

        .. SEEALSO::

            More documentation and additional options are available using the
            iterator
            :class:`sage.rings.polynomial.weil.weil_polynomials.WeilPolynomials`
            directly. In addition, polynomials have a method
            :meth:`is_weil_polynomial` to test whether or not the given
            polynomial is a Weil polynomial.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: R.<T> = ZZ[]
            sage: L = R.weil_polynomials(4, 2)
            sage: len(L)
            35
            sage: L[9]
            T^4 + T^3 + 2*T^2 + 2*T + 4
            sage: all(p.is_weil_polynomial() for p in L)
            True

        Setting multiple leading coefficients::

            sage: R.<T> = QQ[]
            sage: l = R.weil_polynomials(4, 2, lead=((1,0), (2,4), (1,2))); l           # needs sage.libs.flint
            [T^4 + 2*T^3 + 5*T^2 + 4*T + 4,
             T^4 + 2*T^3 + 3*T^2 + 4*T + 4,
             T^4 - 2*T^3 + 5*T^2 - 4*T + 4,
             T^4 - 2*T^3 + 3*T^2 - 4*T + 4]

        We do not require Weil polynomials to be monic. This example generates Weil
        polynomials associated to K3 surfaces over `\GF{2}` of Picard number at least 12::

            sage: R.<T> = QQ[]
            sage: l = R.weil_polynomials(10, 1, lead=2)                                 # needs sage.libs.flint
            sage: len(l)                                                                # needs sage.libs.flint
            4865
            sage: l[len(l)//2]                                                          # needs sage.libs.flint
            2*T^10 + T^8 + T^6 + T^4 + T^2 + 2

        TESTS:

        We check that products of Weil polynomials are also listed as Weil
        polynomials::

            sage: all((f * g) in R.weil_polynomials(6, q) for q in [3, 4]                                               # needs sage.libs.flint
            ....:     for f in R.weil_polynomials(2, q) for g in R.weil_polynomials(4, q))
            True

        We check that irreducible Weil polynomials of degree 6 are CM::

            sage: simples = [f for f in R.weil_polynomials(6, 3) if f.is_irreducible()]                                 # needs sage.libs.flint
            sage: len(simples)                                                                                          # needs sage.libs.flint
            348
            sage: reals = [R([f[3+i] + sum((-3)^j * (i+2*j)/(i+j) * binomial(i+j,j) * f[3+i+2*j]                        # needs sage.libs.flint
            ....:                          for j in range(1, (3+i)//2 + 1))
            ....:          for i in range(4)]) for f in simples]

        Check that every polynomial in this list has 3 real roots between `-2
        \sqrt{3}` and `2 \sqrt{3}`::

            sage: roots = [f.roots(RR, multiplicities=False) for f in reals]                                            # needs sage.libs.flint
            sage: all(len(L) == 3 and all(x^2 <= 12 for x in L) for L in roots)                                         # needs sage.libs.flint
            True

        Finally, check that the original polynomials are reconstructed as CM
        polynomials::

            sage: all(f == T^3*r(T + 3/T) for (f, r) in zip(simples, reals))                                            # needs sage.libs.flint
            True

        A simple check (not sufficient)::

            sage: all(f.number_of_real_roots() == 0 for f in simples)                                                   # needs sage.libs.flint
            True
        """
        R = self.base_ring()
        if not (R is ZZ or R is QQ):
            raise ValueError("Weil polynomials have integer coefficients")
        from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
        return list(WeilPolynomials(d, q, sign, lead, polring=self))

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain
            sage: PolynomialRing_integral_domain._implementation_names_impl(None, ZZ, False)
            ['FLINT', None]
            sage: PolynomialRing_integral_domain._implementation_names_impl(None, ZZ, True)
            [None, 'generic']
            sage: PolynomialRing_integral_domain._implementation_names_impl(None, QQ, False)
            [None, 'generic']
            sage: PolynomialRing_integral_domain._implementation_names_impl(None, QQ, True)
            [None, 'generic']
        """
        if base_ring is ZZ and not sparse:
            defaults = ["FLINT", None]
            if implementation in defaults:
                return defaults
            elif implementation in ["NTL", "generic"]:
                return [implementation]
        elif implementation is None or implementation == "generic":
            return [None, "generic"]
        return NotImplemented

    def _repr_(self):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x', implementation='NTL'); R                           # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
        """
        s = PolynomialRing_commutative._repr_(self)
        return s + self._implementation_repr

    def construction(self):
        """
        Return the construction functor.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x'); R
            Univariate Polynomial Ring in x over Integer Ring
            sage: functor, arg = R.construction(); functor, arg
            (Poly[x], Integer Ring)
            sage: functor.implementation is None
            True

            sage: # needs sage.libs.ntl
            sage: R = PRing(ZZ, 'x', implementation='NTL'); R
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
            sage: functor, arg = R.construction(); functor, arg
            (Poly[x], Integer Ring)
            sage: functor.implementation
            'NTL'
        """
        implementation = None
        # NOTE: This is obviously not a complete solution. The parents
        # don't keep track in a clean way what the implementation is.
        # Issue #31852 is the task of finding a general solution for
        # construction functors of parents with multiple
        # implementations, such as MatrixSpace, Polyhedron, and
        # PolynomialRing.
        if 'NTL' in self._implementation_repr:
            implementation = 'NTL'
        return categories.pushout.PolynomialFunctor(self.variable_name(), sparse=self.is_sparse(),
                                                    implementation=implementation), self.base_ring()


class PolynomialRing_field(PolynomialRing_integral_domain):
    def __init__(self, base_ring, name='x', sparse=False, implementation=None,
                 element_class=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_field as PRing
            sage: R = PRing(QQ, 'x'); R
            Univariate Polynomial Ring in x over Rational Field
            sage: type(R.gen())                                                         # needs sage.libs.flint
            <class 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>
            sage: R = PRing(QQ, 'x', sparse=True); R
            Sparse Univariate Polynomial Ring in x over Rational Field
            sage: type(R.gen())
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_field_with_category.element_class'>

            sage: # needs sage.rings.real_mpfr
            sage: R = PRing(CC, 'x'); R
            Univariate Polynomial Ring in x over Complex Field with 53 bits of precision
            sage: type(R.gen())
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_field_with_category.element_class'>

        Demonstrate that :issue:`8762` is fixed::

            sage: R.<x> = PolynomialRing(GF(next_prime(10^20)), sparse=True)            # needs sage.rings.finite_rings
            sage: x^(10^20)  # this should be fast                                      # needs sage.rings.finite_rings
            x^100000000000000000000
        """
        def _element_class():
            if element_class:
                return element_class
            if sparse:
                from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_sparse_field
                return Polynomial_generic_sparse_field
            if isinstance(base_ring, rational_field.RationalField):
                try:
                    from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
                    return Polynomial_rational_flint
                except ImportError:
                    pass
            elif isinstance(base_ring, NumberField):
                if base_ring.is_absolute():
                    from sage.rings.polynomial.polynomial_number_field import Polynomial_absolute_number_field_dense
                    return Polynomial_absolute_number_field_dense
                else:
                    from sage.rings.polynomial.polynomial_number_field import Polynomial_relative_number_field_dense
                    return Polynomial_relative_number_field_dense
            elif isinstance(base_ring, sage.rings.abc.RealField):
                try:
                    from .polynomial_real_mpfr_dense import PolynomialRealDense
                    return PolynomialRealDense
                except ImportError:
                    pass
            elif isinstance(base_ring, sage.rings.abc.ComplexBallField):
                try:
                    from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
                    return Polynomial_complex_arb
                except ImportError:
                    pass
            from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense_field
            return Polynomial_generic_dense_field

        if category is None:
            cat = PrincipalIdealDomains()
        else:
            cat = category & PrincipalIdealDomains()

        PolynomialRing_integral_domain.__init__(self, base_ring, name=name,
                                                sparse=sparse, implementation=implementation,
                                                element_class=_element_class(), category=cat)

    def _ideal_class_(self, n=0):
        """
        Return the class representing ideals in univariate polynomial rings
        over fields.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: R._ideal_class_()
            <class 'sage.rings.polynomial.ideal.Ideal_1poly_field'>
        """
        from sage.rings.polynomial.ideal import Ideal_1poly_field
        return Ideal_1poly_field

    def divided_difference(self, points, full_table=False):
        r"""
        Return the Newton divided-difference coefficients of the
        Lagrange interpolation polynomial through ``points``.

        INPUT:

        - ``points`` -- list of pairs `(x_0, y_0), (x_1, y_1),
          \dots, (x_n, y_n)` of elements of the base ring of ``self``,
          where `x_i - x_j` is invertible for `i \neq j`.  This method
          converts the `x_i` and `y_i` into the base ring of ``self``.

        - ``full_table`` -- boolean (default: ``False``); if ``True``,
          return the full divided-difference table.  If ``False``,
          only return entries along the main diagonal; these are the
          Newton divided-difference coefficients `F_{i,i}`.

        OUTPUT:

        The Newton divided-difference coefficients of the `n`-th
        Lagrange interpolation polynomial `P_n(x)` that passes through
        the points in ``points`` (see :meth:`lagrange_polynomial`).
        These are the coefficients `F_{0,0}, F_{1,1}, \dots, F_{n,n}`
        in the base ring of ``self`` such that

        .. MATH::

            P_n(x) = \sum_{i=0}^n F_{i,i} \prod_{j=0}^{i-1} (x - x_j)

        EXAMPLES:

        Only return the divided-difference coefficients `F_{i,i}`.
        This example is taken from Example 1, page 121 of [BF2005]_::

            sage: # needs sage.rings.real_mpfr
            sage: points = [(1.0, 0.7651977), (1.3, 0.6200860), (1.6, 0.4554022),
            ....:           (1.9, 0.2818186), (2.2, 0.1103623)]
            sage: R = PolynomialRing(RR, "x")
            sage: R.divided_difference(points)
            [0.765197700000000,
             -0.483705666666666,
             -0.108733888888889,
             0.0658783950617283,
             0.00182510288066044]

        Now return the full divided-difference table::

            sage: # needs sage.rings.real_mpfr
            sage: points = [(1.0, 0.7651977), (1.3, 0.6200860), (1.6, 0.4554022),
            ....:           (1.9, 0.2818186), (2.2, 0.1103623)]
            sage: R = PolynomialRing(RR, "x")
            sage: R.divided_difference(points, full_table=True)
            [[0.765197700000000],
             [0.620086000000000, -0.483705666666666],
             [0.455402200000000, -0.548946000000000, -0.108733888888889],
             [0.281818600000000, -0.578612000000000,
                                -0.0494433333333339, 0.0658783950617283],
             [0.110362300000000, -0.571520999999999, 0.0118183333333349,
                                0.0680685185185209, 0.00182510288066044]]

        The following example is taken from Example 4.12, page 225 of
        [MF1999]_::

            sage: points = [(1, -3), (2, 0), (3, 15), (4, 48), (5, 105), (6, 192)]
            sage: R = PolynomialRing(QQ, "x")
            sage: R.divided_difference(points)
            [-3, 3, 6, 1, 0, 0]
            sage: R.divided_difference(points, full_table=True)
            [[-3],
             [0, 3],
             [15, 15, 6],
             [48, 33, 9, 1],
             [105, 57, 12, 1, 0],
             [192, 87, 15, 1, 0, 0]]
        """
        to_base_ring = self.base_ring()
        points = [tuple(to_base_ring(c) for c in p) for p in points]
        n = len(points)
        F = [[points[i][1]] for i in range(n)]
        for i in range(1, n):
            for j in range(1, i+1):
                numer = F[i][j-1] - F[i-1][j-1]
                denom = points[i][0] - points[i-j][0]
                F[i].append(numer / denom)
        if full_table:
            return F
        else:
            return [F[i][i] for i in range(n)]

    def lagrange_polynomial(self, points, algorithm='divided_difference', previous_row=None):
        r"""
        Return the Lagrange interpolation polynomial through the
        given points.

        INPUT:

        - ``points`` -- list of pairs `(x_0, y_0), (x_1, y_1),
          \dots, (x_n, y_n)` of elements of the base ring of ``self``,
          where `x_i - x_j` is invertible for `i \neq j`.  This method
          converts the `x_i` and `y_i` into the base ring of ``self``.

        - ``algorithm`` -- (default: ``'divided_difference'``) one of
          the following:

          - ``'divided_difference'``: use the method of divided
            differences.

          - ``'neville'``: adapt Neville's method as
            described on page 144 of [BF2005]_ to recursively generate
            the Lagrange interpolation polynomial.  Neville's method
            generates a table of approximating polynomials, where the
            last row of that table contains the `n`-th Lagrange
            interpolation polynomial.  The adaptation implemented by
            this method is to only generate the last row of this
            table, instead of the full table itself.  Generating the
            full table can be memory inefficient.

        - ``previous_row`` -- (default: ``None``) this option is only
          relevant if used with ``algorithm='neville'``.  If provided,
          this should be the last row of the table resulting from a
          previous use of Neville's method.  If such a row is passed,
          then ``points`` should consist of both previous and new
          interpolating points.  Neville's method will then use that
          last row and the interpolating points to generate a new row
          containing an interpolation polynomial for the new points.

        OUTPUT:

        The Lagrange interpolation polynomial through the points
        `(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)`.  This is the
        unique polynomial `P_n` of degree at most `n` in ``self``
        satisfying `P_n(x_i) = y_i` for `0 \le i \le n`.

        EXAMPLES:

        By default, we use the method of divided differences::

            sage: R = PolynomialRing(QQ, 'x')
            sage: f = R.lagrange_polynomial([(0,1), (2,2), (3,-2), (-4,9)]); f
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1
            sage: f(0)
            1
            sage: f(2)
            2
            sage: f(3)
            -2
            sage: f(-4)
            9

            sage: # needs sage.rings.finite_rings
            sage: R = PolynomialRing(GF(2**3, 'a'), 'x')
            sage: a = R.base_ring().gen()
            sage: f = R.lagrange_polynomial([(a^2+a, a), (a, 1), (a^2, a^2+a+1)]); f
            a^2*x^2 + a^2*x + a^2
            sage: f(a^2 + a)
            a
            sage: f(a)
            1
            sage: f(a^2)
            a^2 + a + 1

        Now use a memory efficient version of Neville's method::

            sage: R = PolynomialRing(QQ, 'x')
            sage: R.lagrange_polynomial([(0,1), (2,2), (3,-2), (-4,9)],
            ....:                       algorithm='neville')
            [9,
            -11/7*x + 19/7,
            -17/42*x^2 - 83/42*x + 53/7,
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1]

            sage: # needs sage.rings.finite_rings
            sage: R = PolynomialRing(GF(2**3, 'a'), 'x')
            sage: a = R.base_ring().gen()
            sage: R.lagrange_polynomial([(a^2+a, a), (a, 1), (a^2, a^2+a+1)],
            ....:                       algorithm='neville')
            [a^2 + a + 1, x + a + 1, a^2*x^2 + a^2*x + a^2]

        Repeated use of Neville's method to get better Lagrange
        interpolation polynomials::

            sage: R = PolynomialRing(QQ, 'x')
            sage: p = R.lagrange_polynomial([(0,1), (2,2)], algorithm='neville')
            sage: R.lagrange_polynomial([(0,1), (2,2), (3,-2), (-4,9)],
            ....:                       algorithm='neville', previous_row=p)[-1]
            -23/84*x^3 - 11/84*x^2 + 13/7*x + 1

            sage: # needs sage.rings.finite_rings
            sage: R = PolynomialRing(GF(2**3, 'a'), 'x')
            sage: a = R.base_ring().gen()
            sage: p = R.lagrange_polynomial([(a^2+a, a), (a, 1)], algorithm='neville')
            sage: R.lagrange_polynomial([(a^2+a, a), (a, 1), (a^2, a^2+a+1)],
            ....:                       algorithm='neville', previous_row=p)[-1]
            a^2*x^2 + a^2*x + a^2

        TESTS:

        The value for ``algorithm`` must be either
        ``'divided_difference'`` (default), or ``'neville'``::

            sage: R = PolynomialRing(QQ, 'x')
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm='abc')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'divided_difference' or 'neville'
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm='divided difference')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'divided_difference' or 'neville'
            sage: R.lagrange_polynomial([(0,1),(2,2),(3,-2),(-4,9)], algorithm='')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'divided_difference' or 'neville'

        Make sure that :issue:`10304` is fixed.  The return value
        should always be an element of ``self`` in the case of
        ``divided_difference``, or a list of elements of ``self`` in
        the case of ``neville``::

            sage: R = PolynomialRing(QQ, 'x')
            sage: R.lagrange_polynomial([]).parent() == R
            True
            sage: R.lagrange_polynomial([(2, 3)]).parent() == R
            True
            sage: row = R.lagrange_polynomial([], algorithm='neville')
            sage: all(poly.parent() == R for poly in row)
            True
            sage: row = R.lagrange_polynomial([(2, 3)], algorithm='neville')
            sage: all(poly.parent() == R for poly in row)
            True

        Check that base fields of positive characteristic are treated
        correctly (see :issue:`9787`)::

            sage: R.<x> = GF(101)[]
            sage: R.lagrange_polynomial([[1, 0], [2, 0]])
            0
            sage: R.lagrange_polynomial([[1, 0], [2, 0], [3, 0]])
            0
        """
        # Perhaps we should be slightly stricter on the input and use
        # self.base_ring().coerce here and in the divided_difference()
        # method above.  However, this breaks an example in
        # sage.tests.french_book.nonlinear_doctest where the base ring
        # is CC, but the function values lie in the symbolic ring.
        to_base_ring = self.base_ring()
        points = [[to_base_ring(u) for u in x] for x in points]
        var = self.gen()

        # use the method of divided-difference
        if algorithm == "divided_difference":
            # Evaluate in nested form, similar to Horner's method. This is
            # more efficient than evaluation using the definition of
            # Lagrange interpolation polynomial by means of divided
            # difference.
            n = len(points)
            if n == 0:
                return self.zero()

            F = self.divided_difference(points)
            P = self.coerce(F[n-1])
            for i in range(n-2, -1, -1):
                P *= (var - points[i][0])
                P += F[i]
            return P

            # Evaluate using the definition of Lagrange interpolation
            # polynomial by means of divided difference. This is slow
            # compared to that above, which is in nested form.
#             P = 0
#             for i in range(n):
#                 prod = 1
#                 for j in range(i):
#                     prod *= (var - points[j][0])
#                 P += (F[i] * prod)
#             return P

        # using Neville's method for recursively generating the
        # Lagrange interpolation polynomial
        elif algorithm == "neville":
            if previous_row is None:
                previous_row = []
            N = len(points)
            M = len(previous_row)
            # During the computation, P keeps track of the previous row,
            # and Q keeps track of the current row
            P = previous_row + [None] * (N - M) # use results of previous computation if available
            Q = [None] * N
            for i in range(M, N):
                Q[0] = self.coerce(points[i][1])  # start populating the current row
                for j in range(1, 1 + i):
                    numer = (var - points[i - j][0]) * Q[j - 1] - (var - points[i][0]) * P[j - 1]
                    denom = points[i][0] - points[i - j][0]
                    Q[j] = numer / denom
                P, Q = Q, P # the current row is complete, reuse the old P to hold the next row
            return P # return the last row in the Neville table

#        # use the definition of Lagrange interpolation polynomial
#        elif algorithm == "definition":
#            def Pj(j):
#                denom = 1
#                divis = 1
#                for i in range(len(points)):
#                    if i!=j:
#                        denom *= (var          - points[i][0])
#                        divis *= (points[j][0] - points[i][0])
#            return denom/divis
#
#            P = 0
#            for j in range(len(points)):
#                P += Pj(j)*points[j][1]
#            return P

        else:
            raise ValueError("algorithm must be one of 'divided_difference' or 'neville'")

    @cached_method
    def fraction_field(self):
        """
        Return the fraction field of ``self``.

        EXAMPLES::

            sage: QQbar['x'].fraction_field()
            Fraction Field of Univariate Polynomial Ring in x over Algebraic
            Field

        TESTS:

        Check that :issue:`25449` has been resolved::

            sage: # needs sage.rings.finite_rings
            sage: k = GF(25453)
            sage: F.<x> = FunctionField(k)
            sage: R.<t> = k[]
            sage: t(x)
            x

            sage: # needs sage.rings.finite_rings
            sage: k = GF(55667)
            sage: F.<x> = FunctionField(k)
            sage: R.<t> = k[]
            sage: t(x)
            x

        Fixed :issue:`37374`::

            sage: x = PolynomialRing(GF(37), ['x'], sparse=True).fraction_field().gen()
            sage: type(x.numerator())
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_field_with_category.element_class'>
            sage: (x^8 + 16*x^6 + 4*x^4 + x^2 + 12).numerator() - 1
            x^8 + 16*x^6 + 4*x^4 + x^2 + 11
        """
        from sage.rings.fraction_field import FractionField_1poly_field
        return FractionField_1poly_field(self)


class PolynomialRing_dense_finite_field(PolynomialRing_field):
    """
    Univariate polynomial ring over a finite field.

    EXAMPLES::

        sage: R = PolynomialRing(GF(27, 'a'), 'x')                                      # needs sage.rings.finite_rings
        sage: type(R)                                                                   # needs sage.rings.finite_rings
        <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_finite_field_with_category'>
    """
    def __init__(self, base_ring, name='x', element_class=None, implementation=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_finite_field
            sage: R = PolynomialRing_dense_finite_field(GF(5), implementation='generic')
            sage: type(R(0))
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_finite_field_with_category.element_class'>

            sage: S = PolynomialRing_dense_finite_field(GF(25, 'a'), implementation='NTL')          # needs sage.libs.ntl sage.rings.finite_rings
            sage: type(S(0))                                                                        # needs sage.libs.ntl sage.rings.finite_rings
            <class 'sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX'>

            sage: S = PolynomialRing_dense_finite_field(GF(64), implementation='superfast')         # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: unknown implementation 'superfast' for dense polynomial rings over Finite Field in z6 of size 2^6
        """
        if element_class is None:
            given_implementation = implementation
            for implementation in self._implementation_names(implementation, base_ring):
                if implementation == "NTL":
                    try:
                        from sage.libs.ntl.ntl_ZZ_pEContext import ntl_ZZ_pEContext
                        from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX
                        from sage.rings.polynomial.polynomial_zz_pex import Polynomial_ZZ_pEX
                    except ImportError:
                        if given_implementation:
                            raise
                        continue
                    p = base_ring.characteristic()
                    self._modulus = ntl_ZZ_pEContext(ntl_ZZ_pX(list(base_ring.modulus()), p))
                    element_class = Polynomial_ZZ_pEX
                break
        PolynomialRing_field.__init__(self, base_ring, sparse=False, name=name,
                                      implementation=implementation, element_class=element_class)

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_finite_field
            sage: PolynomialRing_dense_finite_field._implementation_names_impl("NTL", GF(4), False)
            ['NTL', None]
            sage: PolynomialRing_dense_finite_field._implementation_names_impl(None, GF(4), False)
            ['NTL', None]
            sage: PolynomialRing_dense_finite_field._implementation_names_impl("generic", GF(4), False)
            ['generic']
            sage: PolynomialRing_dense_finite_field._implementation_names_impl("FLINT", GF(4), False)
            NotImplemented
            sage: PolynomialRing_dense_finite_field._implementation_names_impl(None, GF(4), True)
            NotImplemented
        """
        if sparse:
            return NotImplemented
        defaults = ["NTL", None]
        if implementation in defaults:
            return defaults
        elif implementation == "generic":
            return [implementation]
        return NotImplemented

    def irreducible_element(self, n, algorithm=None):
        """
        Construct a monic irreducible polynomial of degree `n`.

        INPUT:

        - ``n`` -- integer; degree of the polynomial to construct

        - ``algorithm`` -- string (algorithm to use) or ``None``:

          - ``'random'`` or ``None``:
            try random polynomials until an irreducible one is found

          - ``'uniform_random'``:
            try random polynomials sampled from all monic n-degree
            polynomials until an irreducible one is found.

          - ``'first_lexicographic'``: try polynomials in
            lexicographic order until an irreducible one is found

        OUTPUT: A monic irreducible polynomial of degree `n` in the given finite ring ``self``.

        EXAMPLES::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: f = GF(5^3, 'a')['x'].irreducible_element(2)
            sage: f.degree()
            2
            sage: f.is_irreducible()
            True
            sage: R = GF(19)['x']
            sage: R.irreducible_element(21, algorithm='first_lexicographic')
            x^21 + x + 5
            sage: R = GF(5**2, 'a')['x']
            sage: R.irreducible_element(17, algorithm='first_lexicographic')
            x^17 + a*x + 4*a + 3

        TESTS::

            sage: import collections
            sage: from sage.stats.statistical_tests import sigmas_from_uniform
            sage: R = GF(2**3, 'a')['x']; samples = 10000

        This test verifies that the "uniform_random" generates expected distributions::

            sage: sample_counts = collections.Counter(R.irreducible_element(3, algorithm="uniform_random") for _ in range(samples))
            sage: abs(sigmas_from_uniform(sample_counts.values())) < 5 # physics thinks 5 sigma is worth investigating.
            True

        This test verifies that the "random" generates extremely unlikely distributions::

            sage: sample_counts = collections.Counter(R.irreducible_element(3, algorithm="random") for _ in range(samples))
            sage: abs(sigmas_from_uniform(sample_counts.values())) > 800 # More than 800 sigmas is NOT random
            True

        .. TODO::

            Deprecate the ``random`` approach that isn't uniformly random over the irreducible elements.
            
        AUTHORS:

        - Peter Bruin (June 2013)
        - Jean-Pierre Flori (May 2014)
        """
        if n < 1:
            raise ValueError("degree must be at least 1")

        if algorithm is None:
            algorithm = "random"

        if algorithm == "random":
            while True:
                f = self.gen()**n + self.random_element(degree=(0, n - 1))
                if f.is_irreducible():
                    return f
        elif algorithm == "first_lexicographic":
            for g in self.polynomials(max_degree=n-1):
                f = self.gen()**n + g
                if f.is_irreducible():
                    return f 
        elif algorithm == "uniform_random":
            R = self.base_ring()
            while True:
                g = [R.random_element() for i in range(n)]
                g.append(1)
                f = self(g)
                if f.is_irreducible():
                    return f
        else:
            raise ValueError("no such algorithm for finding an irreducible polynomial: %s" % algorithm)

    def _roth_ruckenstein(self, p, degree_bound, precision):
        r"""
        Return all polynomials which are a solution to the, possibly modular,
        root-finding problem.

        This is the core of Roth-Ruckenstein's algorithm where all conversions,
        checks and parent-extraction have been done.

        INPUT:

        - ``p`` -- a nonzero polynomial over ``F[x][y]``. The polynomial ``p``
          should be first truncated to ``precision``

        - ``degree_bound`` -- a bound on the degree of the roots of ``p`` that
          the algorithm computes

        - ``precision`` -- if given, roots are computed modulo `x^d` where `d` is
          ``precision`` (see below)

        OUTPUT: the list of roots of ``p`` of degree at most ``degree_bound``:

        - If `precision = None` actual roots are computed, i.e. all `f \in F[x]`
          such that `p(f) = 0`.

        - If ``precision = k`` for some integer ``k``, then all `f \in \F[x]` such
          that `Q(f) \equiv 0 \mod x^k` are computed. This set is infinite, thus it
          represented as a list of pairs in `F[x] \times \ZZ_+`, where
          `(f, d)` denotes that `Q(f + x^d h) \equiv 0 \mod x^k` for any `h \in
          F[[x]]`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: F = GF(17)
            sage: Px.<x> = F[]
            sage: Pxy.<y> = Px[]
            sage: p = (y - (x**2 + x + 1)) * (y**2 - x + 1) * (y - (x**3 + 4*x + 16))
            sage: Px._roth_ruckenstein(p, 3, None)
            [x^3 + 4*x + 16, x^2 + x + 1]
            sage: Px._roth_ruckenstein(p, 2, None)
            [x^2 + x + 1]
            sage: Px._roth_ruckenstein(p, 1, 2)
            [(4*x + 16, 2), (2*x + 13, 2), (15*x + 4, 2), (x + 1, 2)]
        """
        def roth_rec(p, lam, k, g):
            r"""
            Recursive core method for Roth-Ruckenstein algorithm.

            INPUT:

            - ``p`` -- the current value of the polynomial
            - ``lam`` -- is the power of x whose coefficient is being computed
            - ``k`` -- the remaining precision to handle (if ``precision`` is given)
            - ``g`` -- the root being computed
            """
            if precision and k <= 0:
                solutions.append((g, lam))
                return
            val = min(c.valuation() for c in p)
            if precision:
                k = k - val
            T = p.map_coefficients(lambda c:c.shift(-val))
            Ty = T.map_coefficients(lambda c:c[0]).change_ring(F)
            if Ty.is_zero() or (precision and k <= 0):
                if precision:
                    solutions.append((g, lam))
                else:
                    solutions.append(g)
                return
            roots = Ty.roots(multiplicities=False)
            for gamma in roots:
                g_new = g + gamma*x**lam
                if lam < degree_bound:
                    Tg = T(x*y + gamma)
                    roth_rec(Tg , lam+1, k, g_new)
                else:
                    if precision:
                        solutions.append((g_new, lam+1))
                    elif p(gamma).is_zero():
                        solutions.append(g_new)
            return

        x = self.gen()
        y = p.parent().gen()
        F = self.base_ring()
        solutions = []
        g = self.zero()

        roth_rec(p, 0, precision, g)
        return solutions

    def _alekhnovich(self, p, degree_bound, precision=None, dc_threshold=None):
        r"""
        Use Alekhnovich's Divide & Conquer variant of Roth-Ruckenstein's
        rootfinding algorithm to find roots modulo-up-to-some-precision of a `Q \in
        F[x][y]` where `F` is a finite field. Supports a mixed strategy with
        Roth-Ruckenstein applied at lowest precision.

        INPUT:

        - ``p`` -- a nonzero polynomial over ``F[x][y]``. The polynomial ``p``
          should be first truncated to ``precision``

        - ``degree_bound`` -- a bound on the degree of the roots of ``p`` that
          the algorithm computes

        - ``precision`` -- if given, roots are computed modulo `x^d` where `d` is
          ``precision`` (see below)

        - ``dc_threshold`` -- if given, the algorithm calls :meth:`_roth_ruckenetein`
          to compute roots of degree at most ``dc_threshold``

        OUTPUT: the list of roots of ``p`` of degree at most ``degree_bound``:

        - If `precision = None` actual roots are computed, i.e. all `f \in F[x]`
          such that `p(f) = 0`.

        - If ``precision = k`` for some integer ``k``, then all `f \in \F[x]` such
          that `Q(f) \equiv 0 \mod x^k` are computed. This set is infinite, thus it
          represented as a list of pairs in `F[x] \times \ZZ_+`, where
          `(f, d)` denotes that `Q(f + x^d h) \equiv 0 \mod x^k` for any `h \in
          F[[x]]`.

        .. NOTE::

            Non-exhaustive testing tends to indicate that ``dc_threhold = None`` is,
            surprisingly, the best strategy. (See the example section.)

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(17)[]
            sage: S.<y> = R[]
            sage: p = (y - 2*x^2 - 3*x - 14) * (y - 3*x + 2) * (y - 1)
            sage: R._alekhnovich(p, 2)
            [3*x + 15, 2*x^2 + 3*x + 14, 1]
            sage: R._alekhnovich(p, 1)
            [3*x + 15, 1]
            sage: R._alekhnovich(p, 1, precision=2)
            [(3*x + 15, 2), (3*x + 14, 2), (1, 2)]

        Example of benchmark to check that `dc_threshold = None` is better::

            sage: # not tested, needs sage.rings.finite_rings
            sage: p = prod(y - R.random_element(20)
            ....:          for _ in range(10)) * S.random_element(10,10)
            sage: %timeit _alekhnovich(R, p, 20, dc_threshold = None)
            1 loop, best of 3: 418 ms per loop
            sage: %timeit _alekhnovich(R, p, 20, dc_threshold = 1)
            1 loop, best of 3: 416 ms per loop
            sage: %timeit _alekhnovich(R, p, 20, dc_threshold = 2)
            1 loop, best of 3: 418 ms per loop
            sage: %timeit _alekhnovich(R, p, 20, dc_threshold = 3)
            1 loop, best of 3: 454 ms per loop
            sage: %timeit _alekhnovich(R, p, 20, dc_threshold = 4)
            1 loop, best of 3: 519 ms per loop

        AUTHORS:

        - Johan Rosenkilde (2015) -- Original implementation
        - Bruno Grenet (August 2016) -- Incorporation into SageMath and polishing
        """
        def alekh_rec(p, k, degree_bound, lvl):
            r"""
            Recursive core method for Alekhnovich algorithm.

            INPUT:

            - ``p`` -- the current value of the polynomial
            - ``k`` -- the number of coefficients left to be computed
            - ``degree_bound`` -- the current degree bound
            - ``lvl`` -- the level in the recursion tree
            """
            if k <= 0:
                return [ (self.zero(),0) ]
            elif degree_bound < 0:
                # The only possible root of (current) p, if any, is y = 0
                if p(0).is_zero() or p(0).valuation() >= k:
                    return [ (self.zero(),0) ]
                else:
                    return []
            elif k == 1 or degree_bound == 0:
                #Either one coefficient left to be computed, or p has only one coefficient
                py = self([c[0] for c in p.list()])  # py = p(x=0, y)
                if py.is_zero():
                    return [ (self.zero(), 0) ]
                roots = py.roots(multiplicities=False)
                return [ (self(r),1) for r in roots ]
            elif k < dc_threshold:
                # Run Roth-Ruckenstein
                return self._roth_ruckenstein(p, degree_bound=degree_bound, precision=k)
            else:
                p = p.map_coefficients(lambda c:c.truncate(k))
                half_roots = alekh_rec(p, k//2, degree_bound, lvl+1)
                whole_roots = []
                for (hi, di) in half_roots:
                    QhatT = p(hi + y*x**di)
                    if not QhatT:
                        whole_roots.append((hi,di))
                    else:
                        val = min(c.valuation() for c in QhatT)
                        Qhat = QhatT.map_coefficients(lambda c:c.shift(-val))
                        sec_half = alekh_rec(Qhat, k-val, degree_bound - di, lvl+1)
                        whole_roots.extend([ (hi + hij.shift(di), di+dij) for (hij, dij) in sec_half ])
                return whole_roots

        x = self.gen()
        y = p.parent().gen()

        # If precision is not given, find actual roots. To be sure, precision then
        # needs to be more than wdeg{1,degree_bound}(Q) since a root might have degree degree_bound.
        if precision is None:
            k = 1 + max( p[i].degree() + degree_bound*i for i in range(1+p.degree()))
        else:
            k = precision

        mod_roots = alekh_rec(p, k, degree_bound, 0)

        if precision is None:
            roots = []
            for hi,_ in mod_roots:
                if p(hi).is_zero():
                    roots.append(hi)
            return roots
        else:
            return mod_roots

    def _roots_univariate_polynomial(self, p, ring=None, multiplicities=False, algorithm=None, degree_bound=None):
        """
        Return the list of roots of ``p``.

        INPUT:

        - ``p`` -- the polynomial whose roots are computed
        - ``ring`` -- the ring to find roots (default: the base ring of ``p``)
        - ``multiplicities`` -- boolean (default: ``True``); currently, roots are only
          computed without their multiplicities
        - ``algorithm`` -- the algorithm to use: either ``'Alekhnovich'`` (default)
          or ``'Roth-Ruckenstein'``
        - ``degree_bound`` -- if not ``None``, return only roots of degree at
          most ``degree_bound``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(13)[]
            sage: S.<y> = R[]
            sage: p = y^2 + (12*x^2 + x + 11)*y + x^3 + 12*x^2 + 12*x + 1
            sage: p.roots(multiplicities=False)
            [x^2 + 11*x + 1, x + 1]
            sage: p.roots(multiplicities=False, degree_bound=1)
            [x + 1]
            sage: p.roots(multiplicities=False, algorithm='Roth-Ruckenstein')
            [x^2 + 11*x + 1, x + 1]

        TESTS:

        Check that :issue:`23639` is fixed::

            sage: R = GF(3)['x']['y']
            sage: R.one().roots(multiplicities=False)
            []
            sage: R.zero().roots(multiplicities=False)
            Traceback (most recent call last):
            ...
            ArithmeticError: roots of 0 are not defined
        """
        if multiplicities:
            raise NotImplementedError("Use multiplicities=False")

        if degree_bound is None:
            l = p.degree()
            if l < 0:
                raise ArithmeticError("roots of 0 are not defined")
            if l == 0:
                return []
            dl = p[l].degree()
            degree_bound = max((p[i].degree() - dl)//(l - i) for i in range(l) if p[i])

        if algorithm is None:
            algorithm = "Alekhnovich"

        if algorithm == "Roth-Ruckenstein":
            return self._roth_ruckenstein(p, degree_bound, None)

        elif algorithm == "Alekhnovich":
            return self._alekhnovich(p, degree_bound)

        else:
            raise ValueError("unknown algorithm '{}'".format(algorithm))


class PolynomialRing_cdvr(PolynomialRing_integral_domain):
    r"""
    A class for polynomial ring over complete discrete valuation rings
    """
    def __init__(self, base_ring, name=None, sparse=False, implementation=None,
                 element_class=None, category=None):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_cdvr

            sage: S.<x> = ZZ[]
            sage: isinstance(S, PolynomialRing_cdvr)
            False

            sage: # needs sage.rings.padics
            sage: S.<x> = Zp(5)[]
            sage: isinstance(S, PolynomialRing_cdvr)
            True
        """
        if element_class is None:
            if sparse:
                from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_sparse_cdvr
                element_class = Polynomial_generic_sparse_cdvr
            else:
                from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense_cdvr
                element_class = Polynomial_generic_dense_cdvr
        PolynomialRing_integral_domain.__init__(self, base_ring, name, sparse,
                                                implementation=implementation,
                                                element_class=element_class, category=category)


class PolynomialRing_cdvf(PolynomialRing_cdvr, PolynomialRing_field):
    """
    A class for polynomial ring over complete discrete valuation fields
    """
    def __init__(self, base_ring, name=None, sparse=False, implementation=None,
                 element_class=None, category=None):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_cdvf

            sage: S.<x> = QQ[]
            sage: isinstance(S, PolynomialRing_cdvf)
            False

            sage: S.<x> = Qp(5)[]                                                       # needs sage.rings.padics
            sage: isinstance(S, PolynomialRing_cdvf)                                    # needs sage.rings.padics
            True
        """
        if element_class is None:
            if sparse:
                from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_sparse_cdvf
                element_class = Polynomial_generic_sparse_cdvf
            else:
                from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense_cdvf
                element_class = Polynomial_generic_dense_cdvf
        PolynomialRing_field.__init__(self, base_ring, name, sparse,
                                      implementation=implementation, element_class=element_class,
                                      category=category)


class PolynomialRing_dense_padic_ring_generic(PolynomialRing_cdvr):
    r"""
    A class for dense polynomial ring over `p`-adic rings
    """
    def __init__(self, base_ring, name=None, implementation=None, element_class=None, category=None):
        PolynomialRing_cdvr.__init__(self, base_ring, sparse=False, name=name,
                                     implementation=implementation, element_class=element_class,
                                     category=category)

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        Only support ``implementation=None`` and ``sparse=False``.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_generic
            sage: PolynomialRing_dense_padic_ring_generic._implementation_names_impl(None, Zp(2), False)                # needs sage.rings.padics
            [None]
            sage: PolynomialRing_dense_padic_ring_generic._implementation_names_impl(None, Zp(2), True)                 # needs sage.rings.padics
            NotImplemented
            sage: PolynomialRing_dense_padic_ring_generic._implementation_names_impl("generic", Zp(2), False)           # needs sage.rings.padics
            NotImplemented
        """
        if implementation is None and not sparse:
            return [None]  # Not a "generic" implementation
        return NotImplemented


class PolynomialRing_dense_padic_field_generic(PolynomialRing_cdvf):
    r"""
    A class for dense polynomial ring over `p`-adic fields
    """
    def __init__(self, base_ring, name=None, implementation=None, element_class=None, category=None):
        PolynomialRing_cdvf.__init__(self, base_ring, sparse=False, name=name,
                                     implementation=implementation, element_class=element_class,
                                     category=category)

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        Only support ``implementation=None`` and ``sparse=False``.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_field_generic
            sage: PolynomialRing_dense_padic_field_generic._implementation_names_impl(None, Qp(2), False)               # needs sage.rings.padics
            [None]
            sage: PolynomialRing_dense_padic_field_generic._implementation_names_impl(None, Qp(2), True)                # needs sage.rings.padics
            NotImplemented
            sage: PolynomialRing_dense_padic_field_generic._implementation_names_impl("generic", Qp(2), False)          # needs sage.rings.padics
            NotImplemented
        """
        if implementation is None and not sparse:
            return [None]  # Not a "generic" implementation
        return NotImplemented


class PolynomialRing_dense_padic_ring_capped_relative(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, implementation=None, element_class=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_capped_relative as PRing
            sage: R = PRing(Zp(13), name='t'); R                                                                        # needs sage.rings.padics
            Univariate Polynomial Ring in t over 13-adic Ring with capped relative precision 20
            sage: type(R.gen())                                                                                         # needs sage.rings.padics
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_padic_ring_capped_relative_with_category.element_class'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.\
                    polynomial_padic_capped_relative_dense import \
                    Polynomial_padic_capped_relative_dense
            element_class = Polynomial_padic_capped_relative_dense
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring, name=name,
                                                         implementation=implementation,
                                                         element_class=element_class, category=category)


class PolynomialRing_dense_padic_ring_capped_absolute(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, implementation=None, element_class=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_capped_absolute as PRing
            sage: R = PRing(Zp(13, type='capped-abs'), name='t'); R                                                     # needs sage.rings.padics
            Univariate Polynomial Ring in t over 13-adic Ring with capped absolute precision 20
            sage: type(R.gen())                                                                                         # needs sage.rings.padics
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_padic_ring_capped_absolute_with_category.element_class'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.polynomial_padic_flat import \
                    Polynomial_padic_flat
            element_class = Polynomial_padic_flat
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring, name=name,
                                                         implementation=implementation,
                                                         element_class=element_class, category=category)


class PolynomialRing_dense_padic_ring_fixed_mod(PolynomialRing_dense_padic_ring_generic):
    def __init__(self, base_ring, name=None, implementation=None, element_class=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_ring_fixed_mod as PRing
            sage: R = PRing(Zp(13, type='fixed-mod'), name='t'); R                                                      # needs sage.rings.padics
            Univariate Polynomial Ring in t over 13-adic Ring of fixed modulus 13^20

            sage: type(R.gen())                                                                                         # needs sage.rings.padics
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_padic_ring_fixed_mod_with_category.element_class'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.polynomial_padic_flat import \
                    Polynomial_padic_flat
            element_class = Polynomial_padic_flat
        PolynomialRing_dense_padic_ring_generic.__init__(self, base_ring, name=name,
                                                         implementation=implementation,
                                                         element_class=element_class, category=category)


class PolynomialRing_dense_padic_field_capped_relative(PolynomialRing_dense_padic_field_generic):
    def __init__(self, base_ring, name=None, implementation=None, element_class=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_padic_field_capped_relative as PRing
            sage: R = PRing(Qp(13), name='t'); R                                                                        # needs sage.rings.padics
            Univariate Polynomial Ring in t over 13-adic Field with capped relative precision 20
            sage: type(R.gen())                                                                                         # needs sage.rings.padics
            <class 'sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_padic_field_capped_relative_with_category.element_class'>
        """
        if element_class is None:
            from sage.rings.polynomial.padics.\
                    polynomial_padic_capped_relative_dense import \
                    Polynomial_padic_capped_relative_dense
            element_class = Polynomial_padic_capped_relative_dense
        PolynomialRing_dense_padic_field_generic.__init__(self, base_ring, name=name,
                                                          implementation=implementation,
                                                          element_class=element_class, category=category)


class PolynomialRing_dense_mod_n(PolynomialRing_commutative):
    def __init__(self, base_ring, name=None, element_class=None,
                 implementation=None, category=None):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_mod_n as PRing
            sage: R = PRing(Zmod(15), 'x'); R
            Univariate Polynomial Ring in x over Ring of integers modulo 15
            sage: type(R.gen())                                                                                         # needs sage.libs.flint
            <class 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

            sage: R = PRing(Zmod(15), 'x', implementation='NTL'); R                                                     # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Ring of integers modulo 15 (using NTL)
            sage: type(R.gen())                                                                                         # needs sage.libs.ntl
            <class 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_zz'>

            sage: R = PRing(Zmod(2**63*3), 'x', implementation='NTL'); R                                                # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Ring of integers modulo 27670116110564327424 (using NTL)
            sage: type(R.gen())                                                                                         # needs sage.libs.ntl
            <class 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ'>

            sage: R = PRing(Zmod(2**63*3), 'x', implementation='FLINT')
            Traceback (most recent call last):
            ...
            ValueError: FLINT does not support modulus 27670116110564327424

            sage: R = PRing(Zmod(2**63*3), 'x'); R                                                                      # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Ring of integers modulo 27670116110564327424 (using NTL)
            sage: type(R.gen())                                                                                         # needs sage.libs.ntl
            <class 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ'>
        """
        if element_class is None:
            self._implementation_repr = ''
            given_implementation = implementation
            for implementation in self._implementation_names(implementation, base_ring):
                if implementation == "FLINT":
                    try:
                        from .polynomial_zmod_flint import Polynomial_zmod_flint as element_class
                    except ImportError:
                        if given_implementation:
                            raise
                        continue
                    self._implementation_repr = ''
                elif implementation == "NTL":
                    modulus = base_ring.order()
                    try:
                        from . import polynomial_modn_dense_ntl as modn_dense_ntl
                    except ImportError:
                        if given_implementation:
                            raise
                        continue
                    if modulus < ZZ(modn_dense_ntl.zz_p_max):
                        element_class = modn_dense_ntl.Polynomial_dense_modn_ntl_zz
                    else:
                        element_class = modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ
                    self._implementation_repr = ' (using NTL)'
                break

        PolynomialRing_commutative.__init__(self, base_ring, name=name, implementation=implementation,
                                            element_class=element_class, category=category)

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_mod_n
            sage: PolynomialRing_dense_mod_n._implementation_names_impl("FLINT", IntegerModRing(10), False)
            ['FLINT', None]
            sage: PolynomialRing_dense_mod_n._implementation_names_impl("NTL", IntegerModRing(10), False)
            ['NTL']
            sage: PolynomialRing_dense_mod_n._implementation_names_impl(None, IntegerModRing(10), False)
            ['FLINT', None]
            sage: PolynomialRing_dense_mod_n._implementation_names_impl("generic", IntegerModRing(10), False)
            NotImplemented
            sage: PolynomialRing_dense_mod_n._implementation_names_impl("FLINT", IntegerModRing(10^30), False)
            Traceback (most recent call last):
            ...
            ValueError: FLINT does not support modulus 1000000000000000000000000000000
            sage: PolynomialRing_dense_mod_n._implementation_names_impl("NTL", IntegerModRing(10^30), False)
            ['NTL', None]
            sage: PolynomialRing_dense_mod_n._implementation_names_impl(None, IntegerModRing(10^30), False)
            ['NTL', None]
            sage: PolynomialRing_dense_mod_n._implementation_names_impl("generic", IntegerModRing(10^30), False)
            NotImplemented
            sage: PolynomialRing_dense_mod_n._implementation_names_impl(None, IntegerModRing(10^30), True)
            NotImplemented
        """
        if sparse:
            return NotImplemented
        modulus = base_ring.order()
        if modulus <= sys.maxsize:
            defaults = ["FLINT", None]
        elif implementation == "FLINT":
            raise ValueError("FLINT does not support modulus %s" % modulus)
        else:
            defaults = ["NTL", None]
        if implementation in defaults:
            return defaults
        elif implementation == "NTL":
            return [implementation]
        return NotImplemented

    @cached_method
    def modulus(self):
        """
        EXAMPLES::

            sage: R.<x> = Zmod(15)[]
            sage: R.modulus()
            15
        """
        return self.base_ring().characteristic()

    def _repr_(self):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_integral_domain as PRing
            sage: R = PRing(ZZ, 'x', implementation='NTL'); R                           # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Integer Ring (using NTL)
        """
        s = PolynomialRing_commutative._repr_(self)
        return s + self._implementation_repr

    def residue_field(self, ideal, names=None):
        """
        Return the residue finite field at the given ideal.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<t> = GF(2)[]
            sage: k.<a> = R.residue_field(t^3 + t + 1); k
            Residue field in a
             of Principal ideal (t^3 + t + 1) of Univariate Polynomial Ring in t
             over Finite Field of size 2 (using GF2X)
            sage: k.list()
            [0, a, a^2, a + 1, a^2 + a, a^2 + a + 1, a^2 + 1, 1]
            sage: R.residue_field(t)
            Residue field of Principal ideal (t) of Univariate Polynomial Ring in t
             over Finite Field of size 2 (using GF2X)
            sage: P = R.irreducible_element(8) * R
            sage: P
            Principal ideal (t^8 + t^4 + t^3 + t^2 + 1) of Univariate Polynomial Ring in t
             over Finite Field of size 2 (using GF2X)
            sage: k.<a> = R.residue_field(P); k
            Residue field in a
             of Principal ideal (t^8 + t^4 + t^3 + t^2 + 1) of Univariate Polynomial Ring in t
             over Finite Field of size 2 (using GF2X)
            sage: k.cardinality()
            256

        Non-maximal ideals are not accepted::

            sage: # needs sage.libs.ntl
            sage: R.residue_field(t^2 + 1)
            Traceback (most recent call last):
            ...
            ArithmeticError: ideal is not maximal
            sage: R.residue_field(0)
            Traceback (most recent call last):
            ...
            ArithmeticError: ideal is not maximal
            sage: R.residue_field(1)
            Traceback (most recent call last):
            ...
            ArithmeticError: ideal is not maximal
        """
        ideal = self.ideal(ideal)
        if not ideal.is_maximal():
            raise ArithmeticError("ideal is not maximal")
        return ideal.residue_field(names)


class PolynomialRing_dense_mod_p(PolynomialRing_dense_finite_field,
                                 PolynomialRing_dense_mod_n,
                                 PolynomialRing_singular_repr):
    def __init__(self, base_ring, name='x', implementation=None, element_class=None, category=None):
        """
        TESTS::

            sage: P = GF(2)['x']; P                                                     # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
            sage: type(P.gen())                                                         # needs sage.libs.ntl
            <class 'sage.rings.polynomial.polynomial_gf2x.Polynomial_GF2X'>

            sage: from sage.rings.polynomial.polynomial_ring import PolynomialRing_dense_mod_p
            sage: P = PolynomialRing_dense_mod_p(GF(5), 'x'); P
            Univariate Polynomial Ring in x over Finite Field of size 5
            sage: type(P.gen())                                                         # needs sage.libs.flint
            <class 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

            sage: # needs sage.libs.ntl
            sage: P = PolynomialRing_dense_mod_p(GF(5), 'x', implementation='NTL'); P
            Univariate Polynomial Ring in x over Finite Field of size 5 (using NTL)
            sage: type(P.gen())
            <class 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>

            sage: P = PolynomialRing_dense_mod_p(GF(9223372036854775837), 'x'); P       # needs sage.libs.ntl sage.rings.finite_rings
            Univariate Polynomial Ring in x over Finite Field of size 9223372036854775837 (using NTL)
            sage: type(P.gen())                                                         # needs sage.libs.ntl sage.rings.finite_rings
            <class 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_mod_p'>

        This caching bug was fixed in :issue:`24264`::

            sage: # needs sage.rings.finite_rings
            sage: p = 2^64 + 13
            sage: A = GF(p^2)
            sage: B = GF(p^3)
            sage: R = A.modulus().parent()
            sage: S = B.modulus().parent()
            sage: R is S
            True
        """
        if element_class is None:
            given_implementation = implementation
            for implementation in self._implementation_names(implementation, base_ring):
                if implementation == "FLINT":
                    try:
                        from .polynomial_zmod_flint import Polynomial_zmod_flint as element_class
                    except ImportError:
                        if given_implementation:
                            raise
                        continue
                    self._implementation_repr = ''
                elif implementation == "NTL":
                    try:
                        from .polynomial_modn_dense_ntl import Polynomial_dense_mod_p as element_class
                    except ImportError:
                        if given_implementation:
                            raise
                        continue
                    self._implementation_repr = ' (using NTL)'
                elif implementation == "GF2X":
                    try:
                        from .polynomial_gf2x import Polynomial_GF2X as element_class
                    except ImportError:
                        if given_implementation:
                            raise
                        continue
                    self._implementation_repr = ' (using GF2X)'
                break

        category = check_default_category(PrincipalIdealDomains(),
                                          category)

        PolynomialRing_dense_mod_n.__init__(self, base_ring, name=name, implementation=implementation,
                                            element_class=element_class, category=category)

        self._has_singular = can_convert_to_singular(self)

    @staticmethod
    def _implementation_names_impl(implementation, base_ring, sparse):
        """
        TESTS::

            sage: # needs sage.libs.ntl
            sage: PolynomialRing(GF(2), 'x', implementation='GF2X')
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
            sage: PolynomialRing(GF(2), 'x', implementation='NTL')
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
            sage: PolynomialRing(GF(2), 'x', implementation=None)
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
            sage: PolynomialRing(GF(3), 'x', implementation='GF2X')
            Traceback (most recent call last):
            ...
            ValueError: GF2X only supports modulus 2
            sage: A = PolynomialRing(Zmod(2), 'x'); A
            Univariate Polynomial Ring in x over Ring of integers modulo 2 (using GF2X)
            sage: A in PrincipalIdealDomains()
            True

            sage: PolynomialRing(GF(2), 'x', implementation="FLINT")                    # needs sage.libs.flint
            Univariate Polynomial Ring in x over Finite Field of size 2
        """
        if sparse:
            return NotImplemented
        modulus = base_ring.characteristic()
        if modulus == 2:
            defaults = ["GF2X", "NTL", None]
        elif implementation == "GF2X":
            raise ValueError("GF2X only supports modulus 2")
        elif modulus <= sys.maxsize:
            defaults = ["FLINT", None]
        elif implementation == "FLINT":
            raise ValueError("FLINT does not support modulus %s" % modulus)
        else:
            defaults = ["NTL", None]
        if implementation in defaults:
            return defaults
        elif implementation in ["NTL", "FLINT"]:
            return [implementation]
        return NotImplemented

    def irreducible_element(self, n, algorithm=None):
        """
        Construct a monic irreducible polynomial of degree `n`.

        INPUT:

        - ``n`` -- integer; the degree of the polynomial to construct

        - ``algorithm`` -- string (algorithm to use) or ``None``;
          currently available options are:

          - ``'adleman-lenstra'``: a variant of the Adleman--Lenstra
            algorithm as implemented in PARI.

          - ``'conway'``: look up the Conway polynomial of degree `n`
            over the field of `p` elements in the database; raise a
            :exc:`RuntimeError` if it is not found.

          - ``'ffprimroot'``: use the :pari:`ffprimroot` function from
            PARI.

          - ``'first_lexicographic'``: return the lexicographically
            smallest irreducible polynomial of degree `n`.

          - ``'minimal_weight'``: return an irreducible polynomial of
            degree `n` with minimal number of non-zero coefficients.
            Only implemented for `p = 2`.

          - ``'primitive'``: return a polynomial `f` such that a root of
            `f` generates the multiplicative group of the finite field
            extension defined by `f`. This uses the Conway polynomial if
            possible, otherwise it uses ``'ffprimroot'``.

          - ``'random'``: try random polynomials until an irreducible
            one is found.

          If ``algorithm`` is ``None``, use `x - 1` in degree 1. In
          degree > 1, the Conway polynomial is used if it is found in
          the database.  Otherwise, the algorithm ``minimal_weight``
          is used if `p = 2`, and the algorithm ``'adleman-lenstra'`` if
          `p > 2`.

        OUTPUT: a monic irreducible polynomial of degree `n` in ``self``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: GF(5)['x'].irreducible_element(2)
            x^2 + 4*x + 2
            sage: GF(5)['x'].irreducible_element(2, algorithm='adleman-lenstra')
            x^2 + x + 1
            sage: GF(5)['x'].irreducible_element(2, algorithm='primitive')
            x^2 + 4*x + 2
            sage: GF(5)['x'].irreducible_element(32, algorithm='first_lexicographic')
            x^32 + 2
            sage: GF(5)['x'].irreducible_element(32, algorithm='conway')
            Traceback (most recent call last):
            ...
            RuntimeError: requested Conway polynomial not in database.
            sage: GF(5)['x'].irreducible_element(32, algorithm='primitive')
            x^32 + ...

        In characteristic 2::

            sage: GF(2)['x'].irreducible_element(33)                                    # needs sage.rings.finite_rings
            x^33 + x^13 + x^12 + x^11 + x^10 + x^8 + x^6 + x^3 + 1
            sage: GF(2)['x'].irreducible_element(33, algorithm='minimal_weight')        # needs sage.rings.finite_rings
            x^33 + x^10 + 1

        In degree 1::

            sage: GF(97)['x'].irreducible_element(1)                                    # needs sage.rings.finite_rings
            x + 96
            sage: GF(97)['x'].irreducible_element(1, algorithm='conway')                # needs sage.rings.finite_rings
            x + 92
            sage: GF(97)['x'].irreducible_element(1, algorithm='adleman-lenstra')       # needs sage.rings.finite_rings
            x

        AUTHORS:

        - Peter Bruin (June 2013)

        - Jeroen Demeyer (September 2014): add "ffprimroot" algorithm,
          see :issue:`8373`.
        """
        from sage.libs.pari.all import pari
        from sage.rings.finite_rings.conway_polynomials import (conway_polynomial,
                                                                exists_conway_polynomial)

        p = self.characteristic()
        n = int(n)
        if n < 1:
            raise ValueError("degree must be at least 1")

        if algorithm is None:
            if n == 1:
                return self((-1,1))  # Polynomial x - 1
            elif exists_conway_polynomial(p, n):
                algorithm = "conway"
            elif p == 2:
                try:
                    from .polynomial_gf2x import GF2X_BuildSparseIrred_list
                except ImportError:
                    algorithm = "adleman-lenstra"
                else:
                    algorithm = "minimal_weight"
            else:
                algorithm = "adleman-lenstra"
        elif algorithm == "primitive":
            if exists_conway_polynomial(p, n):
                algorithm = "conway"
            else:
                algorithm = "ffprimroot"

        if algorithm == "adleman-lenstra":
            return self(pari(p).ffinit(n))
        elif algorithm == "conway":
            return self(conway_polynomial(p, n))
        elif algorithm == "first_lexicographic":
            if p == 2:
                try:
                    from .polynomial_gf2x import GF2X_BuildIrred_list
                except ImportError:
                    pass
                else:
                    return self(GF2X_BuildIrred_list(n))
            else:
                # Fallback to PolynomialRing_dense_finite_field.irreducible_element
                pass
        elif algorithm == "ffprimroot":
            return self(pari(p).ffinit(n).ffgen().ffprimroot().charpoly())
        elif algorithm == "minimal_weight":
            if p == 2:
                from .polynomial_gf2x import GF2X_BuildSparseIrred_list
                return self(GF2X_BuildSparseIrred_list(n))
            else:
                raise NotImplementedError("'minimal_weight' option only implemented for p = 2")
        elif algorithm == "random":
            if p == 2:
                try:
                    from .polynomial_gf2x import GF2X_BuildRandomIrred_list
                except ImportError:
                    pass
                else:
                    return self(GF2X_BuildRandomIrred_list(n))
            else:
                pass

        # No suitable algorithm found, try algorithms from the base class.
        return PolynomialRing_dense_finite_field.irreducible_element(self, n, algorithm)

    @cached_method
    def fraction_field(self):
        """
        Return the fraction field of ``self``.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: R.fraction_field()
            Fraction Field of Univariate Polynomial Ring in t
             over Finite Field of size 5
        """
        try:
            from sage.rings.fraction_field_FpT import FpT
            from sage.rings.polynomial.polynomial_zmod_flint import Polynomial_zmod_flint
        except ImportError:
            pass
        else:
            p = self.base_ring().characteristic()
            if (issubclass(self.element_class, Polynomial_zmod_flint)
                    and 2 < p < FpT.INTEGER_LIMIT):
                return FpT(self)
        return super().fraction_field()


def polygen(ring_or_element, name='x'):
    """
    Return a polynomial indeterminate.

    INPUT:

    - ``polygen(base_ring, name='x')``

    - ``polygen(ring_element, name='x')``

    If the first input is a ring, return a polynomial generator over
    that ring. If it is a ring element, return a polynomial generator
    over the parent of the element.

    EXAMPLES::

        sage: z = polygen(QQ, 'z')
        sage: z^3 + z +1
        z^3 + z + 1
        sage: parent(z)
        Univariate Polynomial Ring in z over Rational Field

    .. NOTE::

       If you give a list or comma-separated string to :func:`polygen`, you'll
       get a tuple of indeterminates, exactly as if you called
       :func:`polygens`.
    """
    if isinstance(ring_or_element, RingElement):
        base_ring = ring_or_element.parent()
    elif ring_or_element in Rings():
        base_ring = ring_or_element
    else:
        raise TypeError("input must be a ring or ring element")
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    t = PolynomialRing(base_ring, name)
    if t.ngens() > 1:
        return t.gens()
    return t.gen()


def polygens(base_ring, names='x', *args):
    """
    Return indeterminates over the given base ring with the given
    names.

    EXAMPLES::

        sage: x,y,z = polygens(QQ,'x,y,z')
        sage: (x+y+z)^2
        x^2 + 2*x*y + y^2 + 2*x*z + 2*y*z + z^2
        sage: parent(x)
        Multivariate Polynomial Ring in x, y, z over Rational Field
        sage: t = polygens(QQ, ['x','yz','abc'])
        sage: t
        (x, yz, abc)

    The number of generators can be passed as a third argument::

        sage: polygens(QQ, 'x', 4)
        (x0, x1, x2, x3)
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    return PolynomialRing(base_ring, names, *args).gens()
