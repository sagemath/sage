"""
Laurent Series Rings

EXAMPLES::

    sage: R = LaurentSeriesRing(QQ, "x")
    sage: R.base_ring()
    Rational Field
    sage: S = LaurentSeriesRing(GF(17)['x'], 'y')
    sage: S
    Laurent Series Ring in y over Univariate Polynomial Ring in x over
    Finite Field of size 17
    sage: S.base_ring()
    Univariate Polynomial Ring in x over Finite Field of size 17

.. SEEALSO::

    * :func:`sage.misc.defaults.set_series_precision`
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2007 Robert Bradshaw <robertwb@math.washington.edu>
#                     2012 David Roe <roed.math@gmail.com>
#                     2014 Peter Bruin <P.J.Bruin@math.leidenuniv.nl>
#                     2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.categories.rings import Rings
from sage.rings.infinity import infinity
from sage.categories.algebras import Algebras
from sage.categories.integral_domains import IntegralDomains
from sage.categories.fields import Fields
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from .laurent_series_ring_element import LaurentSeries
from .ring import CommutativeRing

from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ


def is_LaurentSeriesRing(x):
    """
    Return ``True`` if this is a *univariate* Laurent series ring.

    This is in keeping with the behavior of ``is_PolynomialRing``
    versus ``is_MPolynomialRing``.

    TESTS::

        sage: from sage.rings.laurent_series_ring import is_LaurentSeriesRing
        sage: K.<q> = LaurentSeriesRing(QQ)
        sage: is_LaurentSeriesRing(K)
        True
    """
    return isinstance(x, LaurentSeriesRing)


class LaurentSeriesRing(UniqueRepresentation, CommutativeRing):
    r"""
    Univariate Laurent Series Ring.

    EXAMPLES::

        sage: R = LaurentSeriesRing(QQ, 'x'); R
        Laurent Series Ring in x over Rational Field
        sage: x = R.0
        sage: g = 1 - x + x^2 - x^4 +O(x^8); g
        1 - x + x^2 - x^4 + O(x^8)
        sage: g = 10*x^(-3) + 2006 - 19*x + x^2 - x^4 +O(x^8); g
        10*x^-3 + 2006 - 19*x + x^2 - x^4 + O(x^8)

    You can also use more mathematical notation when the base is a
    field::

        sage: Frac(QQ[['x']])
        Laurent Series Ring in x over Rational Field
        sage: Frac(GF(5)['y'])
        Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5

    When the base ring is a domain, the fraction field is the
    Laurent series ring over the fraction field of the base ring::

        sage: Frac(ZZ[['t']])
        Laurent Series Ring in t over Rational Field

    Laurent series rings are determined by their variable and the base
    ring, and are globally unique::

        sage: K = Qp(5, prec = 5)
        sage: L = Qp(5, prec = 200)
        sage: R.<x> = LaurentSeriesRing(K)
        sage: S.<y> = LaurentSeriesRing(L)
        sage: R is S
        False
        sage: T.<y> = LaurentSeriesRing(Qp(5,prec=200))
        sage: S is T
        True
        sage: W.<y> = LaurentSeriesRing(Qp(5,prec=199))
        sage: W is T
        False

        sage: K = LaurentSeriesRing(CC, 'q')
        sage: K
        Laurent Series Ring in q over Complex Field with 53 bits of precision
        sage: loads(K.dumps()) == K
        True
        sage: P = QQ[['x']]
        sage: F = Frac(P)
        sage: TestSuite(F).run()

    When the base ring `k` is a field, the ring `k((x))` is a CDVF, that is
    a field equipped with a discrete valuation for which it is complete.
    The appropriate (sub)category is automatically set in this case::

        sage: k = GF(11)
        sage: R.<x> = k[[]]
        sage: F = Frac(R)
        sage: F.category()
        Join of Category of complete discrete valuation fields and Category of commutative algebras over (finite enumerated fields and subquotients of monoids and quotients of semigroups) and Category of infinite sets
        sage: TestSuite(F).run()

    TESTS:

    Check if changing global series precision does it right (and
    that :trac:`17955` is fixed)::

        sage: set_series_precision(3)
        sage: R.<x> = LaurentSeriesRing(ZZ)
        sage: 1/(1 - 2*x)
        1 + 2*x + 4*x^2 + O(x^3)
        sage: set_series_precision(5)
        sage: R.<x> = LaurentSeriesRing(ZZ)
        sage: 1/(1 - 2*x)
        1 + 2*x + 4*x^2 + 8*x^3 + 16*x^4 + O(x^5)
        sage: set_series_precision(20)

    Check categories (:trac:`24420`)::

        sage: LaurentSeriesRing(ZZ, 'x').category()
        Category of infinite commutative no zero divisors algebras over (euclidean domains and infinite enumerated sets and metric spaces)
        sage: LaurentSeriesRing(QQ, 'x').category()
        Join of Category of complete discrete valuation fields and Category of commutative algebras over (number fields and quotient fields and metric spaces) and Category of infinite sets
        sage: LaurentSeriesRing(Zmod(4), 'x').category()
        Category of infinite commutative algebras over (finite commutative rings and subquotients of monoids and quotients of semigroups and finite enumerated sets)

    Check coercions (:trac:`24431`)::

        sage: pts = [LaurentSeriesRing,
        ....:        PolynomialRing,
        ....:        PowerSeriesRing,
        ....:        LaurentPolynomialRing]
        sage: LS = LaurentSeriesRing(QQ, 'x')
        sage: LSx = LS.gen()

        sage: for P in pts:
        ....:     x = P(QQ, 'x').gen()
        ....:     assert parent(LSx * x) is LS, "wrong parent for {}".format(P)

        sage: for P in pts:
        ....:     y = P(QQ, 'y').gen()
        ....:     try:
        ....:         LSx * y
        ....:     except TypeError:
        ....:         pass
        ....:     else:
        ....:         print("wrong coercion {}".format(P))
    """
    Element = LaurentSeries

    @staticmethod
    def __classcall__(cls, *args, **kwds):
        r"""
        TESTS::

            sage: L = LaurentSeriesRing(QQ, 'q')
            sage: L is LaurentSeriesRing(QQ, name='q')
            True
            sage: loads(dumps(L)) is L
            True

            sage: L.variable_names()
            ('q',)
            sage: L.variable_name()
            'q'
        """
        from .power_series_ring import PowerSeriesRing, is_PowerSeriesRing

        if not kwds and len(args) == 1 and is_PowerSeriesRing(args[0]):
            power_series = args[0]
        else:
            power_series = PowerSeriesRing(*args, **kwds)

        return UniqueRepresentation.__classcall__(cls, power_series)

    def __init__(self, power_series):
        """
        Initialization

        EXAMPLES::

            sage: K.<q> = LaurentSeriesRing(QQ, default_prec=4); K
            Laurent Series Ring in q over Rational Field
            sage: 1 / (q-q^2)
            q^-1 + 1 + q + q^2 + O(q^3)

            sage: RZZ = LaurentSeriesRing(ZZ, 't')
            sage: RZZ.category()
            Category of infinite commutative no zero divisors algebras over (euclidean domains and infinite enumerated sets and metric spaces)
            sage: TestSuite(RZZ).run()

            sage: R1 = LaurentSeriesRing(Zmod(1), 't')
            sage: R1.category()
            Category of finite commutative algebras over (finite commutative rings and subquotients of monoids and quotients of semigroups and finite enumerated sets)
            sage: TestSuite(R1).run()

            sage: R2 = LaurentSeriesRing(Zmod(2), 't')
            sage: R2.category()
            Join of Category of complete discrete valuation fields and Category of commutative algebras over (finite enumerated fields and subquotients of monoids and quotients of semigroups) and Category of infinite sets
            sage: TestSuite(R2).run()

            sage: R4 = LaurentSeriesRing(Zmod(4), 't')
            sage: R4.category()
            Category of infinite commutative algebras over (finite commutative rings and subquotients of monoids and quotients of semigroups and finite enumerated sets)
            sage: TestSuite(R4).run()

            sage: RQQ = LaurentSeriesRing(QQ, 't')
            sage: RQQ.category()
            Join of Category of complete discrete valuation fields and Category of commutative algebras over (number fields and quotient fields and metric spaces) and Category of infinite sets
            sage: TestSuite(RQQ).run()
        """
        base_ring = power_series.base_ring()
        category = Algebras(base_ring.category())
        if base_ring in Fields():
            category &= CompleteDiscreteValuationFields()
        elif base_ring in IntegralDomains():
            category &= IntegralDomains()
        elif base_ring in Rings().Commutative():
            category = category.Commutative()

        if base_ring.is_zero():
            category = category.Finite()
        else:
            category = category.Infinite()

        self._power_series_ring = power_series
        self._one_element = self.element_class(self, power_series.one())
        CommutativeRing.__init__(self, base_ring,
                names=power_series.variable_names(),
                category=category)

    def base_extend(self, R):
        """
        Return the Laurent series ring over R in the same variable as
        self, assuming there is a canonical coerce map from the base ring
        of self to R.

        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, default_prec=4)
            sage: K.base_extend(QQ['t'])
            Laurent Series Ring in x over Univariate Polynomial Ring in t over Rational Field
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        else:
            raise TypeError("no valid base extension defined")

    def fraction_field(self):
        r"""
        Return the fraction field of this ring of Laurent series.

        If the base ring is a field, then Laurent series are already a field.
        If the base ring is a domain, then the Laurent series over its fraction
        field is returned. Otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: R = LaurentSeriesRing(ZZ, 't', 30).fraction_field()
            sage: R
            Laurent Series Ring in t over Rational Field
            sage: R.default_prec()
            30

            sage: LaurentSeriesRing(Zmod(4), 't').fraction_field()
            Traceback (most recent call last):
            ...
            ValueError: must be an integral domain
        """
        from sage.categories.integral_domains import IntegralDomains
        from sage.categories.fields import Fields
        if self in Fields():
            return self
        elif self in IntegralDomains():
            return LaurentSeriesRing(self.base_ring().fraction_field(),
                    self.variable_names(),
                    self.default_prec())
        else:
            raise ValueError('must be an integral domain')

    def change_ring(self, R):
        """
        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, default_prec=4)
            sage: R = K.change_ring(ZZ); R
            Laurent Series Ring in x over Integer Ring
            sage: R.default_prec()
            4
        """
        return LaurentSeriesRing(R, self.variable_names(),
                                 default_prec=self.default_prec(),
                                 sparse=self.is_sparse())

    def is_sparse(self):
        """
        Return if ``self`` is a sparse implementation.

        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, sparse=True)
            sage: K.is_sparse()
            True
        """
        return self.power_series_ring().is_sparse()

    def is_field(self, proof=True):
        """
        A Laurent series ring is a field if and only if the base ring
        is a field.

        TESTS::

            sage: LaurentSeriesRing(QQ,'t').is_field()
            True
            sage: LaurentSeriesRing(ZZ,'t').is_field()
            False
        """
        return self.base_ring().is_field()

    def is_dense(self):
        """
        EXAMPLES::

            sage: K.<x> = LaurentSeriesRing(QQ, sparse=True)
            sage: K.is_dense()
            False
        """
        return self.power_series_ring().is_dense()

    def _repr_(self):
        """
        EXAMPLES::

            sage: LaurentSeriesRing(QQ,'q') # indirect doctest
            Laurent Series Ring in q over Rational Field
            sage: LaurentSeriesRing(ZZ,'t',sparse=True)
            Sparse Laurent Series Ring in t over Integer Ring
        """
        s = "Laurent Series Ring in %s over %s" % (self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def _element_constructor_(self, x, n=0, prec=infinity):
        r"""
        Construct a Laurent series from `x`.

        INPUT:

        - ``x`` -- object that can be converted into a Laurent series

        - ``n`` -- (default: 0) multiply the result by `t^n`

        - ``prec`` -- (default: ``infinity``) the precision of the series
            as an integer.

        EXAMPLES::

            sage: R.<u> = LaurentSeriesRing(Qp(5, 10))
            sage: S.<t> = LaurentSeriesRing(RationalField())
            sage: R(t + t^2 + O(t^3))
            (1 + O(5^10))*u + (1 + O(5^10))*u^2 + O(u^3)
            sage: R(t + t^2 + O(t^3), prec=2)
            (1 + O(5^10))*u + O(u^2)

        Coercing an element into its own parent produces that element
        again, unless a different ``n`` or ``prec`` is given::

            sage: u is R(u)
            True
            sage: R(u, n=3, prec=7)
            (1 + O(5^10))*u^4 + O(u^7)

        Rational functions are accepted::

            sage: I = sqrt(-1)
            sage: K.<I> = QQ[I]
            sage: P.<t> = PolynomialRing(K)
            sage: L.<u> = LaurentSeriesRing(QQ[I])
            sage: L((t*I)/(t^3+I*2*t))
            1/2 + 1/4*I*u^2 - 1/8*u^4 - 1/16*I*u^6 + 1/32*u^8 +
            1/64*I*u^10 - 1/128*u^12 - 1/256*I*u^14 + 1/512*u^16 +
            1/1024*I*u^18 + O(u^20)

        ::

            sage: L(t*I) / L(t^3+I*2*t)
            1/2 + 1/4*I*u^2 - 1/8*u^4 - 1/16*I*u^6 + 1/32*u^8 +
            1/64*I*u^10 - 1/128*u^12 - 1/256*I*u^14 + 1/512*u^16 +
            1/1024*I*u^18 + O(u^20)

        TESTS:

        Check that :trac:`28993` is fixed::

            sage: from sage.modular.etaproducts import qexp_eta
            sage: qexp_eta(S, prec=30)
            1 - t - t^2 + t^5 + t^7 - t^12 - t^15 + t^22 + t^26 + O(t^30)

        When converting from `R((z))` to `R((z))((w))`, the variable
        `z` is sent to `z` rather than to `w` (see :trac:`7085`)::

            sage: A.<z> = LaurentSeriesRing(QQ)
            sage: B.<w> = LaurentSeriesRing(A)
            sage: B(z)
            z
            sage: z/w
            z*w^-1

        Various conversions from PARI (see also :trac:`2508`)::

            sage: L.<q> = LaurentSeriesRing(QQ, default_prec=10)
            sage: L(pari('1/x'))
            q^-1
            sage: L(pari('polchebyshev(5)'))
            5*q - 20*q^3 + 16*q^5
            sage: L(pari('polchebyshev(5) - 1/x^4'))
            -q^-4 + 5*q - 20*q^3 + 16*q^5
            sage: L(pari('1/polchebyshev(5)'))
            1/5*q^-1 + 4/5*q + 64/25*q^3 + 192/25*q^5 + 2816/125*q^7 + O(q^9)
            sage: L(pari('polchebyshev(5) + O(x^40)'))
            5*q - 20*q^3 + 16*q^5 + O(q^40)
            sage: L(pari('polchebyshev(5) - 1/x^4 + O(x^40)'))
            -q^-4 + 5*q - 20*q^3 + 16*q^5 + O(q^40)
            sage: L(pari('1/polchebyshev(5) + O(x^10)'))
            1/5*q^-1 + 4/5*q + 64/25*q^3 + 192/25*q^5 + 2816/125*q^7 + 8192/125*q^9 + O(q^10)
            sage: L(pari('1/polchebyshev(5) + O(x^10)'), -10)  # Multiply by q^-10
            1/5*q^-11 + 4/5*q^-9 + 64/25*q^-7 + 192/25*q^-5 + 2816/125*q^-3 + 8192/125*q^-1 + O(1)
            sage: L(pari('O(x^-10)'))
            O(q^-10)

        Check that :trac:`30073` is fixed::

            sage: P.<x> = LaurentSeriesRing(QQ)
            sage: P({-3: 1})
            x^-3
        """
        from sage.rings.fraction_field_element import is_FractionFieldElement
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
        from sage.structure.element import parent
        from sage.libs.pari.all import pari_gen

        P = parent(x)
        if isinstance(x, self.element_class) and n == 0 and P is self:
            return x.add_bigoh(prec)  # ok, since Laurent series are immutable (no need to make a copy)
        elif P is self.base_ring():
            # Convert x into a power series; if P is itself a Laurent
            # series ring A((t)), this prevents the implementation of
            # LaurentSeries.__init__() from effectively applying the
            # ring homomorphism A((t)) -> A((t))((u)) sending t to u
            # instead of the one sending t to t.  We cannot easily
            # tell LaurentSeries.__init__() to be more strict, because
            # A((t)) -> B((u)) is expected to send t to u if A admits
            # a coercion to B but A((t)) does not, and this condition
            # would be inefficient to check there.
            x = self.power_series_ring()(x)
        elif isinstance(x, pari_gen):
            t = x.type()
            if t == "t_RFRAC":   # Rational function
                x = self(self.polynomial_ring()(x.numerator())) / \
                    self(self.polynomial_ring()(x.denominator()))
                return (x << n).add_bigoh(prec)
            elif t == "t_SER":   # Laurent series
                n += x._valp()
                bigoh = n + x.length()
                x = self(self.polynomial_ring()(x.Vec()))
                return (x << n).add_bigoh(bigoh)
            else:  # General case, pretend to be a polynomial
                return (self(self.polynomial_ring()(x)) << n).add_bigoh(prec)
        elif (is_FractionFieldElement(x)
              and (x.base_ring() is self.base_ring() or x.base_ring() == self.base_ring())
              and (is_Polynomial(x.numerator()) or is_MPolynomial(x.numerator()))):
            x = self(x.numerator()) / self(x.denominator())
            return (x << n).add_bigoh(prec)
        return self.element_class(self, x, n).add_bigoh(prec)

    def random_element(self, algorithm='default'):
        r"""
        Return a random element of this Laurent series ring.

        The optional ``algorithm`` parameter decides how elements are generated.
        Algorithms currently implemented:

        - ``'default'``: Choose an integer ``shift`` using the standard
          distribution on the integers.  Then choose a list of coefficients
          using the ``random_element`` function of the base ring, and construct
          a new element based on those coefficients, so that the i-th
          coefficient corresponds to the (i+shift)-th power of the uniformizer.
          The amount of coefficients is determined by the ``default_prec``
          of the ring. Note that this method only creates non-exact elements.

        EXAMPLES::

            sage: S.<s> = LaurentSeriesRing(GF(3))
            sage: S.random_element()  # random
            s^-8 + s^-7 + s^-6 + s^-5 + s^-1 + s + s^3 + s^4
            + s^5 + 2*s^6 + s^7 + s^11 + O(s^12)
        """
        if algorithm == 'default':
            shift = ZZ.random_element()
            return self([self.base_ring().random_element()
                         for k in range(self.default_prec())],
                        shift).O(shift + self.default_prec())
        else:
            raise ValueError("algorithm cannot be %s" % algorithm)

    def construction(self):
        r"""
        Return the functorial construction of this Laurent power series ring.

        The construction is given as the completion of the Laurent polynomials.

        EXAMPLES::

            sage: L.<t> = LaurentSeriesRing(ZZ, default_prec=42)
            sage: phi, arg = L.construction()
            sage: phi
            Completion[t, prec=42]
            sage: arg
            Univariate Laurent Polynomial Ring in t over Integer Ring
            sage: phi(arg) is L
            True

        Because of this construction, pushout is automatically available::

            sage: 1/2 * t
            1/2*t
            sage: parent(1/2 * t)
            Laurent Series Ring in t over Rational Field

            sage: QQbar.gen() * t
            I*t
            sage: parent(QQbar.gen() * t)
            Laurent Series Ring in t over Algebraic Field
        """
        from sage.categories.pushout import CompletionFunctor
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        L = LaurentPolynomialRing(self.base_ring(), self._names[0])
        return CompletionFunctor(self._names[0], self.default_prec()), L

    def _coerce_map_from_(self, P):
        """
        Return a coercion map from `P` to ``self``, or True, or None.

        The following rings admit a coercion map to the Laurent series
        ring `A((t))`:

        - any ring that admits a coercion map to `A` (including `A`
          itself);

        - any Laurent series ring, power series ring or polynomial
          ring in the variable `t` over a ring admitting a coercion
          map to `A`.

        EXAMPLES::

            sage: S.<t> = LaurentSeriesRing(ZZ)
            sage: S.has_coerce_map_from(ZZ)
            True
            sage: S.has_coerce_map_from(PolynomialRing(ZZ, 't'))
            True
            sage: S.has_coerce_map_from(LaurentPolynomialRing(ZZ, 't'))
            True
            sage: S.has_coerce_map_from(PowerSeriesRing(ZZ, 't'))
            True
            sage: S.has_coerce_map_from(S)
            True

            sage: S.has_coerce_map_from(QQ)
            False
            sage: S.has_coerce_map_from(PolynomialRing(QQ, 't'))
            False
            sage: S.has_coerce_map_from(LaurentPolynomialRing(QQ, 't'))
            False
            sage: S.has_coerce_map_from(PowerSeriesRing(QQ, 't'))
            False
            sage: S.has_coerce_map_from(LaurentSeriesRing(QQ, 't'))
            False

            sage: R.<t> = LaurentSeriesRing(QQ['x'])
            sage: R.has_coerce_map_from(QQ[['t']])
            True
            sage: R.has_coerce_map_from(QQ['t'])
            True
            sage: R.has_coerce_map_from(ZZ['x']['t'])
            True
            sage: R.has_coerce_map_from(ZZ['t']['x'])
            False
            sage: R.has_coerce_map_from(ZZ['x'])
            True
        """
        A = self.base_ring()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.power_series_ring import is_PowerSeriesRing
        from sage.rings.polynomial.laurent_polynomial_ring import is_LaurentPolynomialRing
        if ((is_LaurentSeriesRing(P) or
             is_LaurentPolynomialRing(P) or
             is_PowerSeriesRing(P) or
             is_PolynomialRing(P))
            and P.variable_name() == self.variable_name()
            and A.has_coerce_map_from(P.base_ring())):
            return True

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(GF(17))
            sage: S.<y> = LaurentSeriesRing(GF(19))
            sage: R.hom([y], S) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: relations do not all (canonically) map to 0 under map determined by images of generators
            sage: f = R.hom(x+x^3,R)
            sage: f(x^2)
            x^2 + 2*x^4 + x^6

        The image of the generator needs to be a unit::

            sage: R.<x> = LaurentSeriesRing(ZZ)
            sage: R._is_valid_homomorphism_(R, [2*x])
            False
        """
        # NOTE: There are no ring homomorphisms from the ring of
        # all formal power series to most rings, e.g, the p-adic
        # field, since you can always (mathematically!) construct
        # some power series that does not converge.
        # NOTE: The above claim is wrong when the base ring is Z.
        # See trac 28486.

        if base_map is None and not codomain.has_coerce_map_from(self.base_ring()):
            return False
        # Note that 0 is not a *ring* homomorphism, and you cannot map to a power series ring
        if is_LaurentSeriesRing(codomain):
            return im_gens[0].valuation() > 0 and im_gens[0].is_unit()
        return False

    def characteristic(self):
        """
        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(GF(17))
            sage: R.characteristic()
            17
        """
        return self.base_ring().characteristic()

    def residue_field(self):
        """
        Return the residue field of this Laurent series field
        if it is a complete discrete valuation field (i.e. if
        the base ring is a field, in which base it is also the
        residue field).

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(GF(17))
            sage: R.residue_field()
            Finite Field of size 17

            sage: R.<x> = LaurentSeriesRing(ZZ)
            sage: R.residue_field()
            Traceback (most recent call last):
            ...
            TypeError: the base ring is not a field
        """
        if not self.base_ring().is_field():
            raise TypeError("the base ring is not a field")
        return self.base_ring()

    def default_prec(self):
        """
        Get the precision to which exact elements are truncated when
        necessary (most frequently when inverting).

        EXAMPLES::

            sage: R.<x> = LaurentSeriesRing(QQ, default_prec=5)
            sage: R.default_prec()
            5
        """
        return self._power_series_ring.default_prec()

    def is_exact(self):
        """
        Laurent series rings are inexact.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.is_exact()
            False
        """
        return False

    @cached_method
    def gen(self, n=0):
        """
        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.gen()
            x
        """
        if n != 0:
            raise IndexError("generator {} not defined".format(n))
        return self.element_class(self, [0, 1])

    def uniformizer(self):
        """
        Return a uniformizer of this Laurent series field if it is
        a discrete valuation field (i.e. if the base ring is actually
        a field). Otherwise, an error is raised.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: R.uniformizer()
            t

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: R.uniformizer()
            Traceback (most recent call last):
            ...
            TypeError: the base ring is not a field
        """
        if not self.base_ring().is_field():
            raise TypeError("the base ring is not a field")
        return self.gen()

    def ngens(self):
        """
        Laurent series rings are univariate.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.ngens()
            1
        """
        return 1

    def polynomial_ring(self):
        r"""
        If this is the Laurent series ring `R((t))`, return the
        polynomial ring `R[t]`.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.polynomial_ring()
            Univariate Polynomial Ring in x over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self.base_ring(), self.variable_name(),
                              sparse=self.is_sparse())

    def laurent_polynomial_ring(self):
        r"""
        If this is the Laurent series ring `R((t))`, return the Laurent
        polynomial ring `R[t,1/t]`.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.laurent_polynomial_ring()
            Univariate Laurent Polynomial Ring in x over Rational Field
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        return LaurentPolynomialRing(self.base_ring(), self.variable_name(),
                                     sparse=self.is_sparse())

    def power_series_ring(self):
        r"""
        If this is the Laurent series ring `R((t))`, return the
        power series ring `R[[t]]`.

        EXAMPLES::

            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.power_series_ring()
            Power Series Ring in x over Rational Field
        """
        return self._power_series_ring
