"""
Power Series Methods

The class ``PowerSeries_poly`` provides additional methods for univariate power series.
"""
from sage.rings.power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element
from sage.rings.infinity import infinity

try:
    from sage.libs.pari.all import pari_gen, PariError
except ImportError:
    pari_gen = ()
    PariError = ()


cdef class PowerSeries_poly(PowerSeries):

    def __init__(self, parent, f=0, prec=infinity, int check=1, is_gen=0):
        """
        EXAMPLES::

            sage: R.<q> = PowerSeriesRing(CC); R                                        # needs sage.rings.real_mpfr
            Power Series Ring in q over Complex Field with 53 bits of precision
            sage: loads(q.dumps()) == q                                                 # needs sage.rings.real_mpfr
            True

            sage: R.<t> = QQ[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: a = f^3; a
            27 - 27*t^3 + O(t^5)
            sage: b = f^-3; b
            1/27 + 1/27*t^3 + O(t^5)
            sage: a*b
            1 + O(t^5)

        Check that :issue:`22216` is fixed::

            sage: R.<T> = PowerSeriesRing(QQ)
            sage: R(pari('1 + O(T)'))                                                   # needs sage.libs.pari
            1 + O(T)
            sage: R(pari('1/T + O(T)'))                                                 # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: series has negative valuation
        """
        R = parent._poly_ring()
        if isinstance(f, Element):
            if (<Element>f)._parent is R:
                pass
            elif (<Element>f)._parent == R.base_ring():
                f = R([f])
            elif isinstance(f, PowerSeries):  # not only PowerSeries_poly
                prec = (<PowerSeries>f)._prec
                f = R(f.polynomial())
            else:
                if f:
                    f = R(f, check=check)
                else:
                    f = R(None)
        elif isinstance(f, pari_gen) and f.type() == 't_SER':
            if f._valp() < 0:
                raise ValueError('series has negative valuation')
            if prec is infinity:
                prec = f.length() + f._valp()
            f = R(f.truncate())
        else:
            if f:
                f = R(f, check=check)
            else: # None is supposed to yield zero
                f = R(None)

        self.__f = f
        if check and not (prec is infinity):
            self.__f = self.__f.truncate(prec)
        PowerSeries.__init__(self, parent, prec, is_gen)

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: hash(t) == hash(R.gen())
            True
            sage: hash(t) != hash(R.one())
            True
        """
        return hash(self.__f)

    def __reduce__(self):
        """
        Used for pickling.

        EXAMPLES::

            sage: A.<z> = RR[[]]
            sage: f = z - z^3 + O(z^10)
            sage: f == loads(dumps(f)) # indirect doctest
            True
        """
        return self.__class__, (self._parent, self.__f, self._prec, self._is_gen)

    def polynomial(self):
        """
        Return the underlying polynomial of ``self``.

        EXAMPLES::

            sage: R.<t> = GF(7)[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: f.polynomial()
            6*t^3 + 3
        """
        return self.__f

    def valuation(self):
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: (5 - t^8 + O(t^11)).valuation()
            0
            sage: (-t^8 + O(t^11)).valuation()
            8
            sage: O(t^7).valuation()
            7
            sage: R(0).valuation()
            +Infinity
        """
        if self.__f == 0:
            return self._prec

        return self.__f.valuation()

    def degree(self):
        """
        Return the degree of the underlying polynomial of ``self``.

        That is, if ``self`` is of the form `f(x) + O(x^n)`, we return
        the degree of `f(x)`. Note that if `f(x)` is `0`, we return `-1`,
        just as with polynomials.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: (5 + t^3 + O(t^4)).degree()
            3
            sage: (5 + O(t^4)).degree()
            0
            sage: O(t^4).degree()
            -1
        """
        return self.__f.degree()

    def __bool__(self):
        """
        Return ``True`` if ``self`` is nonzero, and ``False`` otherwise.

        EXAMPLES::

            sage: R.<t> = GF(11)[[]]
            sage: bool(1 + t + O(t^18))
            True
            sage: bool(R(0))
            False
            sage: bool(O(t^18))
            False
        """
        return not not self.__f

    def __call__(self, *x, **kwds):
        """
        Evaluate the series at `x=a`.

        INPUT:

        - ``x``:

          - a tuple of elements the first of which can be meaningfully
            substituted in ``self``, with the remainder used for substitution
            in the coefficients of ``self``.

          - a dictionary for kwds:value pairs. If the variable name of
            ``self`` is a keyword it is substituted for.  Other keywords
            are used for substitution in the coefficients of ``self``.

        OUTPUT: the value of ``self`` after substitution

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: f = t^2 + t^3 + O(t^6)
            sage: f(t^3)
            t^6 + t^9 + O(t^18)
            sage: f(t=t^3)
            t^6 + t^9 + O(t^18)
            sage: f(f)
            t^4 + 2*t^5 + 2*t^6 + 3*t^7 + O(t^8)
            sage: f(f)(f) == f(f(f))
            True

        The following demonstrates that the problems raised in :issue:`3979`
        and :issue:`5367` are solved::

            sage: [f(t^2 + O(t^n)) for n in [9, 10, 11]]
            [t^4 + t^6 + O(t^11), t^4 + t^6 + O(t^12), t^4 + t^6 + O(t^12)]
            sage: f(t^2)
            t^4 + t^6 + O(t^12)

        It is possible to substitute a series for which only the precision
        is defined::

            sage: f(O(t^5))
            O(t^10)

        or to substitute a polynomial (the result belonging to the power
        series ring over the same base ring)::

            sage: P.<z> = ZZ[]
            sage: g = f(z + z^3); g
            z^2 + z^3 + 2*z^4 + 3*z^5 + O(z^6)
            sage: g.parent()
            Power Series Ring in z over Integer Ring

        A series defined over another ring can be substituted::

            sage: S.<u> = GF(7)[[]]
            sage: f(2*u + u^3 + O(u^5))
            4*u^2 + u^3 + 4*u^4 + 5*u^5 + O(u^6)

        As can a `p`-adic integer as long as the coefficient ring is compatible::

            sage: f(100 + O(5^7))                                                       # needs sage.rings.padics
            5^4 + 3*5^5 + 4*5^6 + 2*5^7 + 2*5^8 + O(5^9)
            sage: f.change_ring(Zp(5))(100 + O(5^7))                                    # needs sage.rings.padics
            5^4 + 3*5^5 + 4*5^6 + 2*5^7 + 2*5^8 + O(5^9)
            sage: f.change_ring(Zp(5))(100 + O(2^7))                                    # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute this value

        To substitute a value it must have valuation at least 1::

            sage: f(0)
            0
            sage: f(1 + t)
            Traceback (most recent call last):
            ...
            ValueError: Can only substitute elements of positive valuation
            sage: f(2 + O(5^3))                                                         # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: Can only substitute elements of positive valuation
            sage: f(t^-2)
            Traceback (most recent call last):
            ...
            ValueError: Can only substitute elements of positive valuation

        Unless, of course, it is being substituted in a series with infinite
        precision, i.e., a polynomial::

            sage: g = t^2 + t^3
            sage: g(1 + t + O(t^2))
            2 + 5*t + O(t^2)
            sage: g(3)
            36

        Arguments beyond the first can refer to the base ring::

            sage: P.<x> = GF(5)[]
            sage: Q.<y> = P[[]]
            sage: h = (1 - x*y)^-1 + O(y^7); h
            1 + x*y + x^2*y^2 + x^3*y^3 + x^4*y^4 + x^5*y^5 + x^6*y^6 + O(y^7)
            sage: h(y^2, 3)
            1 + 3*y^2 + 4*y^4 + 2*y^6 + y^8 + 3*y^10 + 4*y^12 + O(y^14)

        These secondary values can also be specified using keywords::

            sage: h(y=y^2, x=3)
            1 + 3*y^2 + 4*y^4 + 2*y^6 + y^8 + 3*y^10 + 4*y^12 + O(y^14)
            sage: h(y^2, x=3)
            1 + 3*y^2 + 4*y^4 + 2*y^6 + y^8 + 3*y^10 + 4*y^12 + O(y^14)
        """
        P = self.parent()

        if len(kwds) >= 1:
            name = P.variable_name()
            if name in kwds: # a keyword specifies the power series generator
                if x:
                    raise ValueError("must not specify %s keyword and positional argument" % name)
                a = self(kwds[name])
                del kwds[name]
                try:
                    return a(**kwds)
                except TypeError:
                    return a
            elif x:       # both keywords and positional arguments
                a = self(*x)
                try:
                    return a(**kwds)
                except TypeError:
                    return a
            else:                  # keywords but no positional arguments
                return P(self.__f(**kwds)).add_bigoh(self._prec)

        if len(x) == 0:
            return self

        if isinstance(x[0], tuple):
            x = x[0]
        a = x[0]

        s = self._prec
        if s == infinity:
            return self.__f(x)

        Q = a.parent()

        from sage.rings.padics.padic_generic import pAdicGeneric
        padic = isinstance(Q, pAdicGeneric)
        if padic:
            p = Q.prime()

        try:
            t = a.valuation()
        except (TypeError, AttributeError):
            if a.is_zero():
                t = infinity
            else:
                t = 0

        if t == infinity:
            return self[0]

        if t <= 0:
            raise ValueError("Can only substitute elements of positive valuation")

        if not Q.has_coerce_map_from(P.base_ring()):
            from sage.structure.element import canonical_coercion
            try:
                R = canonical_coercion(P.base_ring()(0), Q.base_ring()(0))[0].parent()
                self = self.change_ring(R)
            except TypeError:
                raise ValueError("Cannot substitute this value")

        r = (self - self[0]).valuation()
        if r == s:                 # self is constant + O(x^s)
            if padic:
                from sage.rings.big_oh import O
                return self[0] + O(p**(s*t))
            else:
                return P(self[0]).add_bigoh(s*t)

        try:
            u = a.prec()
        except AttributeError:
            u = a.precision_absolute()
        n = (s - r + 1)*t
        if n < u:
            a = a.add_bigoh(n)
            x = list(x)
            x[0] = a
            x = tuple(x)
        return self.__f(x)

    def _unsafe_mutate(self, i, value):
        """
        Sage assumes throughout that commutative ring elements are immutable.
        This is relevant for caching, etc.  But sometimes you need to change
        a power series and you really know what you're doing.  That's
        when this function is for you.

        ** DO NOT USE THIS ** unless you know what you're doing.

        EXAMPLES::

            sage: R.<t> = GF(7)[[]]
            sage: f = 3 + 6*t^3 + O(t^5)
            sage: f._unsafe_mutate(0, 5)
            sage: f
            5 + 6*t^3 + O(t^5)
            sage: f._unsafe_mutate(2, 1); f
            5 + t^2 + 6*t^3 + O(t^5)

        - Mutating can even bump up the precision::

            sage: f._unsafe_mutate(6, 1); f
            5 + t^2 + 6*t^3 + t^6 + O(t^7)
            sage: f._unsafe_mutate(0, 0); f
            t^2 + 6*t^3 + t^6 + O(t^7)
            sage: f._unsafe_mutate(1, 0); f
            t^2 + 6*t^3 + t^6 + O(t^7)
            sage: f._unsafe_mutate(11,0); f
            t^2 + 6*t^3 + t^6 + O(t^12)

            sage: g = t + O(t^7)
            sage: g._unsafe_mutate(1,0); g
            O(t^7)
        """
        self.__f._unsafe_mutate(i, value)
        self._prec = max(self._prec, i+1)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        This returns 0 for negative coefficients and raises an
        :exc:`IndexError` if trying to access beyond known coefficients.

        If ``n`` is a slice object ``[:k]``, this will return a power
        series of the same precision, whose coefficients are the same
        as ``self`` for those indices in the slice, and 0 otherwise.
        Other kinds of slicing are not allowed.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: f = 3/2 - 17/5*t^3 + O(t^5)
            sage: f[3]
            -17/5
            sage: f[-2]
            0
            sage: f[4]
            0
            sage: f[5]
            Traceback (most recent call last):
            ...
            IndexError: coefficient not known

        Using slices::

            sage: R.<t> = ZZ[[]]
            sage: f = (2-t)^5; f
            32 - 80*t + 80*t^2 - 40*t^3 + 10*t^4 - t^5
            sage: f[:4]
            32 - 80*t + 80*t^2 - 40*t^3
            sage: f = 1 + t^3 - 4*t^4 + O(t^7); f
            1 + t^3 - 4*t^4 + O(t^7)
            sage: f[:4]
            1 + t^3 + O(t^7)

        TESTS::

            sage: f[1:4]
            Traceback (most recent call last):
            ...
            IndexError: polynomial slicing with a start is not defined
        """
        if isinstance(n, slice):
            return PowerSeries_poly(self._parent, self.polynomial()[n],
                                    prec=self._prec, check=False)
        elif n < 0:
            return self.base_ring().zero()
        elif n > self.__f.degree():
            if self._prec > n:
                return self.base_ring().zero()
            else:
                raise IndexError("coefficient not known")
        return self.__f[n]

    def __iter__(self):
        """
        Return an iterator over the coefficients of this power series.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: f = t + 17/5*t^3 + 2*t^4 + O(t^5)
            sage: [a for a in f]
            [0, 1, 0, 17/5, 2]
        """
        return iter(self.__f)

    def __neg__(self):
        """
        Return the negative of this power series.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: f = t + 17/5*t^3 + 2*t^4 + O(t^5)
            sage: -f
            -t - 17/5*t^3 - 2*t^4 + O(t^5)
        """
        return PowerSeries_poly(self._parent, -self.__f,
                                         self._prec, check=False)

    cpdef _add_(self, right_m):
        """
        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: f = x^4 + O(x^5); f
            x^4 + O(x^5)
            sage: g = x^2 + O(x^3); g
            x^2 + O(x^3)
            sage: f+g
            x^2 + O(x^3)

        TESTS:

        In the past this could die with EXC_BAD_ACCESS (:issue:`8029`)::

            sage: # needs sage.rings.real_mpfr
            sage: A.<x> = RR['x']
            sage: B.<t> = PowerSeriesRing(A)
            sage: 1. + O(t)
            1.00000000000000 + O(t)
            sage: 1. + O(t^2)
            1.00000000000000 + O(t^2)
            sage: 1. + O(t^3)
            1.00000000000000 + O(t^3)
            sage: 1. + O(t^4)
            1.00000000000000 + O(t^4)
        """
        cdef PowerSeries_poly right = <PowerSeries_poly>right_m
        return PowerSeries_poly(self._parent, self.__f + right.__f,
                                self.common_prec_c(right), check=True)

    cpdef _sub_(self, right_m):
        """
        Return the difference of two power series.

        EXAMPLES::

            sage: k.<w> = ZZ[]
            sage: R.<t> = k[[]]
            sage: w*t^2 -w*t +13 - (w*t^2 + w*t)
            13 - 2*w*t
        """
        cdef PowerSeries_poly right = <PowerSeries_poly>right_m
        return PowerSeries_poly(self._parent, self.__f - right.__f,
                                self.common_prec_c(right), check=True)

    cpdef _mul_(self, right_r):
        """
        Return the product of two power series.

        EXAMPLES::

            sage: k.<w> = ZZ[[]]
            sage: (1+17*w+15*w^3+O(w^5))*(19*w^10+O(w^12))
            19*w^10 + 323*w^11 + O(w^12)
        """
        prec = self._mul_prec(right_r)
        return PowerSeries_poly(self._parent,
                                self.__f * (<PowerSeries_poly>right_r).__f,
                                prec=prec,
                                check=True)  # check, since truncation may be needed

    cpdef _rmul_(self, Element c):
        """
        Multiply ``self`` on the right by a scalar.

        EXAMPLES::

            sage: R.<t> = GF(7)[[]]
            sage: f = t + 3*t^4 + O(t^11)
            sage: f * GF(7)(3)
            3*t + 2*t^4 + O(t^11)
        """
        return PowerSeries_poly(self._parent, self.__f * c, self._prec, check=False)

    cpdef _lmul_(self, Element c):
        """
        Multiply ``self`` on the left by a scalar.

        EXAMPLES::

            sage: R.<t> = GF(11)[[]]
            sage: f = 1 + 3*t^4 + O(t^120)
            sage: 2 * f
            2 + 6*t^4 + O(t^120)
        """
        return PowerSeries_poly(self._parent, c * self.__f, self._prec, check=False)

    def __lshift__(PowerSeries_poly self, n):
        """
        Shift ``self`` to the left by ``n``, i.e. multiply by `x^n`.

        EXAMPLES::

            sage: R.<t> = QQ[[]]
            sage: f = 1 + t + t^4
            sage: f << 1
            t + t^2 + t^5
        """
        if n:
            return PowerSeries_poly(self._parent, self.__f << n, self._prec + n)
        else:
            return self

    def __rshift__(PowerSeries_poly self, n):
        """
        Shift ``self`` to the right by ``n``, i.e. multiply by `x^{-n}` and
        remove any terms of negative exponent.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<t> = GF(2)[[]]
            sage: f = t + t^4 + O(t^7)
            sage: f >> 1
            1 + t^3 + O(t^6)
            sage: f >> 10
            O(t^0)
        """
        if n:
            return PowerSeries_poly(self._parent, self.__f >> n, max(0,self._prec - n))
        else:
            return self

    def __invert__(self):
        """
        Return the inverse of the power series (i.e., a series `Y` such
        that `XY = 1`).

        The first nonzero coefficient must be a unit in
        the coefficient ring. If the valuation of the series is positive or
        `X` is not a unit, this function will return a
        :class:`sage.rings.laurent_series_ring_element.LaurentSeries`.

        EXAMPLES::

            sage: R.<q> = QQ[[]]
            sage: 1/(1+q + O(q**2))
            1 - q + O(q^2)
            sage: 1/(1+q)
            1 - q + q^2 - q^3 + q^4 - q^5 + q^6 - q^7 + q^8 - q^9 + q^10 - q^11 + q^12 - q^13 + q^14 - q^15 + q^16 - q^17 + q^18 - q^19 + O(q^20)
            sage: prec = R.default_prec(); prec
            20
            sage: 1/(1+q) + O(q^5)
            1 - q + q^2 - q^3 + q^4 + O(q^5)

        ::

            sage: 1/(q + q^2) + O(q^4)
            q^-1 - 1 + q - q^2 + q^3 + O(q^4)
            sage: g = 1/(q + q^2 + O(q^5))
            sage: g; g.parent()
            q^-1 - 1 + q - q^2 + O(q^3)
            Laurent Series Ring in q over Rational Field

        ::

            sage: 1/g
            q + q^2 + O(q^5)
            sage: (1/g).parent()
            Laurent Series Ring in q over Rational Field

        ::

            sage: 1/(2 + q) + O(q^5)
            1/2 - 1/4*q + 1/8*q^2 - 1/16*q^3 + 1/32*q^4 + O(q^5)

        ::

            sage: R.<q> = PowerSeriesRing(QQ, name='q', default_prec=5)
            sage: f = 1 + q + q^2 + O(q^50)
            sage: f/10
            1/10 + 1/10*q + 1/10*q^2 + O(q^50)
            sage: f/(10+q)
            1/10 + 9/100*q + 91/1000*q^2 - 91/10000*q^3 + 91/100000*q^4 + O(q^5)

        ::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: u = 17 + 3*t^2 + 19*t^10 + O(t^12)
            sage: v = ~u; v
            1/17 - 3/289*t^2 + 9/4913*t^4 - 27/83521*t^6 + 81/1419857*t^8 - 1587142/24137569*t^10 + O(t^12)
            sage: u*v
            1 + O(t^12)

        If we try a nonzero, non-unit constant term, we end up in
        the fraction field, i.e. the Laurent series ring::

            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: ~R(2)
            1/2
            sage: parent(~R(2))
            Laurent Series Ring in t over Rational Field

        As for units, we stay in the power series ring::

            sage: ~R(-1)
            -1
            sage: parent(~R(-1))
            Power Series Ring in t over Integer Ring

        However, inversion of non-unit elements must fail when the underlying
        ring is not an integral domain::

            sage: R = IntegerModRing(8)
            sage: P.<s> = R[[]]
            sage: ~P(2)
            Traceback (most recent call last):
            ...
            ValueError: must be an integral domain
        """
        if self.is_one():
            return self
        prec = self.prec()
        if prec is infinity:
            if self.degree() > 0:
                prec = self._parent.default_prec()
            else:
                # constant series
                a = self[0]
                if not a.is_unit():
                    from sage.categories.integral_domains import IntegralDomains
                    if self._parent in IntegralDomains():
                        R = self._parent.fraction_field()
                        return 1 / R(a)
                    else:
                        raise ValueError('must be an integral domain')
                try:
                    a = a.inverse_unit()
                except (AttributeError, NotImplementedError):
                    a = self._parent.base_ring()(~a)
                return self._parent(a, prec=infinity)

        if self.valuation() > 0:
            u = ~self.valuation_zero_part()    # inverse of unit part
            R = self._parent.laurent_series_ring()
            return R(u, -self.valuation())

        return self._parent(self.truncate().inverse_series_trunc(prec), prec=prec)

    def truncate(self, prec=infinity):
        """
        The polynomial obtained from power series by truncation at
        precision ``prec``.

        EXAMPLES::

            sage: R.<I> = GF(2)[[]]
            sage: f = 1/(1+I+O(I^8)); f
            1 + I + I^2 + I^3 + I^4 + I^5 + I^6 + I^7 + O(I^8)
            sage: f.truncate(5)
            I^4 + I^3 + I^2 + I + 1
        """
        if prec is infinity:
            return self.__f
        else:
            return self.__f.truncate(prec)

    cdef _inplace_truncate(self, long prec):
        """
        Truncate ``self`` to precision ``prec`` in place.

        .. NOTE::

            This is very unsafe, since power series are supposed to
            be immutable in Sage. Use at your own risk!
        """
        self.__f = self.__f._inplace_truncate(prec)
        self.prec = prec
        return self

    def truncate_powerseries(self, long prec):
        r"""
        Given input ``prec`` = `n`, returns the power series of degree
        `< n` which is equivalent to ``self`` modulo `x^n`.

        EXAMPLES::

            sage: R.<I> = GF(2)[[]]
            sage: f = 1/(1+I+O(I^8)); f
            1 + I + I^2 + I^3 + I^4 + I^5 + I^6 + I^7 + O(I^8)
            sage: f.truncate_powerseries(5)
            1 + I + I^2 + I^3 + I^4 + O(I^5)
        """
        return PowerSeries_poly(self._parent, self.__f.truncate(prec),
                                min(self._prec, prec), check=False)

    def list(self):
        """
        Return the list of known coefficients for ``self``.

        This is just the list of coefficients of the underlying
        polynomial, so in particular, need not have length equal to
        ``self.prec()``.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: f = 1 - 5*t^3 + t^5 + O(t^7)
            sage: f.list()
            [1, 0, 0, -5, 0, 1]
        """
        return self.__f.list()

    def monomial_coefficients(self):
        """
        Return a dictionary of coefficients for ``self``.

        This is simply a dict for the underlying polynomial, so need
        not have keys corresponding to every number smaller than
        ``self.prec()``.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: f = 1 + t^10 + O(t^12)
            sage: f.monomial_coefficients()
            {0: 1, 10: 1}

        ``dict`` is an alias::

            sage: f.dict()
            {0: 1, 10: 1}
        """
        return self.__f.monomial_coefficients()

    dict = monomial_coefficients

    def _derivative(self, var=None):
        """
        Return the derivative of this power series with respect
        to the variable ``var``.

        If ``var`` is ``None`` or is the generator of this ring, we
        take the derivative with respect to the generator.

        Otherwise, we call ``_derivative(var)`` on each coefficient of
        the series.

        .. SEEALSO::

            ``self.derivative()``

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = 2 + 3*t^2 + t^100000 + O(t^10000000); f
            2 + 3*t^2 + t^100000 + O(t^10000000)
            sage: f._derivative()
            6*t + 100000*t^99999 + O(t^9999999)
            sage: f._derivative(t)
            6*t + 100000*t^99999 + O(t^9999999)

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<y> = PowerSeriesRing(R, sparse=True)
            sage: f = x^3*y^4 + O(y^5)
            sage: f._derivative()
            4*x^3*y^3 + O(y^4)
            sage: f._derivative(y)
            4*x^3*y^3 + O(y^4)
            sage: f._derivative(x)
            3*x^2*y^4 + O(y^5)

        TESTS::

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: x = var('x')                                                          # needs sage.symbolic
            sage: t.derivative(x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to x
        """
        if var is not None and var != self._parent.gen():
            try:
                # call _derivative() recursively on coefficients
                return PowerSeries_poly(self._parent, self.__f._derivative(var),
                                    self.prec(), check=False)
            except AttributeError:
                raise ValueError('cannot differentiate with respect to {}'.format(var))

        # compute formal derivative with respect to generator
        return PowerSeries_poly(self._parent, self.__f._derivative(),
                                self.prec()-1, check=False)

    def integral(self,var=None):
        """
        Return the integral of this power series.

        By default, the integration variable is the variable of the
        power series.

        Otherwise, the integration variable is the optional parameter ``var``.

        .. NOTE::

            The integral is always chosen so the constant term is 0.

        EXAMPLES::

            sage: k.<w> = QQ[[]]
            sage: (1+17*w+15*w^3+O(w^5)).integral()
            w + 17/2*w^2 + 15/4*w^4 + O(w^6)
            sage: (w^3 + 4*w^4 + O(w^7)).integral()
            1/4*w^4 + 4/5*w^5 + O(w^8)
            sage: (3*w^2).integral()
            w^3

        TESTS::

            sage: t = PowerSeriesRing(QQ,'t').gen()
            sage: f = t + 5*t^2 + 21*t^3
            sage: g = f.integral(); g
            1/2*t^2 + 5/3*t^3 + 21/4*t^4
            sage: g.parent()
            Power Series Ring in t over Rational Field

            sage: R.<x> = QQ[]
            sage: t = PowerSeriesRing(R,'t').gen()
            sage: f = x*t +5*t^2
            sage: f.integral()
            1/2*x*t^2 + 5/3*t^3
            sage: f.integral(x)
            1/2*x^2*t + 5*x*t^2
        """
        return PowerSeries_poly(self._parent, self.__f.integral(var),
                                self.prec()+1, check=False)

    def reverse(self, precision=None):
        """
        Return the reverse of `f`, i.e., the series `g` such that `g(f(x)) = x`.

        Given an optional argument ``precision``, return the reverse with given
        precision (note that the reverse can have precision at most
        ``f.prec()``).  If `f` has infinite precision, and the argument
        ``precision`` is not given, then the precision of the reverse defaults
        to the default precision of ``f.parent()``.

        Note that this is only possible if the valuation of ``self`` is exactly
        1.

        ALGORITHM:

        We first attempt to pass the computation to pari; if this fails, we
        use Lagrange inversion.  Using ``sage: set_verbose(1)`` will print
        a message if passing to pari fails.

        If the base ring has positive characteristic, then we attempt to
        lift to a characteristic zero ring and perform the reverse there.
        If this fails, an error is raised.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(QQ)
            sage: f = 2*x + 3*x^2 - x^4 + O(x^5)
            sage: g = f.reverse()
            sage: g
            1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)
            sage: f(g)
            x + O(x^5)
            sage: g(f)
            x + O(x^5)

            sage: A.<t> = PowerSeriesRing(ZZ)
            sage: a = t - t^2 - 2*t^4 + t^5 + O(t^6)
            sage: b = a.reverse(); b
            t + t^2 + 2*t^3 + 7*t^4 + 25*t^5 + O(t^6)
            sage: a(b)
            t + O(t^6)
            sage: b(a)
            t + O(t^6)

            sage: B.<b,c> = PolynomialRing(ZZ)
            sage: A.<t> = PowerSeriesRing(B)
            sage: f = t + b*t^2 + c*t^3 + O(t^4)
            sage: g = f.reverse(); g
            t - b*t^2 + (2*b^2 - c)*t^3 + O(t^4)
            sage: f(g)
            t + O(t^4)
            sage: g(f)
            t + O(t^4)

            sage: A.<t> = PowerSeriesRing(ZZ)
            sage: B.<s> = A[[]]
            sage: f = (1 - 3*t + 4*t^3 + O(t^4))*s + (2 + t + t^2 + O(t^3))*s^2 + O(s^3)
            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(1)
            sage: g = f.reverse(); g
            verbose 1 (<module>) passing to pari failed; trying Lagrange inversion
            (1 + 3*t + 9*t^2 + 23*t^3 + O(t^4))*s + (-2 - 19*t - 118*t^2 + O(t^3))*s^2 + O(s^3)
            sage: set_verbose(0)
            sage: f(g) == g(f) == s
            True

        If the leading coefficient is not a unit, we pass to its fraction
        field if possible::

            sage: A.<t> = PowerSeriesRing(ZZ)
            sage: a = 2*t - 4*t^2 + t^4 - t^5 + O(t^6)
            sage: a.reverse()
            1/2*t + 1/2*t^2 + t^3 + 79/32*t^4 + 437/64*t^5 + O(t^6)

            sage: B.<b> = PolynomialRing(ZZ)
            sage: A.<t> = PowerSeriesRing(B)
            sage: f = 2*b*t + b*t^2 + 3*b^2*t^3 + O(t^4)
            sage: g = f.reverse(); g
            1/(2*b)*t - 1/(8*b^2)*t^2 + ((-3*b + 1)/(16*b^3))*t^3 + O(t^4)
            sage: f(g)
            t + O(t^4)
            sage: g(f)
            t + O(t^4)

        We can handle some base rings of positive characteristic::

            sage: A8.<t> = PowerSeriesRing(Zmod(8))
            sage: a = t - 15*t^2 - 2*t^4 + t^5 + O(t^6)
            sage: b = a.reverse(); b
            t + 7*t^2 + 2*t^3 + 5*t^4 + t^5 + O(t^6)
            sage: a(b)
            t + O(t^6)
            sage: b(a)
            t + O(t^6)

        The optional argument ``precision`` sets the precision of the output::

            sage: R.<x> = PowerSeriesRing(QQ)
            sage: f = 2*x + 3*x^2 - 7*x^3 + x^4 + O(x^5)
            sage: g = f.reverse(precision=3); g
            1/2*x - 3/8*x^2 + O(x^3)
            sage: f(g)
            x + O(x^3)
            sage: g(f)
            x + O(x^3)

        If the input series has infinite precision, the precision of the
        output is automatically set to the default precision of the parent
        ring::

            sage: R.<x> = PowerSeriesRing(QQ, default_prec=20)
            sage: (x - x^2).reverse() # get some Catalan numbers
            x + x^2 + 2*x^3 + 5*x^4 + 14*x^5 + 42*x^6 + 132*x^7 + 429*x^8 + 1430*x^9
             + 4862*x^10 + 16796*x^11 + 58786*x^12 + 208012*x^13 + 742900*x^14
             + 2674440*x^15 + 9694845*x^16 + 35357670*x^17 + 129644790*x^18
             + 477638700*x^19 + O(x^20)
            sage: (x - x^2).reverse(precision=3)
            x + x^2 + O(x^3)

        TESTS::

            sage: R.<x> = PowerSeriesRing(QQ)
            sage: f = 1 + 2*x + 3*x^2 - x^4 + O(x^5)
            sage: f.reverse()
            Traceback (most recent call last):
            ...
            ValueError: Series must have valuation one for reversion.

            sage: Series = PowerSeriesRing(SR, 'x')                                     # needs sage.symbolic
            sage: ser = Series([0, pi]); ser                                            # needs sage.symbolic
            pi*x
            sage: ser.reverse()                                                         # needs sage.symbolic
            1/pi*x + O(x^20)
        """
        if self.valuation() != 1:
            raise ValueError("Series must have valuation one for reversion.")

        f = self

        if f.prec() is infinity and precision is None:
            precision = f.parent().default_prec()
        if precision:
            f = f.add_bigoh(precision)

        out_prec = f.prec()

        if not f[1].is_unit():
            # if leading coefficient is not a unit, attempt passing
            # to fraction field
            try:
                f = f.change_ring(f.base_ring().fraction_field())
            except TypeError:
                raise TypeError("Leading coefficient must be a unit, or base ring must have a fraction field.")

        # set output parent after possibly passing to fraction field,
        # but before possibly lifting to characteristic zero
        out_parent = f.parent()

        # first, try reversion with pari; this is faster than Lagrange inversion
        try:
            f2 = f.__pari__()
            g = f2.serreverse()
            return PowerSeries_poly(f.parent(), g.Vec(-out_prec), out_prec)
        except (TypeError, ValueError, AttributeError, PariError, ImportError):
            # if pari fails, continue with Lagrange inversion
            from sage.misc.verbose import verbose
            verbose("passing to pari failed; trying Lagrange inversion")

        if f.parent().characteristic():
            # over a ring of positive characteristic, attempt lifting to
            # characteristic zero ring
            verbose("parent ring has positive characteristic; attempting lift to characteristic zero")
            base_lift = f.base_ring().lift().codomain()
            verbose("characteristic zero base is "+str(base_lift))
            f_lift = f.change_ring(base_lift)
            verbose("f_lift is "+str(f_lift))
            rev_lift = f_lift.reverse()
            return rev_lift.change_ring(f.base_ring())

        t = f.parent().gen()
        R = f.parent().base_ring()

        h = t/f
        k = 1
        g = 0
        for i in range(1, out_prec):
            k *= h
            g += R(k.padded_list(i)[i - 1]/i)*t**i
        g = g.add_bigoh(out_prec)
        return PowerSeries_poly(out_parent, g, out_prec, check=False)

    def pade(self, m, n):
        r"""
        Return the Padé approximant of ``self`` of index `(m, n)`.

        The Padé approximant of index `(m, n)` of a formal power
        series `f` is the quotient `Q/P` of two polynomials `Q` and `P`
        such that `\deg(Q)\leq m`, `\deg(P)\leq n` and

        .. MATH::

            f(z) - Q(z)/P(z) = O(z^{m+n+1}).

        The formal power series `f` must be known up to order `n + m`.

        See :wikipedia:`Padé\_approximant`

        INPUT:

        - ``m``, ``n`` -- integers, describing the degrees of the polynomials

        OUTPUT: a ratio of two polynomials

        ALGORITHM:

        This method uses the formula as a quotient of two determinants.

        .. SEEALSO::

            * :mod:`sage.matrix.berlekamp_massey`,
            * :meth:`sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint.rational_reconstruction`

        EXAMPLES::

            sage: z = PowerSeriesRing(QQ, 'z').gen()
            sage: exp(z).pade(4, 0)
            1/24*z^4 + 1/6*z^3 + 1/2*z^2 + z + 1
            sage: exp(z).pade(1, 1)
            (-z - 2)/(z - 2)
            sage: exp(z).pade(3, 3)
            (-z^3 - 12*z^2 - 60*z - 120)/(z^3 - 12*z^2 + 60*z - 120)
            sage: log(1-z).pade(4, 4)
            (25/6*z^4 - 130/3*z^3 + 105*z^2 - 70*z)/(z^4 - 20*z^3 + 90*z^2
            - 140*z + 70)
            sage: sqrt(1+z).pade(3, 2)
            (1/6*z^3 + 3*z^2 + 8*z + 16/3)/(z^2 + 16/3*z + 16/3)
            sage: exp(2*z).pade(3, 3)
            (-z^3 - 6*z^2 - 15*z - 15)/(z^3 - 6*z^2 + 15*z - 15)

        TESTS:

        With real coefficients::

            sage: # needs sage.rings.real_mpfr
            sage: R.<z> = RR[[]]
            sage: f = exp(2*z)
            sage: f.pade(3, 3) # abs tol 1e-10
            (-z^3 - 6.0*z^2 - 15.0*z - 15.0)/(z^3 - 6.0*z^2 + 15.0*z - 15.0)

        When precision is too low::

            sage: # needs sage.rings.real_mpfr
            sage: f = z + O(z**6)
            sage: f.pade(4, 4)
            Traceback (most recent call last):
            ...
            ValueError: the precision of the series is not large enough

        Check that :issue:`21212` is fixed::

            sage: QQx.<x> = QQ[[]]
            sage: (1 + x + O(x^100)).pade(2,2)
            x + 1

        Check for correct precision::

            sage: QQx.<x> = QQ[[]]
            sage: (1 + x + O(x^2)).pade(0,1)
            -1/(x - 1)
        """
        if self.precision_absolute() < n + m + 1:
            raise ValueError("the precision of the series is not large enough")
        polyring = self.parent()._poly_ring()
        z = polyring.gen()
        c = self.polynomial()
        u, v = c.rational_reconstruction(z**(n + m + 1), m, n)
        return u / v

    def _symbolic_(self, ring):
        """
        Conversion to symbolic series.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: R.<x> = PowerSeriesRing(QQ)
            sage: s = R([1,2,3,4,5], prec=10); s
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + O(x^10)
            sage: SR(s)
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + Order(x^10)
            sage: SR(s).is_terminating_series()
            False
            sage: SR(s).variables()
            (x,)
            sage: s = R([1,2,3,4,5]); s
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4
            sage: SR(s)
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4
            sage: _.is_terminating_series()
            True

        TESTS:

        Check that :issue:`18094` is fixed::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: SR(R(0).add_bigoh(20))                                                # needs sage.symbolic
            Order(x^20)
        """
        from sage.symbolic.ring import SR
        pex = SR(self.polynomial())
        var = SR.var(self.variable())
        return pex.series(var, self.prec())


def make_powerseries_poly_v0(parent, f, prec, is_gen):
    """
    Return the power series specified by ``f``, ``prec``, and ``is_gen``.

    This function exists for the purposes of pickling. Do not delete
    this function -- if you change the internal representation,
    instead make a new function and make sure that both kinds of
    objects correctly unpickle as the new type.

    EXAMPLES::

        sage: R.<t> = QQ[[]]
        sage: sage.rings.power_series_poly.make_powerseries_poly_v0(R, t, infinity, True)
        t
    """
    return PowerSeries_poly(parent, f, prec, 0, is_gen)


cdef class BaseRingFloorDivAction(Action):
    """
    The floor division action of the base ring on a formal power series.
    """
    cpdef _act_(self, g, x):
        r"""
        Let ``g`` act on ``x`` under ``self``.

        Regardless of whether this is a left or right action, the acting
        element comes first.

        INPUT:

        - ``g`` -- an object with parent ``self.G``
        - ``x`` -- an object with parent ``self.US()``

        .. WARNING::

            This is meant to be a fast internal function, so the
            conditions on the input are not checked!

        EXAMPLES:

        One gets the correct parent with floor division::

            sage: A = ZZ[['t']]
            sage: f = A([3*2**n for n in range(6)]).O(6)
            sage: g = f // 3; g
            1 + 2*t + 4*t^2 + 8*t^3 + 16*t^4 + 32*t^5 + O(t^6)
            sage: g.parent()
            Power Series Ring in t over Integer Ring

        whereas the parent is larger with division::

            sage: parent(f/3)
            Power Series Ring in t over Rational Field

        Floor division in case that the power series is not divisible by the divisor::

            sage: f = A([2**n for n in range(6)]).O(6)
            sage: g = f // 3; g
            t^2 + 2*t^3 + 5*t^4 + 10*t^5 + O(t^6)

        Another example::

            sage: s = polygen(QQ,'s')
            sage: A = s.parent()[['t']]
            sage: f = A([(s+2)*(s+n) for n in range(5)]).O(5)
            sage: g = f // (s + 2); g
            s + (s + 1)*t + (s + 2)*t^2 + (s + 3)*t^3 + (s + 4)*t^4 + O(t^5)
            sage: g.parent()
            Power Series Ring in t over Univariate Polynomial Ring in s
            over Rational Field

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: t // 2
            1/2*t
        """
        cdef PowerSeries_poly elt = <PowerSeries_poly> x
        prec = x.prec()
        P = self.US()
        g = P.base_ring()(g)
        return type(x)(P, elt.__f // g, prec=prec, check=False)
