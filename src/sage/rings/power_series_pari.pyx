r"""
Power series implemented using PARI

EXAMPLES:

This implementation can be selected for any base ring supported by
PARI by passing the keyword ``implementation='pari'`` to the
:func:`~sage.rings.power_series_ring.PowerSeriesRing` constructor::

    sage: R.<q> = PowerSeriesRing(ZZ, implementation='pari'); R
    Power Series Ring in q over Integer Ring
    sage: S.<t> = PowerSeriesRing(CC, implementation='pari'); S                         # needs sage.rings.real_mpfr
    Power Series Ring in t over Complex Field with 53 bits of precision

Note that only the type of the elements depends on the implementation,
not the type of the parents::

    sage: type(R)
    <class 'sage.rings.power_series_ring.PowerSeriesRing_domain_with_category'>
    sage: type(q)
    <class 'sage.rings.power_series_pari.PowerSeries_pari'>
    sage: type(S)                                                                       # needs sage.rings.real_mpfr
    <class 'sage.rings.power_series_ring.PowerSeriesRing_over_field_with_category'>
    sage: type(t)                                                                       # needs sage.rings.real_mpfr
    <class 'sage.rings.power_series_pari.PowerSeries_pari'>

If `k` is a finite field implemented using PARI, this is the default
implementation for power series over `k`::

    sage: k.<c> = GF(5^12)
    sage: type(c)
    <class 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt'>
    sage: A.<x> = k[[]]
    sage: type(x)
    <class 'sage.rings.power_series_pari.PowerSeries_pari'>

.. WARNING::

    Because this implementation uses the PARI interface, the PARI variable
    ordering must be respected in the sense that the variable name of the
    power series ring must have higher priority than any variable names
    occurring in the base ring::

        sage: R.<y> = QQ[]
        sage: S.<x> = PowerSeriesRing(R, implementation='pari'); S
        Power Series Ring in x over Univariate Polynomial Ring in y over Rational Field

    Reversing the variable ordering leads to errors::

        sage: R.<x> = QQ[]
        sage: S.<y> = PowerSeriesRing(R, implementation='pari')
        Traceback (most recent call last):
        ...
        PariError: incorrect priority in gtopoly: variable x <= y

AUTHORS:

- Peter Bruin (December 2013): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013-2017 Peter Bruin <P.J.Bruin@math.leidenuniv.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cypari2.gen cimport Gen as pari_gen
from cypari2.pari_instance cimport get_var
from cypari2.paridecl cimport gel, typ, lg, valp, varn, t_POL, t_SER, t_RFRAC, t_VEC
from sage.libs.pari import pari

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.infinity import infinity


cdef PowerSeries_pari construct_from_pari(parent, pari_gen g):
    r"""
    Fast construction of power series from PARI objects of suitable
    type (series, polynomials, scalars and rational functions).

    The resulting series has the same precision as `g`, unless `g` is
    a rational function, in which case the default precision of
    ``parent`` is used.

    TESTS:

    Check for :issue:`33224`::

        sage: R.<z> = LaurentSeriesRing(QQ, implementation='pari')
        sage: (z^-2).prec()
        +Infinity
    """
    cdef long t = typ(g.g)
    v = parent.variable_name()
    if t == t_SER and varn(g.g) == get_var(v):
        prec = lg(g.g) - 2 + valp(g.g)
    elif t == t_RFRAC and pari.denominator(g) != 1:
        prec = parent.default_prec()
        g = g.Ser(v, prec - g.valuation(v))
    else:
        prec = infinity
    cdef PowerSeries_pari x = PowerSeries_pari.__new__(PowerSeries_pari)
    x._parent = parent
    x._prec = prec
    x.g = g
    return x


cdef class PowerSeries_pari(PowerSeries):
    r"""
    A power series implemented using PARI.

    INPUT:

    - ``parent`` -- the power series ring to use as the parent

    - ``f`` -- object from which to construct a power series

    - ``prec`` -- (default: infinity) precision of the element
      to be constructed

    - ``check`` -- ignored, but accepted for compatibility with
      :class:`~sage.rings.power_series_poly.PowerSeries_poly`
    """
    def __init__(self, parent, f=0, prec=infinity, check=None):
        """
        Initialize ``self``.

        TESTS::

            sage: # needs sage.rings.real_mpfr
            sage: R.<q> = PowerSeriesRing(CC, implementation='pari')
            sage: TestSuite(q).run()
            sage: f = q - q^3 + O(q^10)
            sage: TestSuite(f).run()
        """
        cdef Parent f_parent
        cdef pari_gen g
        cdef long t
        v = parent.variable_name()
        R = parent.base_ring()
        P = parent._poly_ring()

        if isinstance(f, PowerSeries):  # not only PowerSeries_pari
            f_parent = (<PowerSeries>f)._parent
            if f_parent is parent:
                if prec is infinity:
                    prec = (<PowerSeries>f)._prec
                g = f.__pari__()
            elif R.has_coerce_map_from(f_parent):
                g = R.coerce(f).__pari__()
            else:
                if prec is infinity:
                    prec = f.prec()
                g = f.polynomial().change_ring(R).__pari__()
        elif isinstance(f, Polynomial):
            f_parent = (<Polynomial>f)._parent
            if f_parent is P:
                g = f.__pari__()
            elif R.has_coerce_map_from(f_parent):
                g = R.coerce(f).__pari__()
            else:
                g = P.coerce(f).__pari__()
        elif isinstance(f, pari_gen):
            g = f
            t = typ(g.g)
            if t == t_POL:
                g = P(g).__pari__()
            elif t == t_SER and varn(g.g) == get_var(v):
                if valp(g.g) < 0:
                    raise ValueError('series has negative valuation')
                if prec is infinity:
                    prec = lg(g.g) - 2 + valp(g.g)
                g = P(g.Pol(v)).__pari__()
            elif t == t_RFRAC:
                if prec is infinity:
                    prec = parent.default_prec()
                g = P.fraction_field()(g).__pari__()
                g = g.Ser(v, prec - g.valuation(v))
            elif t == t_VEC:
                g = P(g.Polrev(v)).__pari__()
            else:
                g = R(g).__pari__()
        elif isinstance(f, (list, tuple)):
            g = pari([R.coerce(x) for x in f]).Polrev(v)
        else:
            g = R.coerce(f).__pari__()

        if prec is infinity:
            self.g = g
        else:
            if not g:
                self.g = g.Ser(v, prec)
            else:
                self.g = g.Ser(v, prec - g.valuation(v))

        PowerSeries.__init__(self, parent, prec)

    def __hash__(self):
        """
        Return a hash of ``self``.

        TESTS::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: hash(t^2 + 1) == hash(pari(t^2 + 1))
            True
        """
        return hash(self.g)

    def __reduce__(self):
        """
        Used for pickling.

        EXAMPLES::

            sage: A.<z> = PowerSeriesRing(RR, implementation='pari')
            sage: f = z - z^3 + O(z^10)
            sage: z == loads(dumps(z))
            True
            sage: f == loads(dumps(f))
            True
        """
        return PowerSeries_pari, (self._parent, self.g, self._prec, False)

    def __pari__(self):
        """
        Convert ``self`` to a PARI object.

        TESTS::

            sage: R.<t> = PowerSeriesRing(GF(7), implementation='pari')
            sage: (3 - t^3 + O(t^5)).__pari__()
            Mod(3, 7) + Mod(6, 7)*t^3 + O(t^5)
        """
        return self.g

    def polynomial(self):
        """
        Convert ``self`` to a polynomial.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(GF(7), implementation='pari')
            sage: f = 3 - t^3 + O(t^5)
            sage: f.polynomial()
            6*t^3 + 3
        """
        return self._parent._poly_ring()(self.list())

    def valuation(self):
        """
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: (5 - t^8 + O(t^11)).valuation()
            0
            sage: (-t^8 + O(t^11)).valuation()
            8
            sage: O(t^7).valuation()
            7
            sage: R(0).valuation()
            +Infinity
        """
        if not self.g:
            return self._prec
        return self.g.valuation(self._parent.variable_name())

    def __bool__(self):
        """
        Return ``True`` if ``self`` is nonzero, and ``False`` otherwise.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(GF(11), implementation='pari')
            sage: bool(1 + t + O(t^18))
            True
            sage: bool(R(0))
            False
            sage: bool(O(t^18))
            False
        """
        return bool(self.g)

    def __call__(self, *x, **kwds):
        """
        Evaluate ``self`` at `x = a`.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = t^2 + t^3 + O(t^6)
            sage: f(t^3)
            t^6 + t^9 + O(t^18)
            sage: f(t=t^3)
            t^6 + t^9 + O(t^18)
            sage: f(f)
            t^4 + 2*t^5 + 2*t^6 + 3*t^7 + O(t^8)
            sage: f(f)(f) == f(f(f))
            True

        The following demonstrates that the problems raised in
        :issue:`3979` and :issue:`5367` are solved::

            sage: [f(t^2 + O(t^n)) for n in [9, 10, 11]]
            [t^4 + t^6 + O(t^11), t^4 + t^6 + O(t^12), t^4 + t^6 + O(t^12)]
            sage: f(t^2)
            t^4 + t^6 + O(t^12)

        It is possible to substitute a series for which only
        the precision is defined::

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

            sage: S.<u> = PowerSeriesRing(GF(7), implementation='pari')
            sage: f(2*u + u^3 + O(u^5))
            4*u^2 + u^3 + 4*u^4 + 5*u^5 + O(u^6)

        Substituting `p`-adic numbers::

            sage: # needs sage.rings.padics
            sage: f(100 + O(5^7))
            5^4 + 3*5^5 + 4*5^6 + 2*5^7 + 2*5^8 + O(5^9)
            sage: ff = PowerSeriesRing(pAdicRing(5), 't', implementation='pari')(f)
            sage: ff
            (1 + O(5^20))*t^2 + (1 + O(5^20))*t^3 + O(t^6)
            sage: ff(100 + O(5^7))
            5^4 + 3*5^5 + 4*5^6 + 2*5^7 + 2*5^8 + O(5^9)
            sage: ff(100 + O(2^7))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
             '5-adic Ring with capped relative precision 20' and
             '2-adic Ring with capped relative precision 20'

        The argument must have valuation at least 1, unless the series
        is actually a polynomial::

            sage: f(0)
            0
            sage: f(1 + t)
            Traceback (most recent call last):
            ...
            ValueError: can only substitute elements of positive valuation
            sage: f(t^-2)
            Traceback (most recent call last):
            ...
            ValueError: can only substitute elements of positive valuation
            sage: f(2 + O(5^3))                                                         # needs sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: can only substitute elements of positive valuation
            sage: g = t^2 + t^3
            sage: g(1 + t + O(t^2))
            2 + 5*t + O(t^2)
            sage: g(3)
            36

        Substitution of variables belonging to the base ring can be
        done using keywords::

            sage: P.<a> = GF(5)[]
            sage: Q.<x> = PowerSeriesRing(P, implementation='pari')
            sage: h = (1 - a*x)^-1 + O(x^7); h
            1 + a*x + a^2*x^2 + a^3*x^3 + a^4*x^4 + a^5*x^5 + a^6*x^6 + O(x^7)
            sage: h(x^2, a=3)
            1 + 3*x^2 + 4*x^4 + 2*x^6 + x^8 + 3*x^10 + 4*x^12 + O(x^14)
        """
        if len(kwds) >= 1:
            name = self._parent.variable_name()
            if name in kwds:  # the series variable is specified by a keyword
                if len(x):
                    raise ValueError("must not specify %s keyword and positional argument" % name)
                x = [kwds[name]]
                del kwds[name]

        if len(x) != 1:
            raise ValueError("must specify exactly one positional argument")

        a = x[0]

        s = self._prec
        if s is infinity:
            return self.polynomial()(a)

        # Determine the parent of the result.
        P = self._parent
        Q = a.parent()
        if not Q.has_coerce_map_from(P.base_ring()):
            from sage.structure.element import canonical_coercion
            a = canonical_coercion(P.base_ring()(0), a)[1]
            Q = a.parent()

        # The result is defined if the ring Q is complete with respect
        # to an ideal I, and the element a lies in I.  Here we only
        # implement a few special cases.
        from sage.rings.padics.padic_generic import pAdicGeneric
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        from sage.rings.laurent_series_ring import LaurentSeriesRing
        if isinstance(Q, pAdicGeneric):
            # Substitution of p-adic numbers in power series is
            # currently not implemented in PARI (2.8.0-development).
            t = a.valuation()
            if t <= 0:
                raise ValueError("can only substitute elements of positive valuation")
            return Q(self.polynomial()(a)).add_bigoh(t * self._prec)
        elif isinstance(Q, (PowerSeriesRing_generic, LaurentSeriesRing)):
            # In Sage, we want an error to be raised when trying to
            # substitute a series of nonpositive valuation, but PARI
            # (2.8.0-development) does not do this.  For example,
            # subst(1 + O(x), x, 1/y) yields O(y^-1).
            if a.valuation() <= 0:
                raise ValueError("can only substitute elements of positive valuation")
        elif isinstance(Q, PolynomialRing_generic):
            Q = Q.completion(Q.gen())
        elif Q.is_exact() and not a:
            pass
        else:
            raise ValueError('cannot substitute %s in %s' % (a, self))

        if not kwds:
            return Q(self.g(a))
        else:
            kwds[P.variable_name()] = a
            return Q(self.g(**kwds))

    def __getitem__(self, n):
        r"""
        Return the ``n``-th coefficient of ``self``.

        If ``n`` is a slice object, this returns a power series of the
        same precision, whose coefficients are the same as ``self``
        for those indices in the slice, and 0 otherwise.

        Returns 0 for negative coefficients.  Raises an :exc:`IndexError`
        if trying to access beyond known coefficients.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
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
            IndexError: index out of range

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = (2-t)^5; f
            32 - 80*t + 80*t^2 - 40*t^3 + 10*t^4 - t^5
            sage: f[:4]
            32 - 80*t + 80*t^2 - 40*t^3

            sage: f = 1 + t^3 - 4*t^4 + O(t^7); f
            1 + t^3 - 4*t^4 + O(t^7)
            sage: f[:4]
            1 + t^3 + O(t^7)
        """
        cdef long t
        if isinstance(n, slice):
            return PowerSeries_pari(self._parent, self.polynomial()[n],
                                    prec=self._prec)
        if n < 0:
            return self.base_ring().zero()

        t = typ(self.g.g)
        if t == t_POL or t == t_SER:
            h = self.g[n]
        else:
            h = self.g
        return self.base_ring()(h)

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        TESTS::

            sage: R.<t> = PowerSeriesRing(QQ, default_prec=6, implementation='pari')
            sage: ~(R(1-t))
            1 + t + t^2 + t^3 + t^4 + t^5 + O(t^6)
        """
        h = ~self.g
        if h.valuation(self._parent.variable_name()) < 0:
            return self._parent.laurent_series_ring()(h)
        return construct_from_pari(self._parent, h)

    def __neg__(self):
        """
        Return the negative of ``self``.

        TESTS::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = t + 17/5*t^3 + 2*t^4 + O(t^5)
            sage: -f
            -t - 17/5*t^3 - 2*t^4 + O(t^5)
        """
        return construct_from_pari(self._parent, -self.g)

    def __pow__(PowerSeries_pari self, n, m):
        """
        Exponentiation of power series.

        TESTS::

            sage: R.<t> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 3 - t^3 + O(t^5)
            sage: a = f^3; a
            27 - 27*t^3 + O(t^5)
            sage: b = f^-3; b
            1/27 + 1/27*t^3 + O(t^5)
        """
        h = self.g ** n
        if h.valuation(self._parent.variable_name()) < 0:
            return self._parent.laurent_series_ring()(h)
        return construct_from_pari(self._parent, h)

    cpdef _add_(self, right):
        """
        Addition of power series.

        TESTS::

            sage: R.<x> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = x^4 + O(x^5); f
            x^4 + O(x^5)
            sage: g = x^2 + O(x^3); g
            x^2 + O(x^3)
            sage: f+g
            x^2 + O(x^3)
        """
        return construct_from_pari(self._parent, self.g + (<PowerSeries_pari>right).g)

    cpdef _sub_(self, right):
        """
        Subtraction of power series.

        TESTS::

            sage: k.<w> = ZZ[]
            sage: R.<t> = PowerSeriesRing(k, implementation='pari')
            sage: w*t^2 -w*t +13 - (w*t^2 + w*t)
            13 - 2*w*t
        """
        return construct_from_pari(self._parent, self.g - (<PowerSeries_pari>right).g)

    cpdef _mul_(self, right):
        """
        Multiplication of power series.

        TESTS::

            sage: k.<w> = PowerSeriesRing(ZZ, implementation='pari')
            sage: (1+17*w+15*w^3+O(w^5))*(19*w^10+O(w^12))
            19*w^10 + 323*w^11 + O(w^12)
        """
        return construct_from_pari(self._parent, self.g * (<PowerSeries_pari>right).g)

    cpdef _rmul_(self, Element c):
        """
        Right multiplication by a scalar.

        TESTS::

            sage: R.<t> = PowerSeriesRing(GF(7), implementation='pari')
            sage: f = t + 3*t^4 + O(t^11)
            sage: f * GF(7)(3)
            3*t + 2*t^4 + O(t^11)
        """
        return construct_from_pari(self._parent, self.g * c)

    cpdef _lmul_(self, Element c):
        """
        Left multiplication by a scalar.

        TESTS::

            sage: R.<t> = PowerSeriesRing(GF(11), implementation='pari')
            sage: f = 1 + 3*t^4 + O(t^120)
            sage: 2 * f
            2 + 6*t^4 + O(t^120)
        """
        return construct_from_pari(self._parent, c * self.g)

    cpdef _div_(self, right):
        """
        Division of power series.

        TESTS::

            sage: R.<t> = PowerSeriesRing(GF(11), default_prec=8, implementation='pari')
            sage: f = t/(1 - t); f
            t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + O(t^8)
            sage: f.parent()
            Power Series Ring in t over Finite Field of size 11
            sage: g = (1 - t)/t; g
            t^-1 + 10
            sage: g.parent()
            Laurent Series Ring in t over Finite Field of size 11
        """
        h = self.g / (<PowerSeries_pari>right).g
        if h.valuation(self._parent.variable_name()) < 0:
            return self._parent.laurent_series_ring()(h)
        return construct_from_pari(self._parent, h)

    def list(self):
        """
        Return the list of known coefficients for ``self``.

        This is just the list of coefficients of the underlying
        polynomial; it need not have length equal to ``self.prec()``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = 1 - 5*t^3 + t^5 + O(t^7)
            sage: f.list()
            [1, 0, 0, -5, 0, 1]

            sage: # needs sage.rings.padics
            sage: S.<u> = PowerSeriesRing(pAdicRing(5), implementation='pari')
            sage: (2 + u).list()
            [2 + O(5^20), 1 + O(5^20)]
        """
        cdef pari_gen g = self.g
        cdef long vn = get_var(self._parent.variable_name())
        R = self.base_ring()
        if typ(g.g) == t_SER and varn(g.g) == vn:
            g = g.truncate()
        if typ(g.g) == t_POL and varn(g.g) == vn:
            # t_POL has 2 codewords.  Use new_ref instead of g[i] for speed.
            G = g.fixGEN()
            return [R(g.new_ref(gel(G, i))) for i in range(2, lg(G))]
        else:
            return [R(g)]

    def padded_list(self, n=None):
        """
        Return a list of coefficients of ``self`` up to (but not
        including) `q^n`.

        The list is padded with zeroes on the right so that it has
        length `n`.

        INPUT:

        - ``n`` -- nonnegative integer (optional); if `n` is not
          given, it will be taken to be the precision of ``self``,
          unless this is ``+Infinity``, in which case we just
          return ``self.list()``

        EXAMPLES::

            sage: R.<q> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 1 - 17*q + 13*q^2 + 10*q^4 + O(q^7)
            sage: f.list()
            [1, -17, 13, 0, 10]
            sage: f.padded_list(7)
            [1, -17, 13, 0, 10, 0, 0]
            sage: f.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]
            sage: f.padded_list(3)
            [1, -17, 13]
            sage: f.padded_list()
            [1, -17, 13, 0, 10, 0, 0]
            sage: g = 1 - 17*q + 13*q^2 + 10*q^4
            sage: g.list()
            [1, -17, 13, 0, 10]
            sage: g.padded_list()
            [1, -17, 13, 0, 10]
            sage: g.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]
        """
        if n is None:
            if self._prec is infinity:
                return self.list()
            else:
                n = self._prec
        if not n:
            return []

        cdef pari_gen g = self.g
        g.fixGEN()
        cdef long l, m

        R = self.base_ring()
        if typ(g.g) == t_POL and varn(g.g) == get_var(self._parent.variable_name()):
            l = lg(g.g) - 2  # t_POL has 2 codewords
            if n <= l:
                return [R(g.new_ref(gel(g.g, i + 2))) for i in range(n)]
            else:
                return ([R(g.new_ref(gel(g.g, i + 2))) for i in range(l)]
                        + [R.zero()] * (n - l))
        elif typ(g.g) == t_SER and varn(g.g) == get_var(self._parent.variable_name()):
            l = lg(g.g) - 2  # t_SER has 2 codewords
            m = valp(g.g)
            if n <= m:
                return [R.zero()] * n
            elif n <= l + m:
                return ([R.zero()] * m
                        + [R(g.new_ref(gel(g.g, i + 2))) for i in range(n - m)])
            else:
                return ([R.zero()] * m
                        + [R(g.new_ref(gel(g.g, i + 2))) for i in range(l)]
                        + [R.zero()] * (n - l - m))
        else:
            return [R(g)] + [R.zero()] * (n - 1)

    def monomial_coefficients(self, copy=None):
        """
        Return a dictionary of coefficients for ``self``.

        This is simply a dict for the underlying polynomial; it need
        not have keys corresponding to every number smaller than
        ``self.prec()``.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: f = 1 + t^10 + O(t^12)
            sage: f.monomial_coefficients()
            {0: 1, 10: 1}

        ``dict`` is an alias::

            sage: f.dict()
            {0: 1, 10: 1}
        """
        return self.polynomial().monomial_coefficients()

    dict = monomial_coefficients

    def _derivative(self, var=None):
        """
        Return the derivative of ``self`` with respect to the
        variable ``var``.

        If ``var`` is ``None``, the variable of the power series ring
        is used.

        .. SEEALSO::

            :meth:`derivative()`

        EXAMPLES::

            sage: R.<w> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 2 + 3*w^2 + w^10 + O(w^100); f
            2 + 3*w^2 + w^10 + O(w^100)
            sage: f._derivative()
            6*w + 10*w^9 + O(w^99)
            sage: f._derivative(w)
            6*w + 10*w^9 + O(w^99)

            sage: R.<t> = PolynomialRing(ZZ)
            sage: S.<x> = PowerSeriesRing(R, implementation='pari')
            sage: f = t^3*x^4 + O(x^5)
            sage: f._derivative()
            4*t^3*x^3 + O(x^4)
            sage: f._derivative(x)
            4*t^3*x^3 + O(x^4)
            sage: f._derivative(t)
            3*t^2*x^4 + O(x^5)
        """
        if var is None:
            var = self._parent.variable_name()
        return construct_from_pari(self._parent, self.g.deriv(var))

    def integral(self, var=None):
        """
        Return the formal integral of ``self``.

        By default, the integration variable is the variable of the
        power series.  Otherwise, the integration variable is the
        optional parameter ``var``.

        .. NOTE::

            The integral is always chosen so the constant term is 0.

        EXAMPLES::

            sage: k.<w> = PowerSeriesRing(QQ, implementation='pari')
            sage: (1+17*w+15*w^3+O(w^5)).integral()
            w + 17/2*w^2 + 15/4*w^4 + O(w^6)
            sage: (w^3 + 4*w^4 + O(w^7)).integral()
            1/4*w^4 + 4/5*w^5 + O(w^8)
            sage: (3*w^2).integral()
            w^3

        TESTS::

            sage: t = PowerSeriesRing(QQ, 't', implementation='pari').gen()
            sage: f = t + 5*t^2 + 21*t^3
            sage: g = f.integral() ; g
            1/2*t^2 + 5/3*t^3 + 21/4*t^4
            sage: g.parent()
            Power Series Ring in t over Rational Field

            sage: R.<a> = QQ[]
            sage: t = PowerSeriesRing(R, 't', implementation='pari').gen()
            sage: f = a*t +5*t^2
            sage: f.integral()
            1/2*a*t^2 + 5/3*t^3
            sage: f.integral(a)
            1/2*a^2*t + 5*a*t^2
        """
        if var is None:
            var = self._parent.variable_name()
        return construct_from_pari(self._parent, self.g.intformal(var))

    def reverse(self, precision=None):
        r"""
        Return the reverse of ``self``.

        The reverse of a power series `f` is the power series `g` such
        that `g(f(x)) = x`.  This exists if and only if the valuation
        of ``self`` is exactly 1 and the coefficient of `x` is a unit.

        If the optional argument ``precision`` is given, the reverse
        is returned with this precision.  If ``f`` has infinite
        precision and the argument ``precision`` is not given, then
        the reverse is returned with the default precision of
        ``f.parent()``.

        EXAMPLES::

            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 2*x + 3*x^2 - x^4 + O(x^5)
            sage: g = f.reverse()
            sage: g
            1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)
            sage: f(g)
            x + O(x^5)
            sage: g(f)
            x + O(x^5)

            sage: A.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: a = t - t^2 - 2*t^4 + t^5 + O(t^6)
            sage: b = a.reverse(); b
            t + t^2 + 2*t^3 + 7*t^4 + 25*t^5 + O(t^6)
            sage: a(b)
            t + O(t^6)
            sage: b(a)
            t + O(t^6)

            sage: B.<b,c> = PolynomialRing(ZZ)
            sage: A.<t> = PowerSeriesRing(B, implementation='pari')
            sage: f = t + b*t^2 + c*t^3 + O(t^4)
            sage: g = f.reverse(); g
            t - b*t^2 + (2*b^2 - c)*t^3 + O(t^4)
            sage: f(g)
            t + O(t^4)
            sage: g(f)
            t + O(t^4)

            sage: A.<t> = PowerSeriesRing(ZZ, implementation='pari')
            sage: B.<x> = PowerSeriesRing(A, implementation='pari')
            sage: f = (1 - 3*t + 4*t^3 + O(t^4))*x + (2 + t + t^2 + O(t^3))*x^2 + O(x^3)
            sage: g = f.reverse(); g
            (1 + 3*t + 9*t^2 + 23*t^3 + O(t^4))*x + (-2 - 19*t - 118*t^2 + O(t^3))*x^2 + O(x^3)

        The optional argument ``precision`` sets the precision of the output::

            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
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

            sage: R.<x> = PowerSeriesRing(QQ, default_prec=20, implementation='pari')
            sage: (x - x^2).reverse()  # get some Catalan numbers
            x + x^2 + 2*x^3 + 5*x^4 + 14*x^5 + 42*x^6 + 132*x^7 + 429*x^8
             + 1430*x^9 + 4862*x^10 + 16796*x^11 + 58786*x^12 + 208012*x^13
             + 742900*x^14 + 2674440*x^15 + 9694845*x^16 + 35357670*x^17
             + 129644790*x^18 + 477638700*x^19 + O(x^20)
            sage: (x - x^2).reverse(precision=3)
            x + x^2 + O(x^3)

        TESTS::

            sage: R.<x> = PowerSeriesRing(QQ, implementation='pari')
            sage: f = 1 + 2*x + 3*x^2 - x^4 + O(x^5)
            sage: f.reverse()
            Traceback (most recent call last):
            ...
            PariError: domain error in serreverse: valuation != 1
        """
        cdef PowerSeries_pari f
        if self._prec is infinity:
            if precision is None:
                precision = self._parent.default_prec()
            f = self.add_bigoh(precision)
        else:
            if precision is None:
                precision = self._prec
            f = self
        return PowerSeries_pari(self._parent, f.g.serreverse(), precision)
