# sage_setup: distribution = sagemath-categories
r"""
Puiseux Series Ring Element

A Puiseux series is a series of the form

.. MATH::

    p(x) = \sum_{n=N}^{\infty} a_n (x-a)^{n/e},

where the integer :math:`e` is called the *ramification index* of the series
and the number :math:`a` is the *center*. A Puiseux series is essentially a
Laurent series but with fractional exponents.

EXAMPLES:

We begin by constructing the ring of Puiseux series in `x` with coefficients
in the rationals::

    sage: R.<x> = PuiseuxSeriesRing(QQ)

This command also defines ``x`` as the generator of this ring.

When constructing a Puiseux series, the ramification index is automatically
determined from the greatest common divisor of the exponents::

    sage: p = x^(1/2); p
    x^(1/2)
    sage: p.ramification_index()
    2
    sage: q = x^(1/2) + x**(1/3); q
    x^(1/3) + x^(1/2)
    sage: q.ramification_index()
    6

Other arithmetic can be performed with Puiseux Series::

    sage: p + q
    x^(1/3) + 2*x^(1/2)
    sage: p - q
    -x^(1/3)
    sage: p * q
    x^(5/6) + x
    sage: (p / q).add_bigoh(4/3)
    x^(1/6) - x^(1/3) + x^(1/2) - x^(2/3) + x^(5/6) - x + x^(7/6) + O(x^(4/3))

Mind the base ring. However, the base ring can be changed::

    sage: I*q                                                                           # needs sage.rings.number_field
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for *:
    'Number Field in I with defining polynomial x^2 + 1 with I = 1*I' and
    'Puiseux Series Ring in x over Rational Field'
    sage: qz = q.change_ring(ZZ); qz
    x^(1/3) + x^(1/2)
    sage: qz.parent()
    Puiseux Series Ring in x over Integer Ring

Other properties of the Puiseux series can be easily obtained::

    sage: r = (3*x^(-1/5) + 7*x^(2/5) + (1/2)*x).add_bigoh(6/5); r
    3*x^(-1/5) + 7*x^(2/5) + 1/2*x + O(x^(6/5))
    sage: r.valuation()
    -1/5
    sage: r.prec()
    6/5
    sage: r.precision_absolute()
    6/5
    sage: r.precision_relative()
    7/5
    sage: r.exponents()
    [-1/5, 2/5, 1]
    sage: r.coefficients()
    [3, 7, 1/2]

Finally, Puiseux series are compatible with other objects in Sage.
For example, you can perform arithmetic with Laurent series::

    sage: L.<x> = LaurentSeriesRing(ZZ)
    sage: l = 3*x^(-2) + x^(-1) + 2 + x**3
    sage: r + l
    3*x^-2 + x^-1 + 3*x^(-1/5) + 2 + 7*x^(2/5) + 1/2*x + O(x^(6/5))

AUTHORS:

- Chris Swierczewski 2016: initial version on https://github.com/abelfunctions/abelfunctions/tree/master/abelfunctions
- Frédéric Chapoton 2016: integration of code
- Travis Scrimshaw, Sebastian Oehms 2019-2020: basic improvements and completions

REFERENCES:

- :wikipedia:`Puiseux_series`
"""


# ****************************************************************************
#       Copyright (c) 2016 Chris Swierczewski
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import infinity
from sage.rings.laurent_series_ring_element cimport LaurentSeries
from sage.structure.element cimport (Element, AlgebraElement)
from sage.structure.richcmp cimport richcmp


cdef class PuiseuxSeries(AlgebraElement):
    r"""
    A Puiseux series.

    .. MATH::

        \sum_{n=-N}^\infty a_n x^{n/e}

    It is stored as a Laurent series:

    .. MATH::

        \sum_{n=-N}^\infty a_n t^n

    where `t = x^{1/e}`.

    INPUT:

    - ``parent`` -- the parent ring

    - ``f`` -- one of the following types of inputs:

      * instance of :class:`PuiseuxSeries`
      * instance that can be coerced into the Laurent series ring of the parent

    - ``e`` -- integer (default: 1); the ramification index

    EXAMPLES::

        sage: R.<x> = PuiseuxSeriesRing(QQ)
        sage: p = x^(1/2) + x**3; p
        x^(1/2) + x^3
        sage: q = x**(1/2) - x**(-1/2)
        sage: r = q.add_bigoh(7/2); r
        -x^(-1/2) + x^(1/2) + O(x^(7/2))
        sage: r**2
        x^-1 - 2 + x + O(x^3)
    """

    def __init__(self, parent, f, e=1):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + x**3
            sage: TestSuite(p).run()
        """
        AlgebraElement.__init__(self, parent)
        L = parent._laurent_series_ring

        if isinstance(f, PuiseuxSeries):
            if (<PuiseuxSeries>f)._l._parent is L:
                l = (<PuiseuxSeries>f)._l
                e = (<PuiseuxSeries>f)._e
            else:
                l = L((<PuiseuxSeries>f)._l)
                e = L((<PuiseuxSeries>f)._e)
        else:
            l = L(f)

        # --------------------------------------------------------
        # choose a representative for this Puiseux series having
        # minimal ramification index. This is necessary because
        # some methods need it as minimal as possible (for example
        # :meth:`laurent_series' or :meth:`power_series`)
        # --------------------------------------------------------
        exp_list = l.exponents()
        prec     = l.prec()
        if prec is infinity:
            d = gcd(exp_list +[e])
        else:
            d = gcd(exp_list + [e] + [prec])
        if d > 1:
            # ramification index can be reduced dividing by d
            e = e / d
            cf_ori = l.list()
            if cf_ori:
                cf = [cf_ori[d*i] for i in range((len(cf_ori)-1) / d + 1)]
            else:
                cf = cf_ori
            val = l.valuation() / d
            l = l.parent()(cf, n=val)
            if prec != infinity:
                l = l.add_bigoh(prec / d)

        self._l = l
        self._e = long(abs(e))

    def __reduce__(self):
        """
        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(1/2) + x**3-x**(-1/4)
            sage: loads(dumps(p)) == p    # indirect doctest
            True
        """
        return (self._parent, (self._l, self._e))

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(1/3) + x**3
            sage: t = p._im_gens_(QQbar, [2])
            sage: t in QQbar
            True
            sage: f = R.hom([QQbar(2)], check=False)
            sage: t == f(p)
            True
        """
        return self(codomain(im_gens[0]))

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + x**3-x**(-1/4); p
            -x^(-1/4) + x^(1/2) + x^3
            sage: R.zero()
            0

            sage: S.<t> = PuiseuxSeriesRing(Zp(5))                                      # needs sage.rings.padics
            sage: t**(1/2) + 5 * t^(1/3)                                                # needs sage.rings.padics
            (5 + O(5^21))*t^(1/3) + (1 + O(5^20))*t^(1/2)
        """
        laurent = self.laurent_part()
        s = repr(laurent)
        if self.ramification_index() == 1:
            return s

        X = self._parent.variable_name()

        # find a temporary variable name (to avoid multiple transformations)
        Xtemp = '?'
        while Xtemp in s: Xtemp +='?' # if somebody uses '?' in variable_name

        # renaming and generalizing linear term
        s = s.replace('%s' %X, '%s^1' %Xtemp)
        s = s.replace('^1^', '^' )

        # prepare exponent list
        if laurent.prec() is infinity:
            exponents = [ZZ(exp) for exp in set(laurent.exponents())]
        else:
            exponents = [ZZ(exp) for exp in set(laurent.exponents() + [laurent.prec()])]

        # sort exponents such that the largest will be replaced first
        exp_pos = [exp for exp in exponents if exp >= 0]
        exp_pos.sort(reverse=True)
        exp_neg = [exp for exp in exponents if exp < 0]
        exp_neg.sort()
        exponents = exp_neg + exp_pos

        # replacing exponents
        e = ZZ(self.ramification_index())
        for exp_l in exponents:
            exp = exp_l/e
            repl_str = '%s^%s' %(Xtemp, exp_l)
            if exp.is_one():
                s = s.replace(repl_str, '%s' %X)
            elif e.divides(exp_l):
                s = s.replace(repl_str, '%s^%s' %(X, exp))
            else:
                s = s.replace(repl_str, '%s^(%s)' %(X, exp))
        return s

    def __call__(self, x):
        r"""
        Evaluate this Puiseux series.

        INPUT:

        - ``x`` -- element of a ring

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + x**3-x**(-1/4)
            sage: p(16)
            8199/2
            sage: p(pi.n())                                                             # needs sage.symbolic
            32.0276049867404
        """
        # use x.nth_root since x**(1/self._e) returns oo when x = 0
        if isinstance(x, int):
            x = ZZ(x)
        elif isinstance(x, float):
            from sage.rings.complex_mpfr import ComplexField
            x = ComplexField()(x)
        t = x.nth_root(self._e)
        p = self._l.__u.polynomial()
        n = self._l.__n
        return p(t)*t**n

    def _common_ramification_index(self, PuiseuxSeries right):
        r"""
        Return a ramification index common to ``self`` and ``right``.

        In order to perform arithmetic on Puiseux series it is useful to find a
        common ramification index between two operands. That is, given Puiseux
        series :math:`p` and :math:`q` of ramification indices :math:`e` and
        :math:`f` we write both as series :math:`\tilde{f}` and
        :math:`\tilde{g}` in :math:`(x-a)^(1/g)` such that,

        .. MATH::

            f = \tilde{f}((x-a)^M), g = \tilde{g}((x-a)^N).

        INPUT:

        - ``right`` -- a Puiseux series

        OUTPUT:

        - ``g`` -- integer; a ramification index common to ``self`` and
          ``right``
        - ``M``, ``N`` -- integers; scaling factors on ``self`` and ``right``,
          respectively

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/3) + x**2-x**(-1/7)
            sage: q = x^(-1/3) + x**2-x**(1/5)
            sage: p._common_ramification_index(q)
            (105, 5, 7)
            sage: q._common_ramification_index(p)
            (105, 7, 5)
        """
        m = self._e
        n = right._e
        g = lcm(m, n)
        m = g / m
        n = g / n
        return g, m, n

    cpdef _add_(self, right_m):
        """
        Return the sum.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p + q                                        # indirect doctest
            2*x^(1/3) + 3/4*x^(2/5) + x^(1/2) + 3/4*x^(2/3)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_m
        cdef LaurentSeries l, l1, l2
        cdef size_t g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self._l.V(m)
        l2 = right._l.V(n)
        l = l1 + l2
        return type(self)(self._parent, l, g)

    cpdef _sub_(self, right_m):
        """
        Return the difference.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p - q                                        # indirect doctest
            -2*x^(1/3) - 3/4*x^(2/5) + x^(1/2) + 3/4*x^(2/3)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_m
        cdef LaurentSeries l, l1, l2
        cdef size_t g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self._l.V(m)
        l2 = right._l.V(n)
        l = l1 - l2
        return type(self)(self._parent, l, g)

    cpdef _mul_(self, right_r):
        """
        Return the product.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p * q                                        # indirect doctest
            2*x^(5/6) + 3/4*x^(9/10) + 3/2*x + 9/16*x^(16/15)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_r
        cdef LaurentSeries l, l1, l2
        cdef size_t g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self._l.V(m)
        l2 = right._l.V(n)
        l = l1 * l2
        return type(self)(self._parent, l, g)

    cpdef _rmul_(self, Element c):
        """
        Return the right scalar multiplication.

        EXAMPLES::

            sage: P.<y> = PuiseuxSeriesRing(ZZ)
            sage: t = y^(-1/3) + O(y^(0))
            sage: 5*t                                          # indirect doctest
            5*y^(-1/3) + O(1)
        """
        return type(self)(self._parent, self._l._rmul_(c), self._e)

    cpdef _lmul_(self, Element c):
        """
        Return the left scalar multiplication.

        EXAMPLES::

            sage: P.<y> = PuiseuxSeriesRing(Zp(3))                                      # needs sage.rings.padics
            sage: t = y^(2/5) + O(y)                                                    # needs sage.rings.padics
            sage: 5*t                                          # indirect doctest       # needs sage.rings.padics
            (2 + 3 + O(3^20))*y^(2/5) + O(y)
        """
        return type(self)(self._parent, self._l._lmul_(c), self._e)

    cpdef _div_(self, right_r):
        """
        Return the quotient.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p / q
            1/2*x^(1/6) - 3/16*x^(7/30) + 9/128*x^(3/10) + 3/8*x^(1/3)
             - 27/1024*x^(11/30) - 9/64*x^(2/5) + 81/8192*x^(13/30)
             + 27/512*x^(7/15) - 243/65536*x^(1/2) - 81/4096*x^(8/15)
             + 729/524288*x^(17/30) + 243/32768*x^(3/5) - 2187/4194304*x^(19/30)
             - 729/262144*x^(2/3) + 6561/33554432*x^(7/10)
             + 2187/2097152*x^(11/15) - 19683/268435456*x^(23/30)
             - 6561/16777216*x^(4/5) + O(x^(5/6))
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_r
        cdef LaurentSeries l, l1, l2
        cdef size_t g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self._l.V(m)
        l2 = right._l.V(n)
        l = l1 / l2
        return type(self)(self._parent, l, g)

    def __pow__(_self, r, dummy):
        """
        Return the power.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p ** 3
            x^(3/2) + 9/4*x^(5/3) + 27/16*x^(11/6) + 27/64*x^2
        """
        cdef PuiseuxSeries self = _self
        cdef LaurentSeries l
        cdef size_t e

        r = QQ(r)
        numer = r.numerator()
        denom = r.denominator()

        # if the exponent is integral then do normal exponentiation
        if denom == 1:
            l = self._l ** int(numer)
            e = self._e
        # otherwise, we only exponentiate by a rational number if there is a
        # single term in the Puiseux series
        #
        # (I suppose we could use Taylor series expansions in the general case)
        else:
            if not self.is_monomial():
                raise ValueError('can only exponentiate single term by rational')
            l = self._l.V(numer)
            e = self._e * int(denom)
        return type(self)(self._parent, l, e)

    cpdef _richcmp_(self, right_r, int op):
        r"""
        Comparison of ``self`` and ``right``.

        We say two approximate Puiseux series are equal, if they agree for
        all coefficients up to the *minimum* of the precisions of each.

        Comparison is done in dictionary order going from lowest degree
        to highest degree coefficients with respect to the corresponding
        Laurent series. That means that comparison is performed for
        corresponding `LaurentSeries` instances obtained for the common
        ramification index.

        See :meth:`power_series_ring_element._richcmp_` and
        :meth:`_laurent_series_ring_element._richcmp_` for more
        information.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: (p < q, p >= q, p == 0, q != 0)
            (True, False, False, True)
            sage: p2 = x^(1/2) + 3/5 * x^(2/3)
            sage: (p2 < p, p2 >= p)
            (True, False)
            sage: p3 = p2.add_bigoh(2/3); p3
            x^(1/2) + O(x^(2/3))
            sage: (p3 == p2, p3 == p, p3 < q, p3 > q)
            (True, True, True, False)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_r
        if self._e == right._e:
            return richcmp(self._l, right._l, op)

        # If both have different ramification indices they must be different as
        # Puiseux series (by the normalization performed in the python constructor).
        # We use the ramification index to order them.
        return richcmp(self._e, right._e, op)

    def __lshift__(self, r):
        """
        Apply :meth:`shift` using the operator `<<`.

        EXAMPLES::

            sage: P.<y> = LaurentPolynomialRing(ZZ)
            sage: R.<x> = PuiseuxSeriesRing(P)
            sage: p = y*x**(-1/3) + 2*y^(-2)*x**(1/2)
            sage: p << 1/3                             # indirect doctest
            y + (2*y^-2)*x^(5/6)
        """
        return self.shift(r)

    def __rshift__(self, r):
        """
        Apply :meth:`shift` with negative argument using the operator `>>`.

        EXAMPLES::

            sage: P.<y> = LaurentPolynomialRing(ZZ)
            sage: R.<x> = PuiseuxSeriesRing(P)
            sage: p = y*x**(-1/3) + 2*y^(-2)*x**(1/2)
            sage: p >> 1/3                             # indirect doctest
            y*x^(-2/3) + (2*y^-2)*x^(1/6)
        """
        return self.shift(-r)

    def __bool__(self):
        """
        Return whether ``self`` is not zero.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.is_zero()                     # indirect doctest
            False
            sage: R.zero() != 0
            False
        """
        return bool(self._l)

    def __hash__(self):
        """
        Return a hash of ``self``.

        EXAMPLES::

            sage: from operator import xor
            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: hash(p) == xor(hash(p.laurent_part()), 2)  # indirect doctest
            True
        """
        return hash(self._l) ^ self._e

    def __getitem__(self, r):
        r"""
        Return the coefficient with exponent ``r``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p[-7/2]
            1
            sage: p[0]
            3
            sage: p[1/2]
            5
            sage: p[3]
            -7
            sage: p[100]
            0
        """
        if isinstance(r, slice):
            start, stop, step = r.start, r.stop, r.step
            n = slice(start * self._e, stop * self._e, step * self._e)
            return PuiseuxSeries(self._parent, self._l[start:stop:step], self._e)
        else:
            n = int(r * self._e)
            return self._l[n]

    def __iter__(self):
        """
        Return an iterator over the coefficients.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: list(p)
            [1, 0, 0, 0, 0, 0, 0, 3, 5, 0, 0, 0, 0, -7]
        """
        return iter(self._l)

    def __copy__(self):
        """
        Since this is immutable, return ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(3/4) + 2*x^(4/5) + 3* x^(5/6)
            sage: p2 = copy(p); p2
            x^(3/4) + 2*x^(4/5) + 3*x^(5/6)
            sage: p == p2
            True
            sage: p is p2
            True
        """
        return self

    def laurent_part(self):
        """
        Return the underlying Laurent series.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.laurent_part()
            x^3 + 3/4*x^4
        """
        return self._l

    def ramification_index(self):
        """
        Return the ramification index.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.ramification_index()
            6
        """
        return self._e

    def valuation(self):
        r"""
        Return the valuation of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.valuation()
            -7/2

        TESTS::

            sage: R.zero().valuation()
            +Infinity
        """
        return self._l.valuation() / QQ(self._e)

    def add_bigoh(self, prec):
        r"""
        Return the truncated series at chosen precision ``prec``.

        INPUT:

        - ``prec`` -- the precision of the series as a rational number

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.add_bigoh(2)
            x^(-7/2) + 3 + 5*x^(1/2) + O(x^2)
            sage: p.add_bigoh(0)
            x^(-7/2) + O(1)
            sage: p.add_bigoh(-1)
            x^(-7/2) + O(x^-1)

        .. NOTE::

            The precision passed to the method is adapted to the common
            ramification index::

                sage: R.<x> = PuiseuxSeriesRing(ZZ)
                sage: p = x**(-1/3) + 2*x**(1/5)
                sage: p.add_bigoh(1/2)
                x^(-1/3) + 2*x^(1/5) + O(x^(7/15))
        """
        if prec is infinity or prec >= self.prec():
            return self

        l_prec = int(prec * self._e)
        l = self._l.add_bigoh(l_prec)
        return PuiseuxSeries(self._parent, l, self._e)

    def change_ring(self, R):
        r"""
        Return ``self`` over a the new ring ``R``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: q = p.change_ring(QQ); q
            x^(-7/2) + 3 + 5*x^(1/2) - 7*x^3
            sage: q.parent()
            Puiseux Series Ring in x over Rational Field
        """
        return self._parent.change_ring(R)(self)

    def is_unit(self):
        r"""
        Return whether ``self`` is a unit.

        A Puiseux series is a unit if and only if its leading coefficient is.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.is_unit()
            True
            sage: q = 4 * x^(-7/2) + 3 * x**4
            sage: q.is_unit()
            False
        """
        return self._l.is_unit()

    def is_zero(self):
        """
        Return whether ``self`` is zero.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.is_zero()
            False
            sage: R.zero().is_zero()
            True
        """
        return self._l.is_zero()

    def is_monomial(self):
        r"""
        Return whether ``self`` is a monomial.

        This is ``True`` if and only if ``self`` is `x^p` for
        some rational `p`.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.is_monomial()
            False
            sage: q = x**(11/13)
            sage: q.is_monomial()
            True
            sage: q = 4*x**(11/13)
            sage: q.is_monomial()
            False
        """
        return self._l.is_monomial()

    def list(self):
        r"""
        Return the list of coefficients indexed by the exponents of the
        the corresponding Laurent series.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(3/4) + 2*x^(4/5) + 3* x^(5/6)
            sage: p.list()
            [1, 0, 0, 2, 0, 3]
        """
        return self._l.list()

    def coefficients(self):
        r"""
        Return the list of coefficients.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(3/4) + 2*x^(4/5) + 3* x^(5/6)
            sage: p.coefficients()
            [1, 2, 3]
        """
        return self._l.coefficients()

    def exponents(self):
        r"""
        Return the list of exponents.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(3/4) + 2*x^(4/5) + 3* x^(5/6)
            sage: p.exponents()
            [3/4, 4/5, 5/6]
        """
        return [QQ(n) /  self._e for n in self._l.exponents()]

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: R.<t> = PuiseuxSeriesRing(QQ)
            sage: f = t^2 + t^3 + O(t^10)
            sage: f[2] = 5
            Traceback (most recent call last):
            ...
            IndexError: Puiseux series are immutable
        """
        raise IndexError('Puiseux series are immutable')

    def degree(self):
        r"""
        Return the degree of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: P.<y> = PolynomialRing(GF(5))
            sage: R.<x> = PuiseuxSeriesRing(P)
            sage: p = 3*y*x**(-2/3) + 2*y**2*x**(1/5); p
            3*y*x^(-2/3) + 2*y^2*x^(1/5)
            sage: p.degree()
            1/5
        """
        return self._l.degree() / self._e

    def shift(self, r):
        r"""
        Return this Puiseux series multiplied by `x^r`.

        EXAMPLES::

            sage: P.<y> = LaurentPolynomialRing(ZZ)
            sage: R.<x> = PuiseuxSeriesRing(P)
            sage: p = y*x**(-1/3) + 2*y^(-2)*x**(1/2); p
            y*x^(-1/3) + (2*y^-2)*x^(1/2)
            sage: p.shift(3)
            y*x^(8/3) + (2*y^-2)*x^(7/2)
        """
        cdef LaurentSeries l = self._l.shift(r * self._e)
        return PuiseuxSeries(self._parent, l, self._e)

    def truncate(self, r):
        r"""
        Return the Puiseux series of degree `< r`.

        This is equivalent to ``self`` modulo `x^r`.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = (x**(-1/3) + 2*x**3)**2; p
            x^(-2/3) + 4*x^(8/3) + 4*x^6
            sage: q = p.truncate(5); q
            x^(-2/3) + 4*x^(8/3)
            sage: q == p.add_bigoh(5)
            True
        """
        l = self._l.truncate(r * self._e)
        return PuiseuxSeries(self._parent, l, self._e)

    def prec(self):
        r"""
        Return the precision of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = (x**(-1/3) + 2*x**3)**2; p
            x^(-2/3) + 4*x^(8/3) + 4*x^6
            sage: q = p.add_bigoh(5); q
            x^(-2/3) + 4*x^(8/3) + O(x^5)
            sage: q.prec()
            5
        """
        if self._l.prec() is infinity:
            return infinity
        return self._l.prec() / self._e

    precision_absolute = prec

    def precision_relative(self):
        r"""
        Return the relative precision of the series.

        The relative precision of the Puiseux series is the difference
        between its absolute precision and its valuation.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = PuiseuxSeriesRing(GF(3))
            sage: p = (x**(-1/3) + 2*x**3)**2; p
            x^(-2/3) + x^(8/3) + x^6
            sage: q = p.add_bigoh(7); q
            x^(-2/3) + x^(8/3) + x^6 + O(x^7)
            sage: q.precision_relative()
            23/3
        """
        if self.is_zero():
            return 0
        return self.prec() - self.valuation()

    def common_prec(self, p):
        r"""
        Return the minimum precision of `p` and ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = (x**(-1/3) + 2*x**3)**2
            sage: q5 = p.add_bigoh(5); q5
            x^(-2/3) + 4*x^(8/3) + O(x^5)
            sage: q7 = p.add_bigoh(7); q7
            x^(-2/3) + 4*x^(8/3) + 4*x^6 + O(x^7)
            sage: q5.common_prec(q7)
            5
            sage: q7.common_prec(q5)
            5
        """
        if self.prec() is infinity:
            return p.prec()
        elif p.prec() is infinity:
            return self.prec()
        return min(self.prec(), p.prec())

    def variable(self):
        r"""
        Return the variable of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.variable()
            'x'
        """
        return self._parent.variable_name()

    def laurent_series(self):
        r"""
        If ``self`` is a Laurent series, return it as a Laurent series.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x**(1/2) - x**(-1/2)
            sage: p.laurent_series()
            Traceback (most recent call last):
            ...
            ArithmeticError: self is not a Laurent series
            sage: q = p**2
            sage: q.laurent_series()
            x^-1 - 2 + x
        """
        if self._e != 1:
            raise ArithmeticError('self is not a Laurent series')
        return self._l

    def power_series(self):
        r"""
        If ``self`` is a power series, return it as a power series.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PuiseuxSeriesRing(QQbar)
            sage: p = x**(3/2) - QQbar(I)*x**(1/2)
            sage: p.power_series()
            Traceback (most recent call last):
            ...
            ArithmeticError: self is not a power series
            sage: q = p**2
            sage: q.power_series()
            -x - 2*I*x^2 + x^3
        """
        try:
            l = self.laurent_series()
            return l.power_series()
        except ArithmeticError:
            raise ArithmeticError('self is not a power series')

    def inverse(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: 1/p
            x^(7/2) - 3*x^7 - 5*x^(15/2) + 7*x^10 + 9*x^(21/2) + 30*x^11 +
            25*x^(23/2) + O(x^(27/2))
        """
        return self.__invert__()
