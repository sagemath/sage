r"""
Lazy Laurent Series

A lazy Laurent series is a Laurent series whose coefficients are computed as
demanded or needed. Unlike the usual Laurent series in Sage, lazy Laurent
series do not have precisions because a lazy Laurent series knows (can be
computed, lazily) all its coefficients.

EXAMPLES:

Generating functions are Laurent series over the integer ring::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)

This defines the generating function of Fibonacci sequence::

    sage: def coeff(s, i):
    ....:     if i in [0, 1]:
    ....:         return 1
    ....:     else:
    ....:         return s.coefficient(i - 1) + s.coefficient(i - 2)
    sage: f = L.series(coeff, valuation=0); f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...

The 100th element of Fibonacci sequence can be obtained from the generating
function::

    sage: f.coefficient(100)
    573147844013817084101

Coefficients are computed and cached only when necessary::

    sage: f._cache[100]
    573147844013817084101
    sage: f._cache[101]
    Traceback (most recent call last):
    ...
    KeyError: 101

You can do arithmetic with lazy power series::

    sage: f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: f^-1
    1 - z - z^2 + ...
    sage: f + f^-1
    2 + z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: g = (f + f^-1)*(f - f^-1); g
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + ...

You may need to change the base ring::

    sage: h = g.change_ring(QQ)
    sage: h.parent()
    Lazy Laurent Series Ring in z over Rational Field
    sage: h
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + ...
    sage: h^-1
    1/4*z^-1 - 3/8 + 1/16*z - 17/32*z^2 + 5/64*z^3 - 29/128*z^4 + 165/256*z^5 + ...
    sage: _.valuation()
    -1

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .infinity import infinity

from sage.structure.element import ModuleElement
from sage.structure.richcmp import op_EQ, op_NE

from sage.arith.power import generic_power

from .lazy_laurent_series_operator import (
    LazyLaurentSeriesOperator_mul,
    LazyLaurentSeriesOperator_div,
    LazyLaurentSeriesOperator_add,
    LazyLaurentSeriesOperator_sub,
    LazyLaurentSeriesOperator_neg,
    LazyLaurentSeriesOperator_inv,
    LazyLaurentSeriesOperator_scale,
    LazyLaurentSeriesOperator_apply,
    LazyLaurentSeriesOperator_change_ring,
    LazyLaurentSeriesOperator_truncate
)


class LazyLaurentSeries(ModuleElement):
    r"""
    Return a lazy Laurent series.

    INPUT:

    - ``coefficient`` -- Python function that computes coefficients

    - ``valuation`` -- integer; approximate valuation of the series

    - ``constant`` -- either ``None`` or pair of an element of the base ring and an integer

    Let the coefficient of index `i` mean the coefficient of the term of the
    series with exponent `i`.

    Python function ``coefficient`` returns the value of the coefficient of
    index `i` from input `s` and `i` where `s` is the series itself.

    Let ``valuation`` be `n`. All coefficients of index below `n` are zero.  If
    ``constant`` is ``None``, then the ``coefficient`` function is responsible
    to compute the values of all coefficients of index `\ge n`. If ``constant``
    is a pair `(c,m)`, then the ``coefficient`` function is responsible to
    compute the values of all coefficients of index `\ge n` and `< m` and all
    the coefficients of index `\ge m` is the constant `c`.

    EXAMPLES::

        sage: L = LazyLaurentSeriesRing(ZZ, 'z')
        sage: L.series(lambda s, i: i, valuation=-3, constant=(-1,3))
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + ...

    ::

        sage: def coeff(s, i):
        ....:     if i in [0, 1]:
        ....:         return 1
        ....:     else:
        ....:         return s.coefficient(i - 1) + s.coefficient(i - 2)
        sage: f = L.series(coeff, valuation=0); f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: f.coefficient(100)
        573147844013817084101

    Lazy Laurent series is picklable::

        sage: z = L.gen()
        sage: f = 1/(1 - z - z^2)
        sage: f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: g = loads(dumps(f))
        sage: g
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: g == f
        True
    """
    def __init__(self, parent, coefficient=None, valuation=0, constant=None):
        """
        Initialize.

        TESTS::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: z = L.gen()
            sage: TestSuite(z).run()
        """
        ModuleElement.__init__(self, parent)

        self._coefficient_function = coefficient
        self._approximate_valuation = valuation
        self._constant = constant

        self._cache = dict() # cache of known coefficients

    def _richcmp_(self, other, op):
        """
        Compare ``self` with ``other`` with respect to the comparison operator ``op``.

        Equality is verified if the corresponding coefficients of both series
        can be checked for equality without computing coefficients
        indefinitely.  Otherwise an exception is raised to declare that
        equality is not decidable.

        Inequality is not defined for lazy Laurent series.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: z + z^2 == z^2 + z
            True
            sage: z + z^2 != z^2 + z
            False
            sage: z + z^2 > z^2 + z
            False
            sage: z + z^2 < z^2 + z
            False
        """
        if op is op_EQ:
            if self._constant is None:
                if other._constant is None:
                    n = min(self._approximate_valuation, other._approximate_valuation)
                    m = max(self._approximate_valuation, other._approximate_valuation)
                    for i in range(n, m):
                        if self.coefficient(i) != other.coefficient(i):
                            return False
                    if self._coefficient_function == other._coefficient_function:
                        return True
                    raise ValueError("undecidable as lazy Laurent series")
                else:
                    raise ValueError("undecidable as lazy Laurent series")
            elif other._constant is None:
                raise ValueError("undecidable as lazy Laurent series")

            sc, sm = self._constant
            oc, om = other._constant

            if sc != oc:
                return False

            n = self._approximate_valuation
            m = max(sm, om)

            for i in range(n, m):
                if self.coefficient(i) != other.coefficient(i):
                    return False

            return True

        if op is op_NE:
            return not (self == other)

        return False

    def __hash__(self):
        """
        Return the hash of ``self``

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5)
            sage: g = (1 + f)/(1 - f)^2
            sage: {g: 1}
            {z^5 - 2*z^6 + z^7 + 5*z^9 - 11*z^10 + z^11 + ...: 1}
        """
        return hash((type(self), self._coefficient_function,
                     self._approximate_valuation, self._constant))

    def __bool__(self):
        """
        Test whether ``self`` is not zero.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(GF(2))
            sage: (z-z).is_zero()
            True
            sage: f = 1/(1 - z)
            sage: f.is_zero()
            False
        """
        if self._constant is None:
            for a in self._cache:
                if a:
                    return True
            if self.coefficient(self._approximate_valuation):
                return True
            raise ValueError("undecidable as lazy Laurent series")

        sc, sm = self._constant

        if sc:
            return True

        for i in range(self._approximate_valuation, sm):
            if self.coefficient(i):
                return True

        return False

    # for Python 2 compatibility
    __nonzero__ = __bool__

    def _repr_(self):
        """
        Return the string representation of this Laurent series.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: -1/(1 + 2*z)
            -1 + 2*z - 4*z^2 + 8*z^3 - 16*z^4 + 32*z^5 - 64*z^6 + ...
        """
        atomic_repr = self.base_ring()._repr_option('element_is_atomic')
        X = self.parent().variable_name()

        n = self.valuation()

        if self._constant is None:
            m = n + 7 # long enough
        elif self._constant[0] != 0:
            m = self._constant[1] + 3
        else:
            m = self._constant[1]

        s = ' '
        first = True
        while n < m:
            x = repr(self.coefficient(n))
            if x != '0':
                if not first:
                    s += ' + '
                if not atomic_repr and n > 0 and (x[1:].find('+') != -1 or x[1:].find('-') != -1):
                    x = '({})'.format(x)
                if n > 1 or n < 0:
                    var = '*{}^{}'.format(X,n)
                elif n == 1:
                    var = '*{}'.format(X)
                else:  # n == 0
                    var = ''
                s += '{}{}'.format(x,var)
                first = False
            n += 1

        s = s.replace(" + -", " - ").replace(" 1*"," ").replace(" -1*", " -")[1:]

        if not s:  # zero series
            s = '0'

        if self._constant is None or self._constant[1] > m or self._constant[0] != 0:
            s += ' + {}'.format('...')

        return s

    def __getitem__(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = z/(1 - 2*z^3)
            sage: [f[n] for n in range(20)]
            [0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 16, 0, 0, 32, 0, 0, 64]
        """
        return self.coefficient(n)

    def coefficient(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: def g(s, i):
            ....:     if i == 0:
            ....:         return 1
            ....:     else:
            ....:         return sum(s.coefficient(j)*s.coefficient(i - 1 -j) for j in [0..i-1])
            sage: e = L.series(g, valuation=0)
            sage: e.coefficient(10)
            16796
            sage: e
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...
        """
        R = self.base_ring()

        if self._approximate_valuation == infinity:
            return R.zero()
        elif n < self._approximate_valuation:
            return R.zero()
        elif self._constant is not None and n >= self._constant[1]:
            return self._constant[0]

        try:
            c = self._cache[n]
        except KeyError:
            c = R(self._coefficient_function(self, n))
            self._cache[n] = c

        return c

    def valuation(self):
        """
        Return the valuation of the series.

        This method determines the valuation of the series by looking for a
        nonzero coefficient. Hence if the series happens to be zero, then it
        may run forever.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: s = 1/(1 - z) - 1/(1 - 2*z)
            sage: s.valuation()
            1
            sage: t = z - z
            sage: t.valuation()
            +Infinity
        """
        if self._constant is None:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n in cache:
                    if cache[n]:
                        self._approximate_valuation = n
                        return n
                    n += 1
                else:
                    if self.coefficient(n) != 0:
                        self._approximate_valuation = n
                        return n
                    n += 1
        else:
            n = self._approximate_valuation
            m = self._constant[1]
            while n <= m:
                if self.coefficient(n) != 0:
                    self._approximate_valuation = n
                    return n
                n += 1
            return infinity

    def prec(self):
        """
        Return the precision of the series, which is infinity.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z)
            sage: f.prec()
            +Infinity
        """
        return infinity

    def polynomial(self, degree=None, name=None):
        """
        Return the polynomial or Laurent polynomial if the series is actually so.

        INPUT:

        - ``degree`` -- ``None`` or an integer

        - ``name`` -- name of the variable; if it is ``None``, the name of the variable
          of the series is used

        OUTPUT: a Laurent polynomial if the valuation of the series is negative or
        a polynomial otherwise.

        If ``degree`` is not ``None``, the terms of the series of degree
        greater than ``degree`` are truncated first. If ``degree`` is ``None``
        and the series is not a polynomial or a Laurent polynomial, a
        ``ValueError`` is raised.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,0,0,2,0,0,0,3], 5); f
            z^5 + 2*z^8 + 3*z^12
            sage: f.polynomial()
            3*z^12 + 2*z^8 + z^5

        ::

            sage: g = L.series([1,0,0,2,0,0,0,3], -5); g
            z^-5 + 2*z^-2 + 3*z^2
            sage: g.polynomial()
            z^-5 + 2*z^-2 + 3*z^2

        ::

            sage: z = L.gen()
            sage: f = (1 + z)/(z^3 - z^5)
            sage: f
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + ...
            sage: f.polynomial(5)
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + z^4 + z^5
            sage: f.polynomial(0)
            z^-3 + z^-2 + z^-1 + 1
            sage: f.polynomial(-5)
            0
        """
        if degree is None:
            if self._constant is None or not self._constant[0].is_zero():
                raise ValueError("not a polynomial")
            m = self._constant[1]
        else:
            m = degree + 1

        S = self.parent()

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            from sage.rings.all import LaurentPolynomialRing
            R = LaurentPolynomialRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self.coefficient(i) for i in range(n,m)]).shift(n)
        else:
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(S.base_ring(), name=name)
            return R([self.coefficient(i) for i in range(m)])

    def approximate_series(self, prec, name=None):
        """
        Return the Laurent series with absolute precision ``prec`` approximated
        from this series.

        INPUT:

        - ``prec`` -- an integer

        - ``name`` -- name of the variable; if it is ``None``, the name of the variable
          of the series is used

        OUTPUT: a Laurent series with absolute precision ``prec``

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: f = (z - 2*z^3)^5/(1 - 2*z)
            sage: f
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + 32*z^10 - 16*z^11 + ...
            sage: g = f.approximate_series(10)
            sage: g
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + O(z^10)
            sage: g.parent()
            Power Series Ring in z over Integer Ring

        ::

            sage: h = (f^-1).approximate_series(3)
            sage: h
            z^-5 - 2*z^-4 + 10*z^-3 - 20*z^-2 + 60*z^-1 - 120 + 280*z - 560*z^2 + O(z^3)
            sage: h.parent()
            Laurent Series Ring in z over Integer Ring
        """
        S = self.parent()

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            from sage.rings.all import LaurentSeriesRing
            R = LaurentSeriesRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self.coefficient(i) for i in range(n,prec)], n).add_bigoh(prec)
        else:
            from sage.rings.all import PowerSeriesRing
            R = PowerSeriesRing(S.base_ring(), name=name)
            return R([self.coefficient(i) for i in range(prec)]).add_bigoh(prec)

    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
        """

        R = self.parent()

        op = LazyLaurentSeriesOperator_mul(self, other)

        a = self._approximate_valuation + other._approximate_valuation

        c = None
        if self._constant is not None and other._constant is not None:
            if self._constant[0] == 0 and other._constant[0] == 0:
                c = (self._constant[0], self._constant[1] + other._constant[1] - 1)

        return R.element_class(R, coefficient=op, valuation=a, constant=c)

    def _rmul_(self, scalar):
        """
        Return the scalar multiplication of this series by ``scalar``.

        INPUT:

        - ``scalar`` -- an element of the base ring

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: 2*z
            2*z
            sage: -1*z
            -z
            sage: 0*z
            0
        """
        R = self.parent()

        if scalar.is_zero():
            return R.zero()

        op = LazyLaurentSeriesOperator_scale(self, scalar)

        a = self._approximate_valuation

        if self._constant is not None:
            c = (scalar * self._constant[0], self._constant[1])
        else:
            c = None

        return R.element_class(R, coefficient=op, valuation=a, constant=c)

    def _add_(self, other):
        """
        Return the sum of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
        """
        R = self.parent()

        op = LazyLaurentSeriesOperator_add(self, other)

        a = min(self._approximate_valuation, other._approximate_valuation)

        if self._constant is not None and other._constant is not None:
            c = (self._constant[0] + other._constant[0],
                 max(self._constant[1], other._constant[1]))
        else:
            c = None

        return R.element_class(R, coefficient=op, valuation=a, constant=c)

    def _sub_(self, other):
        """
        Return the series of this series minus ``other`` series.

        INPUT:

        - ``other`` -- other series

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: z - z
            0
        """
        R = self.parent()

        op = LazyLaurentSeriesOperator_sub(self, other)

        a = min(self._approximate_valuation, other._approximate_valuation)

        if self._constant is not None and other._constant is not None:
            c = (self._constant[0] - other._constant[0],
                 max(self._constant[1], other._constant[1]))
        else:
            c = None

        return R.element_class(R, coefficient=op, valuation=a, constant=c)

    def _neg_(self):
        """
        Return the negative of this series.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: -(1 - z)
            -1 + z
        """
        R = self.parent()

        op = LazyLaurentSeriesOperator_neg(self)

        a = self._approximate_valuation

        if self._constant is not None:
            c = (-self._constant[0], self._constant[1])
        else:
            c = None

        return R.element_class(R, coefficient=op, valuation=a, constant=c)

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: ~(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        v = self.valuation()

        if v == infinity:
            raise ZeroDivisionError('cannot invert zero')

        R = self.parent()

        op = LazyLaurentSeriesOperator_inv(self)

        return R.element_class(R, coefficient=op, valuation=-v, constant=None)

    def __pow__(self, n):
        """
        Return the `n`-th power of the series.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: (1 - z)^-1
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: (1 - z)^0
            1
            sage: (1 - z)^3
            1 - 3*z + 3*z^2 - z^3
            sage: (1 - z)^-3
            1 + 3*z + 6*z^2 + 10*z^3 + 15*z^4 + 21*z^5 + 28*z^6 + ...
        """
        if n == 0:
            return self.parent().one()

        return generic_power(self, n)

    def _div_(self, other):
        """
        Return ``self`` divided by ``other``.

        INPUT:

        - ``other`` -- nonzero series

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: z/(1 - z)
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
        """
        R = self.parent()

        if other.is_zero():
            raise ZeroDivisionError("division by zero series")

        if self.is_zero():
            return R.zero()

        op = LazyLaurentSeriesOperator_div(self, other)

        a = self.valuation() - other.valuation()

        return R.element_class(R, coefficient=op, valuation=a, constant=None)

    def apply_to_coefficients(self, function):
        """
        Return the series with ``function`` applied to each coefficient of this series.

        INPUT:

        - ``function`` -- Python function

        Python function ``function`` returns a new coefficient for input coefficient.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: s = z/(1 - 2*z)
            sage: t = s.apply_to_coefficients(lambda c: c + 1)
            sage: s
            z + 2*z^2 + 4*z^3 + 8*z^4 + 16*z^5 + 32*z^6 + 64*z^7 + ...
            sage: t
            2*z + 3*z^2 + 5*z^3 + 9*z^4 + 17*z^5 + 33*z^6 + 65*z^7 + ...
        """
        R = self.parent()

        op = LazyLaurentSeriesOperator_apply(self, function)

        a = self._approximate_valuation

        if self._constant:
            c = (function(self._constant[0]), self._constant[1])
        else:
            c = None

        return R.element_class(R, coefficient=op, valuation=a, constant=c)

    def change_ring(self, ring):
        """
        Return this series with coefficients converted to elements of ``ring``.

        INPUT:

        - ``ring`` -- a ring

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: s = 2 + z
            sage: t = s.change_ring(QQ)
            sage: t^-1
            1/2 - 1/4*z + 1/8*z^2 - 1/16*z^3 + 1/32*z^4 - 1/64*z^5 + 1/128*z^6 + ...
        """
        R = self.parent()

        op = LazyLaurentSeriesOperator_change_ring(self, ring)

        a = self._approximate_valuation

        if self._constant:
            c = (ring(self._constant[0]), self._constant[1])
        else:
            c = None

        from .lazy_laurent_series_ring import LazyLaurentSeriesRing
        Q = LazyLaurentSeriesRing(ring, names=R.variable_name())
        return Q.element_class(Q, coefficient=op, valuation=a, constant=c)

    def truncate(self, d):
        """
        Return this series with its terms of degree >= ``d`` truncated.

        INPUT:

        - ``d`` -- integer

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: alpha = 1/(1-z)
            sage: alpha
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: beta = alpha.truncate(5)
            sage: beta
            1 + z + z^2 + z^3 + z^4
            sage: alpha - beta
            z^5 + z^6 + z^7 + z^8 + z^9 + z^10 + z^11 + ...
        """
        R = self.parent()

        op = LazyLaurentSeriesOperator_truncate(self, d)
        a = self._approximate_valuation
        c = (self.base_ring().zero(), d)

        return R.element_class(R, coefficient=op, valuation=a, constant=c)
