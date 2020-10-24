r"""
Lazy Laurent Series Operators

This module implements operators internally used to construct lazy Laurent
series. The job of an operator attached to a series is to compute the `n`-th
coefficient of the series if `n` is not less than the valuation of the series
and the `n`-th coefficient is not declared to be a constant.

If a new operator is added to this module, an example of how it is used should be
added below.

EXAMPLES::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: f = 1/(1 - 2*z)
    sage: g = 1/(1 + z^2)

Constructors::

    sage: L(1)
    1

::

    sage: L.series([1,2,3,4], -10)
    z^-10 + 2*z^-9 + 3*z^-8 + 4*z^-7

::

    sage: L.gen()
    z

::

    sage: P.<x> = LaurentPolynomialRing(ZZ)
    sage: p = (1 + 1/x)^3 + (1 + x)^4
    sage: L(p)
    z^-3 + 3*z^-2 + 3*z^-1 + 2 + 4*z + 6*z^2 + 4*z^3 + z^4

Unary operators::

    sage: -f
    -1 - 2*z - 4*z^2 - 8*z^3 - 16*z^4 - 32*z^5 - 64*z^6 + ...

::

    sage: ~f
    1 - 2*z + ...

Binary operators::

    sage: f + g
    2 + 2*z + 3*z^2 + 8*z^3 + 17*z^4 + 32*z^5 + 63*z^6 + ...

::

    sage: f - g
    2*z + 5*z^2 + 8*z^3 + 15*z^4 + 32*z^5 + 65*z^6 + 128*z^7 + ...

::

    sage: f * g
    1 + 2*z + 3*z^2 + 6*z^3 + 13*z^4 + 26*z^5 + 51*z^6 + ...

::

    sage: f / g
    1 + 2*z + 5*z^2 + 10*z^3 + 20*z^4 + 40*z^5 + 80*z^6 + ...

Transformers::

    sage: 2*f
    2 + 4*z + 8*z^2 + 16*z^3 + 32*z^4 + 64*z^5 + 128*z^6 + ...

::

    sage: f.change_ring(GF(3))
    1 + 2*z + z^2 + 2*z^3 + z^4 + 2*z^5 + z^6 + ...

::

    sage: f.apply_to_coefficients(lambda c: c^2)
    1 + 4*z + 16*z^2 + 64*z^3 + 256*z^4 + 1024*z^5 + 4096*z^6 + ...

::

    sage: f.truncate(5)
    1 + 2*z + 4*z^2 + 8*z^3 + 16*z^4

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


class LazyLaurentSeriesOperator(object):
    """
    Base class for operators computing coefficients of a lazy Laurent series.

    Subclasses of this class are used to implement arithmetic operations for
    lazy Laurent series. These classes are not to be used directly by the user.
    """
    def __ne__(self, other):
        """
        Test inequality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f != g
            False
        """
        return not (self == other)


class LazyLaurentSeriesBinaryOperator(LazyLaurentSeriesOperator):
    """
    Abstract base class for binary operators.

    INPUT:

    - ``left`` -- series on the left side of the binary operator

    - ``right`` -- series on the right side of the binary operator

    """
    def __init__(self, left, right):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
            sage: f = 1/(1 - z) - 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        self._left = left
        self._right = right

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: {f: 1}
            {2 + 2*z^2 + 2*z^4 + 2*z^6 + ...: 1}
        """
        return hash((type(self), self._left, self._right))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
        """
        return (isinstance(other, type(self)) and
                self._left == other._left and self._right == other._right)


class LazyLaurentSeriesUnaryOperator(LazyLaurentSeriesOperator):
    """
    Abstract base class for unary operators.

    INPUT:

    - ``series`` -- series upon which the operator operates

    """
    def __init__(self, series):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = -1/(1 - z)
            sage: f
            -1 - z - z^2 - z^3 - z^4 - z^5 - z^6 + ...
            sage: loads(dumps(f)) == f
            True
        """
        self._series = series

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - z)
            sage: {f: 1}
            {1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
            sage: f = ~(1 - z)
            sage: g = ~(1 - z)
            sage: f == g
            True
        """
        return isinstance(other, type(self)) and self._series == other._series


class LazyLaurentSeriesOperator_add(LazyLaurentSeriesBinaryOperator):
    """
    Operator for addition.
    """
    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = z + z^2 + z
            sage: f.coefficient(1)
            2
        """
        return self._left.coefficient(n) + self._right.coefficient(n)

class LazyLaurentSeriesOperator_sub(LazyLaurentSeriesBinaryOperator):
    """
    Operator for subtraction.
    """
    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1 + 3*z - z
            sage: f.coefficient(1)
            2
        """
        return self._left.coefficient(n) - self._right.coefficient(n)

class LazyLaurentSeriesOperator_mul(LazyLaurentSeriesBinaryOperator):
    """
    Operator for multiplication.
    """
    def __init__(self, left, right):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) * 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        LazyLaurentSeriesBinaryOperator.__init__(self, left, right)
        self._zero = left.base_ring().zero()

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = (1 + z)*(1 - z)
            sage: f.coefficient(2)
            -1
        """
        c = self._zero
        for k in range(self._left._approximate_valuation, n - self._right._approximate_valuation + 1):
            c += self._left.coefficient(k) * self._right.coefficient(n-k)
        return c

class LazyLaurentSeriesOperator_neg(LazyLaurentSeriesUnaryOperator):
    """
    Operator for negation.

    INPUT:

    - ``series`` -- a lazy Laurent series

    """
    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = -(1 + z)
            sage: f.coefficient(1)
            -1
        """
        return -self._series.coefficient(n)

class LazyLaurentSeriesOperator_inv(LazyLaurentSeriesUnaryOperator):
    """
    Operator for inversion.

    INPUT:

    - ``series`` -- a lazy Laurent series

    """
    def __init__(self, series):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - z)
            sage: f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: loads(dumps(f)) == f
            True
        """
        LazyLaurentSeriesUnaryOperator.__init__(self, series)

        self._v = series.valuation()
        self._ainv = ~series.coefficient(self._v)
        self._zero = series.base_ring().zero()

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - z)
            sage: f.coefficient(2)
            1
        """
        v = self._v
        if n == -v:
            return self._ainv
        c = self._zero
        for k in range(-v, n):
            c += s.coefficient(k) * self._series.coefficient(n + v - k)
        return -c * self._ainv

class LazyLaurentSeriesOperator_div(LazyLaurentSeriesBinaryOperator):
    """
    Operator for division.
    """
    def __init__(self, left, right):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = (1 - z)/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        LazyLaurentSeriesBinaryOperator.__init__(self, left, right)

        lv = left.valuation()
        rv = right.valuation()

        self._lv = lv
        self._rv = rv
        self._ainv = ~right.coefficient(rv)

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = (1 + z)/(1 - z)
            sage: f.coefficient(2)
            2
        """
        lv = self._lv
        rv = self._rv

        if n == lv - rv:
            return self._left.coefficient(lv)/self._right.coefficient(rv)
        c = self._left.coefficient(n + rv)
        for k in range(lv - rv, n):
            c -= s.coefficient(k) * self._right.coefficient(n + rv - k)
        return c * self._ainv

class LazyLaurentSeriesOperator_scale(LazyLaurentSeriesOperator):
    """
    Operator for scalar multiplication of ``series`` with ``scalar``.

    INPUT:

    - ``series`` -- a lazy Laurent series

    - ``scalar`` -- an element of the base ring

    """
    def __init__(self, series, scalar):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: g = 2*z
            sage: loads(dumps(g)) == g
            True
        """
        self._series = series
        self._scalar = scalar

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 2*(z + z^2)
            sage: f.coefficient(2)
            2
            sage: f
            2*z + 2*z^2
        """
        return self._scalar * self._series.coefficient(n)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 2*(z + z^2)
            sage: {f: 1}
            {2*z + 2*z^2: 1}
        """
        return hash((type(self), self._series, self._scalar))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 2*z
            sage: g = 2*z
            sage: f == g
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_scale)
                and self._series == other._series and self._scalar == other._scalar)

class LazyLaurentSeriesOperator_change_ring(LazyLaurentSeriesOperator):
    """
    Operator for changing the base ring of the ``series`` to ``ring``.

    INPUT:

    - ``series`` -- a lazy Laurent series

    - ``ring`` -- a ring

    """
    def __init__(self, series, ring):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - 2*z)
            sage: g = f.change_ring(GF(3))
            sage: g
            1 + 2*z + z^2 + 2*z^3 + z^4 + 2*z^5 + z^6 + ...
            sage: loads(dumps(g)) == g
            True
        """
        self._series = series
        self._ring = ring

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - 2*z)
            sage: f
            1 + 2*z + 4*z^2 + 8*z^3 + 16*z^4 + 32*z^5 + 64*z^6 + ...
            sage: g = f.change_ring(GF(2))
            sage: g.coefficient(3)
            0
            sage: g
            1 + ...
        """
        return self._ring(self._series.coefficient(n))

    def __hash__(self):
        """
        Return the hash of ``self``.

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - 2*z)
            sage: g = f.change_ring(GF(2))
            sage: {g: 1}
            {1 + ...: 1}
        """
        return hash((type(self), self._series, self._ring))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z)
            sage: g = 1/(1 - z)
            sage: f == g
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_change_ring)
                and self._series == other._series and self._ring == other._ring)

class LazyLaurentSeriesOperator_apply(LazyLaurentSeriesOperator):
    """
    Operator for applying a function.

    INPUT:

    - ``series`` -- a lazy Laurent series

    - ``function`` -- a Python function to apply to each coefficient of the series

    """
    def __init__(self, series, function):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 + z)
            sage: g = f.apply_to_coefficients(abs)
            sage: g
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: loads(dumps(g))
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        self._series = series
        self._function = function

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 + z)
            sage: f
            1 - z + z^2 - z^3 + z^4 - z^5 + z^6 + ...
            sage: f.apply_to_coefficients(lambda c: c if c >= 0 else -c)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        return self._function(self._series.coefficient(n))

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 + z).apply_to_coefficients(lambda c: c if c >= 0 else -c)
            sage: {f: 1}
            {1 - z + z^2 - z^3 + z^4 - z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series, self._function))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 + z)
            sage: g = f.apply_to_coefficients(abs)
            sage: h = f.apply_to_coefficients(abs)
            sage: g == h
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_apply)
                and self._series == other._series and self._function is other._function)


class LazyLaurentSeriesOperator_truncate(LazyLaurentSeriesOperator):
    """
    Operator for truncation.

    INPUT:

    - ``series`` -- a lazy Laurent series

    - ``d`` -- an integer; the series is truncated the terms of degree `> d`

    """
    def __init__(self, series, d):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 + z)
            sage: g = f.truncate(4)
            sage: loads(dumps(g)) == g
            True
        """
        self._series = series
        self._d = d

        self._zero = series.base_ring().zero()

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 + z)
            sage: f
            1 - z + z^2 - z^3 + z^4 - z^5 + z^6 + ...
            sage: f.truncate(4)
            1 - z + z^2 - z^3
        """
        if n <= self._d:
            return self._series.coefficient(n)
        else:
            return self._zero

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 + z)
            sage: {f: 1}
            {1 - z + z^2 - z^3 + z^4 - z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series, self._d))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 + z)
            sage: g = f.truncate(4)
            sage: h = f.truncate(4)
            sage: g == h
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_truncate)
                and self._series == other._series and self._d == other._d)

class LazyLaurentSeriesOperator_gen(LazyLaurentSeriesOperator):
    """
    Operator for the generator element.

    INPUT:

    - ``ring`` -- a lazy Laurent series ring

    """
    def __init__(self, ring):
        """
        Initialize.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: loads(dumps(z)) == z
            True
        """
        self._ring = ring

        self._one = ring.base_ring().one()
        self._zero = ring.base_ring().zero()

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.gen()
            z
        """
        return self._one if n == 1 else self._zero

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: {z: 1}
            {z: 1}
        """
        return hash((type(self), self._ring))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z1 = L.gen()
            sage: z2 = L.gen()
            sage: z1 == z2
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_gen) and self._ring == other._ring)

class LazyLaurentSeriesOperator_constant(LazyLaurentSeriesOperator):
    """
    Operator for the generator element.

    INPUT:

    - ``ring`` -- a lazy Laurent series ring

    - ``constant`` -- a constant of the base ring of ``ring``

    """
    def __init__(self, ring, constant):
        """
        Initialize.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L(10)
            sage: loads(dumps(f)) == f
            True
        """
        self._ring = ring
        self._constant = constant

        self._zero = ring.base_ring().zero()

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L(10); f
            10
            sage: f.coefficient(0)
            10
            sage: f.coefficient(1)
            0
        """
        return self._constant if n == 0 else self._zero

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L(10)
            sage: {f: 1}
            {10: 1}
        """
        return hash((type(self), self._ring, self._constant))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z1 = L(10)
            sage: z2 = L(10)
            sage: z1 == z2
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_constant)
                and self._ring == other._ring and self._constant == other._constant)

class LazyLaurentSeriesOperator_list(LazyLaurentSeriesOperator):
    """
    Operator for the series defined by a list.

    INPUT:

    - ``l`` -- list

    - ``v`` -- integer

    """
    def __init__(self, ring, l, v):
        """
        Initialize.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5)
            sage: loads(dumps(f)) == f
            True
        """
        self._ring = ring
        self._list = tuple([ring.base_ring()(e) for e in l])
        self._valuation = v

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5)
            sage: f
            z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2
        """
        return self._ring.base_ring()(self._list[n - self._valuation])

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5)
            sage: {f: 1}
            {z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2: 1}
        """
        return hash((type(self), self._ring, self._list, self._valuation))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f1 = L.series([1,2,3,4], -5)
            sage: f2 = L.series([1,2,3,4,0], -5)
            sage: f1 == f2
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_list) and
                self._ring == other._ring and self._list == other._list and
                self._valuation == other._valuation)

class LazyLaurentSeriesOperator_polynomial(LazyLaurentSeriesOperator):
    """
    Operator for the series coerced from a polynomial or a Laurent polynomial.

    INPUT:

    - ``ring`` -- a lazy Laurent series ring

    - ``poly`` -- a polynomial or a Laurent polynomial

    """
    def __init__(self, ring, poly):
        """
        Initialize.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: P.<x> = ZZ[]
            sage: p = 1 + 2*x + x^10
            sage: f = L(p)
            sage: loads(dumps(f)) == f
            True
        """
        self._ring = ring
        self._poly = poly

    def __call__(self, s, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: P.<x> = ZZ[]
            sage: p = (1 + 2*x + 3*x)^3
            sage: f = L(p)
            sage: f
            1 + 15*z + 75*z^2 + 125*z^3
        """
        return self._poly[n]

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: P.<x> = ZZ[]
            sage: p = (1 + 2*x)^3
            sage: f = L(p)
            sage: {f: 1}
            {1 + 6*z + 12*z^2 + 8*z^3: 1}
        """
        return hash((type(self), self._ring, self._poly))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: P.<x> = ZZ[]
            sage: p = (1 + 2*x)^3
            sage: f = L(p)
            sage: g = L(p)
            sage: f == g
            True
        """
        return (isinstance(other, LazyLaurentSeriesOperator_polynomial) and
                self._ring == other._ring and self._poly == other._poly)
