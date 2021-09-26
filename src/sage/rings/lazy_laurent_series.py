r"""
Lazy Laurent Series

A lazy Laurent series is a Laurent series whose coefficients are
computed on demand.  Therefore, unlike the usual Laurent series in
Sage, lazy Laurent series have infinite precision.

EXAMPLES:

Laurent series over the integer ring are particularly useful as
generating functions for sequences arising in combinatorics. ::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)

The generating function of the Fibonacci sequence is::

    sage: f = 1 / (1 - z - z^2)
    sage: f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + O(z^7)

In principle, we can now compute any coefficient of `f`::

    sage: f.coefficient(100)
    573147844013817084101

Which coefficients are actually computed depends on the type of
implementation.  For the sparse implementation, only the coefficients
which are needed are computed. ::

    sage: s = L(lambda n: n, valuation=0); s
    z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
    sage: s.coefficient(10)
    10
    sage: s._coeff_stream._cache
    {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 10: 10}

Using the dense implementation, all coefficients up to the required
coefficient are computed. ::

    sage: L.<x> = LazyLaurentSeriesRing(ZZ, sparse=False)
    sage: s = L(lambda n: n, valuation=0); s
    x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5 + 6*x^6 + O(x^7)
    sage: s.coefficient(10)
    10
    sage: s._coeff_stream._cache
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

We can do arithmetic with lazy power series::

    sage: f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + O(z^7)
    sage: f^-1
    1 - z - z^2 + O(z^7)
    sage: f + f^-1
    2 + z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + O(z^7)
    sage: g = (f + f^-1)*(f - f^-1); g
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + O(z^7)

We can change the base ring::

    sage: h = g.change_ring(QQ)
    sage: h.parent()
    Lazy Laurent Series Ring in z over Rational Field
    sage: h
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + O(z^7)
    sage: hinv = h^-1; hinv
    1/4*z^-1 - 3/8 + 1/16*z - 17/32*z^2 + 5/64*z^3 - 29/128*z^4 + 165/256*z^5 + O(z^6)
    sage: hinv.valuation()
    -1

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version
- Tejasvi Chebrolu, Martin Rubey, Travis Scrimshaw (2021-08):
  refactored and expanded functionality

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


from sage.rings.infinity import infinity
from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.structure.richcmp import op_EQ, op_NE
from sage.arith.power import generic_power
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.data_structures.stream import (
    Stream_add,
    Stream_cauchy_mul,
    Stream_sub,
    Stream_cauchy_compose,
    Stream_lmul,
    Stream_rmul,
    Stream_neg,
    Stream_cauchy_invert,
    Stream_map_coefficients,
    Stream_zero,
    Stream_exact,
    Stream_uninitialized,
    Stream_shift
)

class LazyModuleElement(Element):
    r"""
    A lazy sequence with a module structure given by term-wise
    addition and scalar multiplication.

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: M = L(lambda n: n, valuation=0)
        sage: N = L(lambda n: 1, valuation=0)
        sage: M[:10]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: N[:10]
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    Two sequences can be added::

        sage: O = M + N
        sage: O[0:10]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    Two sequences can be subtracted::

        sage: P = M - N
        sage: P[:10]
        [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]

    A sequence can be multiplied by a scalar::

        sage: Q = 2 * M
        sage: Q[:10]
        [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

    The negation of a sequence can also be found::

        sage: R = -M
        sage: R[:10]
        [0, -1, -2, -3, -4, -5, -6, -7, -8, -9]
    """
    def __init__(self, parent, coeff_stream):
        """
        Initialize the series.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: TestSuite(L.an_element()).run()

        """
        Element.__init__(self, parent)
        self._coeff_stream = coeff_stream

    def __getitem__(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer; the exponent

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: f = z / (1 - 2*z^3)
            sage: [f[n] for n in range(20)]
            [0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 16, 0, 0, 32, 0, 0, 64]
            sage: f[0:20]
            [0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 16, 0, 0, 32, 0, 0, 64]

            sage: M = L(lambda n: n, valuation=0)
            sage: [M[n] for n in range(20)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: n, valuation=0)
            sage: [M[n] for n in range(20)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        """
        R = self.parent()._coeff_ring
        if isinstance(n, slice):
            if n.stop is None:
                raise NotImplementedError("cannot list an infinite set")
            start = n.start if n.start is not None else self._coeff_stream._approximate_order
            step = n.step if n.step is not None else 1
            return [R(self._coeff_stream[k]) for k in range(start, n.stop, step)]
        return R(self._coeff_stream[n])

    coefficient = __getitem__

    def map_coefficients(self, func, ring=None):
        r"""
        Return the series with ``func`` applied to each nonzero
        coefficient of ``self``.

        INPUT:

        - ``func`` -- function that takes in a coefficient and returns
          a new coefficient

        EXAMPLES:

        Dense Implementation::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: s = z/(1 - 2*z^2)
            sage: t = s.map_coefficients(lambda c: c + 1)
            sage: s
            z + 2*z^3 + 4*z^5 + 8*z^7 + O(z^8)
            sage: t
            2*z + 3*z^3 + 5*z^5 + 9*z^7 + O(z^8)
            sage: m = L(lambda n: n, valuation=0); m
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: m.map_coefficients(lambda c: c + 1)
            2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + O(z^7)

        Sparse Implementation::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: m = L(lambda n: n, valuation=0); m
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: m.map_coefficients(lambda c: c + 1)
            2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + O(z^7)

        An example where the series is known to be exact::

            sage: f = z + z^2 + z^3
            sage: f.map_coefficients(lambda c: c + 1)
            2*z + 2*z^2 + 2*z^3

        TESTS::

            sage: from sage.data_structures.stream import Stream_zero
            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: s = L(0).map_coefficients(lambda c: c + 1); s
            0
            sage: isinstance(s._coeff_stream, Stream_zero)
            True

        """
        P = self.parent()
        coeff_stream = self._coeff_stream
        if isinstance(coeff_stream, Stream_zero):
            return self
        BR = P.base_ring()
        if isinstance(coeff_stream, Stream_exact):
            initial_coefficients = [func(i) if i else 0
                                    for i in coeff_stream._initial_coefficients]
            c = func(coeff_stream._constant) if coeff_stream._constant else 0
            if not any(initial_coefficients) and not c:
                return P.zero()
            coeff_stream = Stream_exact(initial_coefficients,
                                                   self._coeff_stream._is_sparse,
                                                   order=coeff_stream._approximate_order,
                                                   degree=coeff_stream._degree,
                                                   constant=BR(c))
            return P.element_class(P, coeff_stream)
        coeff_stream = Stream_map_coefficients(self._coeff_stream, func, P._coeff_ring)
        return P.element_class(P, coeff_stream)

    def truncate(self, d):
        r"""
        Return this series with its terms of degree >= ``d`` truncated.

        INPUT:

        - ``d`` -- integer; the degree from which the series is truncated

        EXAMPLES:

        Dense Implementation::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: alpha = 1/(1-z)
            sage: alpha
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: beta = alpha.truncate(5)
            sage: beta
            1 + z + z^2 + z^3 + z^4
            sage: alpha - beta
            z^5 + z^6 + O(z^7)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: M.truncate(4)
            z + 2*z^2 + 3*z^3

        Sparse Implementation::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: M.truncate(4)
            z + 2*z^2 + 3*z^3

        Series which are known to be exact can also be truncated::

            sage: M = z + z^2 + z^3 + z^4
            sage: M.truncate(4)
            z + z^2 + z^3
        """
        P = self.parent()
        coeff_stream = self._coeff_stream
        v = coeff_stream._approximate_order
        initial_coefficients = [coeff_stream[i] for i in range(v, d)]
        return P.element_class(P, Stream_exact(initial_coefficients, P._sparse,
                                                          order=v))

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

    def valuation(self):
        r"""
        Return the valuation of ``self``.

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
            sage: M = L(lambda n: n^2, valuation=0)
            sage: M.valuation()
            1
            sage: (M - M).valuation()
            +Infinity
        """
        return self._coeff_stream.order()

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` with ``other`` with respect to the comparison
        operator ``op``.

        Equality is verified if the corresponding coefficients of both series
        can be checked for equality without computing coefficients
        indefinitely.  Otherwise an exception is raised to declare that
        equality is not decidable.

        Inequality is not defined for lazy Laurent series.

        INPUT:

        - ``other`` -- another Laurent series
        - ``op`` -- comparison operator

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: z + z^2 == z^2 + z
            True
            sage: z + z^2 != z^2 + z
            False
            sage: z + z^2 > z^2 + z
            False
            sage: z + z^2 < z^2 + z
            False

            sage: fz = L(lambda n: 0, valuation=0)
            sage: L.zero() == fz
            Traceback (most recent call last):
            ...
            ValueError: undecidable
            sage: fz == L.zero()
            Traceback (most recent call last):
            ...
            ValueError: undecidable
        """
        if op is op_EQ:
            if isinstance(self._coeff_stream, Stream_zero):  # self == 0
                if isinstance(other._coeff_stream, Stream_zero):
                    return True
                if other._coeff_stream.is_nonzero():
                    return False
            # other == 0 but self likely != 0
            elif (isinstance(other._coeff_stream, Stream_zero)
                  and self._coeff_stream.is_nonzero()):
                return False

            if (not isinstance(self._coeff_stream, Stream_exact)
                or not isinstance(other._coeff_stream, Stream_exact)):
                # One of the lazy laurent series is not known to eventually be constant
                # Implement the checking of the caches here.
                n = min(self._coeff_stream._approximate_order, other._coeff_stream._approximate_order)
                m = max(self._coeff_stream._approximate_order, other._coeff_stream._approximate_order)
                for i in range(n, m):
                    if self[i] != other[i]:
                        return False
                if self._coeff_stream == other._coeff_stream:
                    return True
                if self._coeff_stream != other._coeff_stream:
                    return False
                raise ValueError("undecidable")

            # Both are Stream_exact, which implements a full check
            return self._coeff_stream == other._coeff_stream

        if op is op_NE:
            return not (self == other)

        return False

    def __hash__(self):
        """
        Return the hash of ``self``

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L([1,2,3,4], valuation=-5)
            sage: hash(f) == hash(f)
            True
            sage: g = (1 + f)/(1 - f)^2
            sage: {g: 1}
            {z^5 - 2*z^6 + z^7 + 5*z^9 - 11*z^10 + z^11 + O(z^12): 1}
        """
        return hash(self._coeff_stream)

    def __bool__(self):
        """
        Test whether ``self`` is not zero.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(GF(2))
            sage: bool(z-z)
            False
            sage: f = 1/(1 - z)
            sage: bool(f)
            True
            sage: M = L(lambda n: n, valuation=0); M
            z + z^3 + z^5 + O(z^7)
            sage: M.is_zero()
            False
            sage: M = L(lambda n: 2*n if n < 10 else 1, valuation=0); M
            O(z^7)
            sage: bool(M)
            Traceback (most recent call last):
            ...
            ValueError: undecidable as lazy Laurent series
            sage: M[15]
            1
            sage: bool(M)
            True

            sage: L.<z> = LazyLaurentSeriesRing(GF(2), sparse=True)
            sage: M = L(lambda n: 2*n if n < 10 else 1, valuation=0); M
            O(z^7)
            sage: bool(M)
            Traceback (most recent call last):
            ...
            ValueError: undecidable as lazy Laurent series
            sage: M[15]
            1
            sage: bool(M)
            True
        """
        if isinstance(self._coeff_stream, Stream_zero):
            return False
        if isinstance(self._coeff_stream, Stream_exact):
            return True
        if self.parent()._sparse:
            cache = self._coeff_stream._cache
            if any(cache[a] for a in cache):
                return True
        else:
            if any(a for a in self._coeff_stream._cache):
                return True
        if self[self._coeff_stream._approximate_order]:
            return True
        raise ValueError("undecidable as lazy Laurent series")

    def define(self, s):
        r"""
        Define an equation by ``self = s``.

        INPUT:

        - ``s`` -- a Laurent polynomial

        EXAMPLES:

        We begin by constructing the Catalan numbers::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: C = L(None, valuation=0)
            sage: C.define(1 + z*C^2)
            sage: C
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + O(z^7)

        The Catalan numbers but with a valuation 1::

            sage: B = L(None, valuation=1)
            sage: B.define(z + B^2)
            sage: B
            z + z^2 + 2*z^3 + 5*z^4 + 14*z^5 + 42*z^6 + 132*z^7 + O(z^8)

        We can define multiple series that are linked::

            sage: s = L(None, valuation=0)
            sage: t = L(None, valuation=0)
            sage: s.define(1 + z*t^3)
            sage: t.define(1 + z*s^2)
            sage: s[:9]
            [1, 1, 3, 9, 34, 132, 546, 2327, 10191]
            sage: t[:9]
            [1, 1, 2, 7, 24, 95, 386, 1641, 7150]

        A bigger example::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: A = L(None, valuation=5)
            sage: B = L(None, valuation=0)
            sage: C = L(None, valuation=2)
            sage: A.define(z^5 + B^2)
            sage: B.define(z^5 + C^2)
            sage: C.define(z^2 + C^2 + A^2)
            sage: A[0:15]
            [0, 0, 0, 0, 0, 1, 0, 0, 1, 2, 5, 4, 14, 10, 48]
            sage: B[0:15]
            [0, 0, 0, 0, 1, 1, 2, 0, 5, 0, 14, 0, 44, 0, 138]
            sage: C[0:15]
            [0, 0, 1, 0, 1, 0, 2, 0, 5, 0, 15, 0, 44, 2, 142]

        Counting binary trees::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: s = L(None, valuation=1)
            sage: s.define(z + (s^2+s(z^2))/2)
            sage: [s[i] for i in range(9)]
            [0, 1, 1, 1, 2, 3, 6, 11, 23]

        The `q`-Catalan numbers::

            sage: R.<q> = ZZ[]
            sage: L.<z> = LazyLaurentSeriesRing(R)
            sage: s = L(None, valuation=0)
            sage: s.define(1+z*s*s(q*z))
            sage: s
            1 + z + (q + 1)*z^2 + (q^3 + q^2 + 2*q + 1)*z^3
             + (q^6 + q^5 + 2*q^4 + 3*q^3 + 3*q^2 + 3*q + 1)*z^4
             + (q^10 + q^9 + 2*q^8 + 3*q^7 + 5*q^6 + 5*q^5 + 7*q^4 + 7*q^3 + 6*q^2 + 4*q + 1)*z^5
             + (q^15 + q^14 + 2*q^13 + 3*q^12 + 5*q^11 + 7*q^10 + 9*q^9 + 11*q^8
                + 14*q^7 + 16*q^6 + 16*q^5 + 17*q^4 + 14*q^3 + 10*q^2 + 5*q + 1)*z^6 + O(z^7)

        We count unlabeled ordered trees by total number of nodes
        and number of internal nodes::

            sage: R.<q> = QQ[]
            sage: Q.<z> = LazyLaurentSeriesRing(R)
            sage: leaf = z
            sage: internal_node = q * z
            sage: L = Q(constant=1, degree=1)
            sage: T = Q(None, valuation=1)
            sage: T.define(leaf + internal_node * L(T))
            sage: [T[i] for i in range(6)]
            [0, 1, q, q^2 + q, q^3 + 3*q^2 + q, q^4 + 6*q^3 + 6*q^2 + q]

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: s = L(None, valuation=0)
            sage: s.define(1 + z*s^3)
            sage: s[:10]
            [1, 1, 3, 12, 55, 273, 1428, 7752, 43263, 246675]

            sage: e = L(None, valuation=0)
            sage: e.define(1 + z*e)
            sage: e.define(1 + z*e)
            Traceback (most recent call last):
            ...
            ValueError: series already defined
            sage: z.define(1 + z^2)
            Traceback (most recent call last):
            ...
            ValueError: series already defined
        """
        if not isinstance(self._coeff_stream, Stream_uninitialized) or self._coeff_stream._target is not None:
            raise ValueError("series already defined")
        self._coeff_stream._target = s._coeff_stream

    # === module structure ===

    def _add_(self, other):
        """
        Return the sum of ``self`` and ``other``.

        INPUT:

        - ``other`` -- other series

        EXAMPLES:

        Dense series can be added::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: m = L(lambda n: 1 + n, valuation=0)
            sage: n = L(lambda n: -n, valuation=0)
            sage: s = m + n
            sage: s[0:10]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        Sparse series can be added::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: m = L(lambda n: 1 + n, valuation=0)
            sage: n = L(lambda n: -n, valuation=0)
            sage: s = m + n
            sage: s[0:10]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        Series which are known to be exact can be added::

            sage: m = L(1)
            sage: n = L([0, 1])
            sage: s = m + n
            sage: s[0:10]
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]

        Adding zero gives the same series::

            sage: m = L(lambda n: 1 + n, valuation=0)
            sage: m + 0 is 0 + m is m
            True
        """
        P = self.parent()
        left = self._coeff_stream
        right = other._coeff_stream
        if isinstance(left, Stream_zero):
            return other
        if isinstance(right, Stream_zero):
            return self
        if (isinstance(left, Stream_exact)
            and isinstance(right, Stream_exact)):
            approximate_order = min(left.order(), right.order())
            degree = max(left._degree, right._degree)
            initial_coefficients = [left[i] + right[i] for i in range(approximate_order, degree)]
            constant = left._constant + right._constant
            if not any(initial_coefficients) and not constant:
                return P.zero()
            coeff_stream = Stream_exact(initial_coefficients, P._sparse,
                                                   constant=constant,
                                                   degree=degree,
                                                   order=approximate_order)
            return P.element_class(P, coeff_stream)
        return P.element_class(P, Stream_add(self._coeff_stream, other._coeff_stream))

    def _sub_(self, other):
        """
        Return the series of this series minus ``other`` series.

        INPUT:

        - ``other`` -- other series

        EXAMPLES:

        Dense series can be subtracted::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: m = L(lambda n: 1 + n, valuation=0)
            sage: n = L(lambda n: -n, valuation=0)
            sage: d = m - n
            sage: d[0:10]
            [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]

        Sparse series can be subtracted::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: m = L(lambda n: 1 + n, valuation=0)
            sage: n = L(lambda n: -n, valuation=0)
            sage: d = m - n
            sage: d[0:10]
            [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]

        Series which are known to be exact can be subtracted::

            sage: m = L.one()
            sage: n = L([0, 1])
            sage: d = m - n
            sage: d[0:10]
            [1, -1, 0, 0, 0, 0, 0, 0, 0, 0]

            sage: m = L([1, 0, 1])
            sage: n = L([0, 0, 1])
            sage: d = m - L.one() - n
            sage: d
            0

        Subtraction involving 0::

            sage: m = L(lambda n: 1 + n, valuation=0)
            sage: m - 0 is m
            True
            sage: 0 - m == -m
            True

            sage: A.<t> = LazyLaurentSeriesRing(QQ)
            sage: B.<z> = LazyLaurentSeriesRing(A)
            sage: 1 - z
            1 - z
        """
        right = other._coeff_stream
        if isinstance(right, Stream_zero):
            return self
        left = self._coeff_stream
        if isinstance(left, Stream_zero):
            return -other
        P = self.parent()
        if (isinstance(left, Stream_exact) and isinstance(right, Stream_exact)):
            approximate_order = min(left.order(), right.order())
            degree = max(left._degree, right._degree)
            initial_coefficients = [left[i] - right[i] for i in range(approximate_order, degree)]
            constant = left._constant - right._constant
            if not any(initial_coefficients) and not constant:
                return P.zero()
            coeff_stream = Stream_exact(initial_coefficients, P._sparse,
                                                   constant=constant,
                                                   degree=degree,
                                                   order=approximate_order)
            return P.element_class(P, coeff_stream)
        if left == right:
            return P.zero()
        return P.element_class(P, Stream_sub(self._coeff_stream, other._coeff_stream))

    def _acted_upon_(self, scalar, self_on_left):
        r"""
        Scalar multiplication for ``self`` by ``scalar``.

        INPUT:

        - ``scalar`` -- an element of the base ring
        - ``self_on_left`` -- boolean; if ``True``, compute ``self * scalar``

        EXAMPLES:

        Dense series can be multiplied with a scalar::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: M = L(lambda n: 1 + n, valuation=0)
            sage: O = M * 2
            sage: O[0:10]
            [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
            sage: type(O._coeff_stream)
            <class 'sage.data_structures.stream.Stream_lmul'>
            sage: M * 1 is M
            True
            sage: M * 0 == 0
            True
            sage: O = 2 * M
            sage: type(O._coeff_stream)
            <class 'sage.data_structures.stream.Stream_lmul'>
            sage: O[0:10]
            [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
            sage: 1 * M is M
            True
            sage: 0 * M == 0
            True

        Different scalars potentially give different series::

            sage: 2 * M == 3 * M
            Traceback (most recent call last):
            ...
            ValueError: undecidable

        Sparse series can be multiplied with a scalar::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: 1 + n, valuation=0)
            sage: O = M * 2
            sage: O[0:10]
            [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
            sage: type(O._coeff_stream)
            <class 'sage.data_structures.stream.Stream_lmul'>
            sage: M * 1 is M
            True
            sage: M * 0 == 0
            True
            sage: O = 2 * M
            sage: type(O._coeff_stream)
            <class 'sage.data_structures.stream.Stream_lmul'>
            sage: O[0:10]
            [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
            sage: 1 * M is M
            True
            sage: 0 * M == 0
            True

        Series which are known to be exact can be multiplied with a scalar
        and remain exact::

            sage: N = L([0, 1], degree=5, constant=3)
            sage: O = N * -1
            sage: O[0:10]
            [0, -1, 0, 0, 0, -3, -3, -3, -3, -3]
            sage: N * 1 is N
            True
            sage: N * 0 == 0
            True
            sage: O = -1 * N
            sage: O[0:10]
            [0, -1, 0, 0, 0, -3, -3, -3, -3, -3]
            sage: 1 * N is N
            True
            sage: 0 * N == 0
            True
        """
        # With the current design, the coercion model does not have
        # enough information to detect a priori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        P = self.parent()
        R = P.base_ring()
        if isinstance(scalar, Element) and scalar.parent() is not R:
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if R.has_coerce_map_from(scalar.parent()):
                scalar = R(scalar)
            else:
                return None

        coeff_stream = self._coeff_stream
        if isinstance(coeff_stream, Stream_zero):
            return self

        if not scalar:
            return P.zero()
        if scalar == R.one():
            return self
        if scalar == -R.one():
            return -self

        if isinstance(coeff_stream, Stream_exact):
            v = coeff_stream.order()
            init_coeffs = coeff_stream._initial_coefficients
            if self_on_left:
                c = coeff_stream._constant * scalar
                initial_coefficients = [val * scalar for val in init_coeffs]
            else:
                c = scalar * coeff_stream._constant
                initial_coefficients = [scalar * val for val in init_coeffs]
            return P.element_class(P, Stream_exact(initial_coefficients, P._sparse,
                                                   order=v, constant=c,
                                                   degree=coeff_stream._degree))
        if self_on_left or R.is_commutative():
            return P.element_class(P, Stream_lmul(coeff_stream, scalar))
        return P.element_class(P, Stream_rmul(coeff_stream, scalar))

    def _neg_(self):
        """
        Return the negative of ``self``.

        EXAMPLES:

        Dense series can be negated::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: m = L(lambda n: n, valuation=0)
            sage: n = L(lambda n: -n, valuation=0)
            sage: -n
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: -m
            -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 + O(z^7)
            sage: -(-m) == m
            True

        Sparse series can be negated::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: m = L(lambda n: n, valuation=0)
            sage: n = L(lambda n: -n, valuation=0)
            sage: -n
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: -m
            -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 + O(z^7)
            sage: -(-m) == m
            True

        The negation of an exact series is exact::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: -z
            -z
            sage: -L.one()
            -1
            sage: -(-L.one()) == L.one()
            True

            sage: L([1, 2, 3], constant=2) - L([0, 1], degree=5, constant=2)
            1 + z + 3*z^2 + 2*z^3 + 2*z^4

            sage: -L(0)
            0
        """
        P = self.parent()
        coeff_stream = self._coeff_stream
        if isinstance(coeff_stream, Stream_zero):
            return self
        if isinstance(coeff_stream, Stream_exact):
            initial_coefficients = [-v for v in coeff_stream._initial_coefficients]
            constant = -coeff_stream._constant
            coeff_stream = Stream_exact(initial_coefficients, P._sparse,
                                        constant=constant,
                                        degree=coeff_stream._degree,
                                        order=coeff_stream.order())
            return P.element_class(P, coeff_stream)
        # -(-f) = f
        if isinstance(coeff_stream, Stream_neg):
            return P.element_class(P, coeff_stream._series)
        return P.element_class(P, Stream_neg(coeff_stream))


class LazyCauchyProductSeries(LazyModuleElement):
    r"""
    A class for series where multiplication is the Cauchy product.

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: f = 1 / (1 - z)
        sage: f
        1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
        sage: f * (1 - z)
        1 + O(z^7)

        sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
        sage: f = 1 / (1 - z)
        sage: f
        1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
    """
    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
            sage: M = L(lambda n: n, valuation=0)
            sage: M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: N = M * (1 - M)
            sage: N
            z + z^2 - z^3 - 6*z^4 - 15*z^5 - 29*z^6 + O(z^7)

            sage: p = (1 - z)*(1 + z^2)^3 * z^-2
            sage: p
            z^-2 - z^-1 + 3 - 3*z + 3*z^2 - 3*z^3 + z^4 - z^5
            sage: M = L(lambda n: n, valuation=-2, degree=5, constant=2)
            sage: M
            -2*z^-2 - z^-1 + z + 2*z^2 + 3*z^3 + 4*z^4 + 2*z^5 + 2*z^6 + 2*z^7 + O(z^8)
            sage: M * p
            -2*z^-4 + z^-3 - 5*z^-2 + 4*z^-1 - 2 + 7*z + 5*z^2 + 5*z^3
             + 7*z^4 - 2*z^5 + 4*z^6 - 5*z^7 + z^8 - 2*z^9
            sage: M * p == p * M
            True

            sage: q = (1 - 2*z)*(1 + z^2)^3 * z^-2
            sage: q * M
            -2*z^-4 + 3*z^-3 - 4*z^-2 + 10*z^-1 + 11*z + 2*z^2 - 3*z^3
             - 6*z^4 - 22*z^5 - 14*z^6 - 27*z^7 - 16*z^8 - 20*z^9
             - 16*z^10 - 16*z^11 - 16*z^12 + O(z^13)
            sage: q * M == M * q
            True

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: N = L(lambda n: 1, valuation=0); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: M * N
            z + 3*z^2 + 6*z^3 + 10*z^4 + 15*z^5 + 21*z^6 + O(z^7)

            sage: L.one() * M is M
            True
            sage: M * L.one() is M
            True

        Multiplication of series with eventually constant
        coefficients may yield another such series::

            sage: L.<z> = LazyLaurentSeriesRing(SR)
            sage: var("a b c d e u v w")
            (a, b, c, d, e, u, v, w)
            sage: s = a/z^2 + b*z + c*z^2 + d*z^3 + e*z^4
            sage: t = L([u, v], constant=w, valuation=-1)
            sage: s1 = s.approximate_series(44)
            sage: t1 = t.approximate_series(44)
            sage: s1 * t1 - (s * t).approximate_series(42)
            O(z^42)
        """
        P = self.parent()
        left = self._coeff_stream
        right = other._coeff_stream

        # Check some trivial products
        if isinstance(left, Stream_zero) or isinstance(right, Stream_zero):
            return P.zero()
        if isinstance(left, Stream_exact) and left._initial_coefficients == (P._coeff_ring.one(),) and left.order() == 0:
            return other  # self == 1
        if isinstance(right, Stream_exact) and right._initial_coefficients == (P._coeff_ring.one(),) and right.order() == 0:
            return self  # right == 1

        # The product is exact if and only if both factors are exact
        # and one of the factors has eventually 0 coefficients:
        # (p + a x^d/(1-x))(q + b x^e/(1-x))
        # = p q + (a x^d q + b x^e p)/(1-x) + a b x^(d+e)/(1-x)^2
        if (isinstance(left, Stream_exact)
            and isinstance(right, Stream_exact)
            and not (left._constant and right._constant)):
            il = left._initial_coefficients
            ir = right._initial_coefficients
            initial_coefficients = [sum(il[k]*ir[n-k]
                                        for k in range(max(n - len(ir) + 1, 0),
                                                       min(len(il) - 1, n) + 1))
                                    for n in range(len(il) + len(ir) - 1)]
            lv = left.order()
            rv = right.order()
            # The coefficients of the series (a * x^d * q)/(1-x) are
            # eventually equal to `a * q(1)`, and its initial
            # coefficients are the cumulative sums of the
            # coefficients of q.
            if right._constant:
                d = right._degree
                c = left._constant # this is zero
                # left._constant must be 0 and thus len(il) >= 1
                for k in range(len(il)-1):
                    c += il[k] * right._constant
                    initial_coefficients[d - rv + k] += c
                c += il[-1] * right._constant
            elif left._constant:
                d = left._degree
                c = right._constant # this is zero
                # left._constant must be 0 and thus len(il) >= 1
                for k in range(len(ir)-1):
                    c += left._constant * ir[k]
                    initial_coefficients[d - lv + k] += c
                c += left._constant * ir[-1]
            else:
                c = left._constant # this is zero
            coeff_stream = Stream_exact(initial_coefficients, P._sparse,
                                                   order=lv+rv, constant=c)
            return P.element_class(P, coeff_stream)

        return P.element_class(P, Stream_cauchy_mul(left, right))

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES:

        Lazy Laurent series that have a dense implementation can be inverted::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: ~(1 - z)
            1 + z + z^2 + O(z^3)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: P = ~M; P
            z^-1 - 2 + z + O(z^6)

        Lazy Laurent series that have a sparse implementation can be inverted::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: P = ~M; P
            z^-1 - 2 + z + O(z^6)

            sage: ~(~(1 - z))
            1 - z

        Lazy Laurent series that are known to be exact can be inverted::

            sage: ~z
            z^-1

        TESTS::

            sage: L.<x> = LazyLaurentSeriesRing(QQ)
            sage: g = L([2], valuation=-1, constant=1); g
            2*x^-1 + 1 + x + x^2 + O(x^3)
            sage: g*g^-1
            1 + O(x^7)

        """
        P = self.parent()
        coeff_stream = self._coeff_stream
        # the inverse is exact if and only if coeff_stream corresponds to one of
        # cx^d/(1-x) ... (c, ...)
        # cx^d       ... (c, 0, ...)
        # cx^d (1-x) ... (c, -c, 0, ...)
        if isinstance(coeff_stream, Stream_exact):
            initial_coefficients = coeff_stream._initial_coefficients
            if not initial_coefficients:
                i = ~coeff_stream._constant
                v = -coeff_stream.order()
                c = P._coeff_ring.zero()
                coeff_stream = Stream_exact((i, -i), P._sparse,
                                                       order=v, constant=c)
                return P.element_class(P, coeff_stream)
            if len(initial_coefficients) == 1 and not coeff_stream._constant:
                i = ~initial_coefficients[0]
                v = -coeff_stream.order()
                c = P._coeff_ring.zero()
                coeff_stream = Stream_exact((i,), P._sparse,
                                                       order=v, constant=c)
                return P.element_class(P, coeff_stream)
            if (len(initial_coefficients) == 2
                and not (initial_coefficients[0] + initial_coefficients[1])
                and  not coeff_stream._constant):
                v = -coeff_stream.order()
                c = ~initial_coefficients[0]
                coeff_stream = Stream_exact((), P._sparse,
                                                       order=v, constant=c)
                return P.element_class(P, coeff_stream)

        # (f^-1)^-1 = f
        if isinstance(coeff_stream, Stream_cauchy_invert):
            return P.element_class(P, coeff_stream._series)
        return P.element_class(P, Stream_cauchy_invert(coeff_stream))


class LazyLaurentSeries(LazyCauchyProductSeries):
    r"""
    A Laurent series where the coefficients are computed lazily.

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)

    We can build a series from a function and specify if the series
    eventually takes a constant value::

        sage: f = L(lambda i: i, valuation=-3, constant=-1, degree=3)
        sage: f
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + O(z^6)
        sage: f[-2]
        -2
        sage: f[10]
        -1
        sage: f[-5]
        0

        sage: f = L(lambda i: i, valuation=-3)
        sage: f
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 + 3*z^3 + O(z^4)
        sage: f[20]
        20

    Anything that converts into a polynomial can be input, where
    we can also specify the valuation or if the series eventually
    takes a constant value::

        sage: L([-5,2,0,5])
        -5 + 2*z + 5*z^3
        sage: L([-5,2,0,5], constant=6)
        -5 + 2*z + 5*z^3 + 6*z^4 + 6*z^5 + 6*z^6 + O(z^7)
        sage: L([-5,2,0,5], degree=6, constant=6)
        -5 + 2*z + 5*z^3 + 6*z^6 + 6*z^7 + 6*z^8 + O(z^9)
        sage: L([-5,2,0,5], valuation=-2, degree=3, constant=6)
        -5*z^-2 + 2*z^-1 + 5*z + 6*z^3 + 6*z^4 + 6*z^5 + O(z^6)
        sage: L([-5,2,0,5], valuation=5)
        -5*z^5 + 2*z^6 + 5*z^8
        sage: L({-2:9, 3:4}, constant=2, degree=5)
        9*z^-2 + 4*z^3 + 2*z^5 + 2*z^6 + 2*z^7 + O(z^8)

    We can also perform arithmetic::

        sage: f = 1 / (1 - z - z^2)
        sage: f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + O(z^7)
        sage: f.coefficient(100)
        573147844013817084101
        sage: f = (z^-2 - 1 + 2*z) / (z^-1 - z + 3*z^2)
        sage: f
        z^-1 - z^2 - z^4 + 3*z^5 + O(z^6)

    However, we may not always be able to know when a result is
    exactly a polynomial::

        sage: f * (z^-1 - z + 3*z^2)
        z^-2 - 1 + 2*z + O(z^5)

    TESTS::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
        sage: f = 1 / (1 - z - z^2)
        sage: TestSuite(f).run()

        sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
        sage: f = 1 / (1 - z - z^2)
        sage: TestSuite(f).run()
    """

    def change_ring(self, ring):
        r"""
        Return ``self`` with coefficients converted to elements of ``ring``.

        INPUT:

        - ``ring`` -- a ring

        EXAMPLES:

        Dense Implementation::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: s = 2 + z
            sage: t = s.change_ring(QQ)
            sage: t^-1
            1/2 - 1/4*z + 1/8*z^2 - 1/16*z^3 + 1/32*z^4 - 1/64*z^5 + 1/128*z^6 + O(z^7)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: N = M.change_ring(QQ)
            sage: N.parent()
            Lazy Laurent Series Ring in z over Rational Field
            sage: M.parent()
            Lazy Laurent Series Ring in z over Integer Ring

        Sparse Implementation::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: M.parent()
            Lazy Laurent Series Ring in z over Integer Ring
            sage: N = M.change_ring(QQ)
            sage: N.parent()
            Lazy Laurent Series Ring in z over Rational Field
            sage: M ^-1
            z^-1 - 2 + z + O(z^6)
        """
        from .lazy_laurent_series_ring import LazyLaurentSeriesRing
        Q = LazyLaurentSeriesRing(ring, names=self.parent().variable_names())
        return Q.element_class(Q, self._coeff_stream)

    def __call__(self, g):
        r"""
        Return the composition of ``self`` with ``g``.

        Given two Laurent Series `f` and `g` over the same base ring, the
        composition `(f \circ g)(z) = f(g(z))` is defined if and only if:

        - `g = 0` and `val(f) >= 0`,
        - `g` is non-zero and `f` has only finitely many non-zero coefficients,
        - `g` is non-zero and `val(g) > 0`.

        INPUT:

        - ``g`` -- other series

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = z^2 + 1 + z
            sage: f(0)
            1
            sage: f(L(0))
            1
            sage: f(f)
            3 + 3*z + 4*z^2 + 2*z^3 + z^4
            sage: g = z^-3/(1-2*z); g
            z^-3 + 2*z^-2 + 4*z^-1 + 8 + 16*z + 32*z^2 + 64*z^3 + O(z^4)
            sage: f(g)
            z^-6 + 4*z^-5 + 12*z^-4 + 33*z^-3 + 82*z^-2 + 196*z^-1 + 457 + O(z)
            sage: g^2 + 1 + g
            z^-6 + 4*z^-5 + 12*z^-4 + 33*z^-3 + 82*z^-2 + 196*z^-1 + 457 + O(z)

            sage: f = z^-2 + z + 4*z^3
            sage: f(f)
            4*z^-6 + 12*z^-3 + z^-2 + 48*z^-1 + 12 + O(z)
            sage: f^-2 + f + 4*f^3
            4*z^-6 + 12*z^-3 + z^-2 + 48*z^-1 + 12 + O(z)
            sage: f(g)
            4*z^-9 + 24*z^-8 + 96*z^-7 + 320*z^-6 + 960*z^-5 + 2688*z^-4 + 7169*z^-3 + O(z^-2)
            sage: g^-2 + g + 4*g^3
            4*z^-9 + 24*z^-8 + 96*z^-7 + 320*z^-6 + 960*z^-5 + 2688*z^-4 + 7169*z^-3 + O(z^-2)

            sage: f = z^-3 + z^-2 + 1 / (1 + z^2); f
            z^-3 + z^-2 + 1 - z^2 + O(z^4)
            sage: g = z^3 / (1 + z - z^3); g
            z^3 - z^4 + z^5 - z^7 + 2*z^8 - 2*z^9 + O(z^10)
            sage: f(g)
            z^-9 + 3*z^-8 + 3*z^-7 - z^-6 - 4*z^-5 - 2*z^-4 + z^-3 + O(z^-2)
            sage: g^-3 + g^-2 + 1 / (1 + g^2)
            z^-9 + 3*z^-8 + 3*z^-7 - z^-6 - 4*z^-5 - 2*z^-4 + z^-3 + O(z^-2)

            sage: f = z^-3
            sage: g = z^-2 + z^-1
            sage: g^(-3)
            z^6 - 3*z^7 + 6*z^8 - 10*z^9 + 15*z^10 - 21*z^11 + 28*z^12 + O(z^13)
            sage: f(g)
            z^6 - 3*z^7 + 6*z^8 - 10*z^9 + 15*z^10 - 21*z^11 + 28*z^12 + O(z^13)

            sage: f = z^2 + z^3
            sage: g = z^-3 + z^-2
            sage: f^-3 + f^-2
            z^-6 - 3*z^-5 + 7*z^-4 - 12*z^-3 + 18*z^-2 - 25*z^-1 + 33 + O(z)
            sage: g(f)
            z^-6 - 3*z^-5 + 7*z^-4 - 12*z^-3 + 18*z^-2 - 25*z^-1 + 33 + O(z)
            sage: g^2 + g^3
            z^-9 + 3*z^-8 + 3*z^-7 + 2*z^-6 + 2*z^-5 + z^-4
            sage: f(g)
            z^-9 + 3*z^-8 + 3*z^-7 + 2*z^-6 + 2*z^-5 + z^-4

            sage: f = L(lambda n: n, valuation=0); f
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: f(z^2)
            z^2 + 2*z^4 + 3*z^6 + O(z^7)

            sage: f = L(lambda n: n, valuation=-2); f
            -2*z^-2 - z^-1 + z + 2*z^2 + 3*z^3 + 4*z^4 + O(z^5)
            sage: f3 = f(z^3); f3
            -2*z^-6 - z^-3 + O(z)
            sage: [f3[i] for i in range(-6,13)]
            [-2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4]

        We compose a Laurent polynomial with a generic element::

            sage: R.<x> = QQ[]
            sage: f = z^2 + 1 + z^-1
            sage: g = x^2 + x + 3
            sage: f(g)
            (x^6 + 3*x^5 + 12*x^4 + 19*x^3 + 37*x^2 + 28*x + 31)/(x^2 + x + 3)
            sage: f(g) == g^2 + 1 + g^-1
            True

        We compose with another lazy Laurent series::

            sage: LS.<y> = LazyLaurentSeriesRing(QQ)
            sage: f = z^2 + 1 + z^-1
            sage: fy = f(y); fy
            y^-1 + 1 + y^2
            sage: fy.parent() is LS
            True
            sage: g = y - y
            sage: f(g)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: the valuation of the series must be nonnegative

            sage: g = 1 - y
            sage: f(g)
            3 - y + 2*y^2 + y^3 + y^4 + y^5 + O(y^6)
            sage: g^2 + 1 + g^-1
            3 - y + 2*y^2 + y^3 + y^4 + y^5 + O(y^6)

            sage: f = L(lambda n: n, valuation=0); f
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: f(0)
            0
            sage: f(y)
            y + 2*y^2 + 3*y^3 + 4*y^4 + 5*y^5 + 6*y^6 + O(y^7)
            sage: fp = f(y - y)
            sage: fp == 0
            True
            sage: fp.parent() is LS
            True

            sage: f = z^2 + 3 + z
            sage: f(y - y)
            3

        With both of them sparse::

            sage: L.<z> = LazyLaurentSeriesRing(QQ, sparse=True)
            sage: LS.<y> = LazyLaurentSeriesRing(QQ, sparse=True)
            sage: f = L(lambda n: 1, valuation=0); f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: f(y^2)
            1 + y^2 + y^4 + y^6 + O(y^7)

            sage: fp = f - 1 + z^-2; fp
            z^-2 + z + z^2 + z^3 + z^4 + O(z^5)
            sage: fpy = fp(y^2); fpy
            y^-4 + y^2 + O(y^3)
            sage: fpy.parent() is LS
            True
            sage: [fpy[i] for i in range(-4,11)]
            [1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]

            sage: g = LS(valuation=2, constant=1); g
            y^2 + y^3 + y^4 + O(y^5)
            sage: fg = f(g); fg
            1 + y^2 + y^3 + 2*y^4 + 3*y^5 + 5*y^6 + O(y^7)
            sage: 1 + g + g^2 + g^3 + g^4 + g^5 + g^6
            1 + y^2 + y^3 + 2*y^4 + 3*y^5 + 5*y^6 + O(y^7)

            sage: h = LS(lambda n: 1 if n % 2 else 0, valuation=2); h
            y^3 + y^5 + y^7 + O(y^9)
            sage: fgh = fg(h); fgh
            1 + y^6 + O(y^7)
            sage: [fgh[i] for i in range(0, 15)]
            [1, 0, 0, 0, 0, 0, 1, 0, 2, 1, 3, 3, 6, 6, 13]
            sage: t = 1 + h^2 + h^3 + 2*h^4 + 3*h^5 + 5*h^6
            sage: [t[i] for i in range(0, 15)]
            [1, 0, 0, 0, 0, 0, 1, 0, 2, 1, 3, 3, 6, 6, 13]

        We look at mixing the sparse and the dense::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = L(lambda n: 1, valuation=0); f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: g = LS(lambda n: 1, valuation=1); g
            y + y^2 + y^3 + y^4 + y^5 + y^6 + y^7 + O(y^8)
            sage: f(g)
            1 + y + 2*y^2 + 4*y^3 + 8*y^4 + 16*y^5 + 32*y^6 + O(y^7)

            sage: f = z^-2 + 1 + z
            sage: g = 1/(y*(1-y)); g
            y^-1 + 1 + y + y^2 + y^3 + y^4 + y^5 + O(y^6)
            sage: f(g)
            y^-1 + 2 + y + 2*y^2 - y^3 + 2*y^4 + y^5 + O(y^6)
            sage: g^-2 + 1 + g
            y^-1 + 2 + y + 2*y^2 - y^3 + 2*y^4 + y^5 + O(y^6)

            sage: f = z^-3 + z^-2 + 1
            sage: g = 1/(y^2*(1-y)); g
            y^-2 + y^-1 + 1 + y + y^2 + y^3 + y^4 + O(y^5)
            sage: f(g)
            1 + y^4 - 2*y^5 + 2*y^6 + O(y^7)
            sage: g^-3 + g^-2 + 1
            1 + y^4 - 2*y^5 + 2*y^6 + O(y^7)
            sage: z(y)
            y

        We look at cases where the composition does not exist.
        `g = 0` and `val(f) < 0`::

            sage: g = L(0)
            sage: f = z^-1 + z^-2
            sage: f.valuation() < 0
            True
            sage: f(g)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: the valuation of the series must be nonnegative

        `g \neq 0` and `val(g) \leq 0` and `f` has infinitely many
        non-zero coefficients`::

            sage: g = z^-1 + z^-2
            sage: g.valuation() <= 0
            True
            sage: f = L(lambda n: n, valuation=0)
            sage: f(g)
            Traceback (most recent call last):
            ...
            ValueError: can only compose with a positive valuation series

            sage: f = L(lambda n: n, valuation=1)
            sage: f(1 + z)
            Traceback (most recent call last):
            ...
            ValueError: can only compose with a positive valuation series
        """
        # f = self and compute f(g)
        P = g.parent()

        # f = 0
        if isinstance(self._coeff_stream, Stream_zero):
            return self

        # g = 0 case
        if ((not isinstance(g, LazyModuleElement) and not g)
            or (isinstance(g, LazyModuleElement)
                and isinstance(g._coeff_stream, Stream_zero))):
            if self._coeff_stream._approximate_order >= 0:
                return P(self[0])
            # Perhaps we just don't yet know if the valuation is non-negative
            if any(self._coeff_stream[i] for i in range(self._coeff_stream._approximate_order, 0)):
                raise ZeroDivisionError("the valuation of the series must be nonnegative")
            self._coeff_stream._approximate_order = 0
            return P(self[0])

        # f has finite length and f != 0
        if isinstance(self._coeff_stream, Stream_exact) and not self._coeff_stream._constant:
            # constant polynomial
            R = self.parent()._laurent_poly_ring
            poly = self._coeff_stream._polynomial_part(R)
            if poly.is_constant():
                return self
            if not isinstance(g, LazyModuleElement):
                return poly(g)
            # g also has finite length, compose the polynomials
            if isinstance(g._coeff_stream, Stream_exact) and not g._coeff_stream._constant:
                R = P._laurent_poly_ring
                g_poly = g._coeff_stream._polynomial_part(R)
                try:
                    ret = poly(g_poly)
                except (ValueError, TypeError):  # the result is not a Laurent polynomial
                    ret = None
                if ret is not None and ret.parent() is R:
                    val = ret.valuation()
                    deg = ret.degree() + 1
                    initial_coefficients = [ret[i] for i in range(val, deg)]
                    coeff_stream = Stream_exact(initial_coefficients,
                             self._coeff_stream._is_sparse, constant=P.base_ring().zero(),
                             degree=deg, order=val)
                    return P.element_class(P, coeff_stream)

            # Return the sum since g is not known to be finite or we do not get a Laurent polynomial
            # TODO: Optimize when f has positive valuation
            ret = P.zero()
            # We build this iteratively so each power can benefit from the caching
            # Equivalent to P.sum(poly[i] * g**i for i in range(poly.valuation(), poly.degree()+1))
            # We could just do "return poly(g)" if we don't care about speed
            d = poly.degree()
            v = poly.valuation()
            if d >= 0:
                ind = max(0, v)
                gp = P.one() if ind == 0 else g ** ind
                for i in range(ind, d):
                    if poly[i]:
                        ret += poly[i] * gp
                    gp *= g
                ret += poly[d] * gp
            if v < 0:
                gi = ~g
                ind = min(d, -1)
                gp = gi if ind == -1 else gi ** -ind
                for i in range(ind, v, -1):
                    if poly[i]:
                        ret += poly[i] * gp
                    gp *= gi
                ret += poly[v] * gp
            return ret

        # g != 0 and val(g) > 0
        if not isinstance(g, LazyModuleElement):
            # TODO: Implement case for a regular (Laurent)PowerSeries element
            #   as we can use the (default?) order given
            try:
                g = self.parent()(g)
            except (TypeError, ValueError):
                raise NotImplementedError("can only compose with a lazy series")
        # Perhaps we just don't yet know if the valuation is positive
        if g._coeff_stream._approximate_order <= 0:
            if any(g._coeff_stream[i] for i in range(g._coeff_stream._approximate_order, 1)):
                raise ValueError("can only compose with a positive valuation series")
            g._coeff_stream._approximate_order = 1

        coeff_stream = Stream_cauchy_compose(self._coeff_stream, g._coeff_stream)
        return P.element_class(P, coeff_stream)

    compose = __call__

    def _div_(self, other):
        r"""
        Return ``self`` divided by ``other``.

        INPUT:

        - ``other`` -- nonzero series

        EXAMPLES:

        Lazy Laurent series that have a dense implementation can be divided::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: z/(1 - z)
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + O(z^8)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: N = L(lambda n: 1, valuation=0); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: P = M / N; P
            z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)

        Lazy Laurent series that have a sparse implementation can be divided::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: N = L(lambda n: 1, valuation=0); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: P = M / N; P
            z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)

        Lazy Laurent series that are known to be exact can be divided::

            M = z^2 + 2*z + 1
            N = z + 1
            O = M / N; O
            z + 1

        An example over the ring of symmetric functions::

            sage: e = SymmetricFunctions(QQ).e()
            sage: R.<z> = LazyLaurentSeriesRing(e)
            sage: 1 / (1 - e[1]*z)
            e[] + e[1]*z + e[1, 1]*z^2 + e[1, 1, 1]*z^3 + e[1, 1, 1, 1]*z^4
             + e[1, 1, 1, 1, 1]*z^5 + e[1, 1, 1, 1, 1, 1]*z^6 + O(e[]*z^7)
        """
        if isinstance(other._coeff_stream, Stream_zero):
            raise ZeroDivisionError("cannot divide by 0")

        P = self.parent()
        left = self._coeff_stream
        if isinstance(left, Stream_zero):
            return P.zero()
        right = other._coeff_stream
        if (isinstance(left, Stream_exact)
                and isinstance(right, Stream_exact)):
            if not left._constant and not right._constant:
                R = P._laurent_poly_ring
                pl = left._polynomial_part(R)
                pr = right._polynomial_part(R)
                try:
                    ret = pl / pr
                    ret = P._laurent_poly_ring(ret)
                    initial_coefficients = [ret[i] for i in range(ret.valuation(), ret.degree() + 1)]
                    return P.element_class(P, Stream_exact(initial_coefficients, P._sparse,
                             order=ret.valuation(), constant=left._constant))
                except (TypeError, ValueError, NotImplementedError):
                    # We cannot divide the polynomials, so the result must be a series
                    pass

        return P.element_class(P, Stream_cauchy_mul(left, Stream_cauchy_invert(right)))

    def __pow__(self, n):
        """
        Return the ``n``-th power of the series.

        INPUT:

        - ``n`` -- integer; the power to which to raise the series

        EXAMPLES:

        Lazy Laurent series that have a dense implementation can be
        raised to the power ``n``::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
            sage: (1 - z)^-1
            1 + z + z^2 + O(z^3)
            sage: (1 - z)^0
            1
            sage: (1 - z)^3
            1 - 3*z + 3*z^2 - z^3
            sage: (1 - z)^-3
            1 + 3*z + 6*z^2 + 10*z^3 + 15*z^4 + 21*z^5 + 28*z^6 + O(z^7)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: M^2
            z^2 + 4*z^3 + 10*z^4 + 20*z^5 + 35*z^6 + O(z^7)

        We can create a really large power of a monomial, even with
        the dense implementation::

            sage: z^1000000
            z^1000000

        Lazy Laurent series that have a sparse implementation can be
        raised to the power ``n``::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: M = L(lambda n: n, valuation=0); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
            sage: M^2
            z^2 + 4*z^3 + 10*z^4 + 20*z^5 + 35*z^6 + O(z^7)

        Lazy Laurent series that are known to be exact can be raised
        to the power ``n``::

            sage: z^2
            z^2
            sage: (1 - z)^2
            1 - 2*z + z^2
            sage: (1 + z)^2
            1 + 2*z + z^2

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L(0)^10
            0

        """
        if n == 0:
            return self.parent().one()

        cs = self._coeff_stream
        if (isinstance(cs, Stream_exact)
            and not cs._constant and n in ZZ
            and (n > 0 or len(cs._initial_coefficients) == 1)):
            P = self.parent()
            ret = cs._polynomial_part(P._laurent_poly_ring) ** ZZ(n)
            val = ret.valuation()
            deg = ret.degree() + 1
            initial_coefficients = [ret[i] for i in range(val, deg)]
            return P.element_class(P, Stream_exact(initial_coefficients, P._sparse,
                            constant=cs._constant, degree=deg, order=val))

        return generic_power(self, n)

    def approximate_series(self, prec, name=None):
        """
        Return the Laurent series with absolute precision ``prec`` approximated
        from this series.

        INPUT:

        - ``prec`` -- an integer
        - ``name`` -- name of the variable; if it is ``None``, the name of
          the variable of the series is used

        OUTPUT: a Laurent series with absolute precision ``prec``

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: f = (z - 2*z^3)^5/(1 - 2*z)
            sage: f
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + 32*z^10 - 16*z^11 + O(z^12)
            sage: g = f.approximate_series(10)
            sage: g
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + O(z^10)
            sage: g.parent()
            Power Series Ring in z over Integer Ring
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
            return R([self[i] for i in range(n, prec)], n).add_bigoh(prec)
        else:
            from sage.rings.all import PowerSeriesRing
            R = PowerSeriesRing(S.base_ring(), name=name)
            return R([self[i] for i in range(prec)]).add_bigoh(prec)

    def polynomial(self, degree=None, name=None):
        r"""
        Return ``self`` as a Laurent polynomial if ``self`` is actually so.

        INPUT:

        - ``degree`` -- ``None`` or an integer
        - ``name`` -- name of the variable; if it is ``None``, the name of
          the variable of the series is used

        OUTPUT:

        A Laurent polynomial if the valuation of the series is negative or
        a polynomial otherwise.

        If ``degree`` is not ``None``, the terms of the series of degree
        greater than ``degree`` are truncated first. If ``degree`` is ``None``
        and the series is not a polynomial or a Laurent polynomial, a
        ``ValueError`` is raised.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = L([1,0,0,2,0,0,0,3], valuation=5); f
            z^5 + 2*z^8 + 3*z^12
            sage: f.polynomial()
            3*z^12 + 2*z^8 + z^5

        TESTS::

            sage: g = L([1,0,0,2,0,0,0,3], valuation=-5); g
            z^-5 + 2*z^-2 + 3*z^2
            sage: g.polynomial()
            z^-5 + 2*z^-2 + 3*z^2
            sage: z = L.gen()
            sage: f = (1 + z)/(z^3 - z^5)
            sage: f
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + O(z^4)
            sage: f.polynomial(5)
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + z^4 + z^5
            sage: f.polynomial(0)
            z^-3 + z^-2 + z^-1 + 1
            sage: f.polynomial(-5)
            0
            sage: M = L(lambda n: n^2, valuation=0)
            sage: M.polynomial(3)
            9*z^3 + 4*z^2 + z
            sage: M = L(lambda n: n^2, valuation=0)
            sage: M.polynomial(5)
            25*z^5 + 16*z^4 + 9*z^3 + 4*z^2 + z

            sage: f = 1/(1 + z)
            sage: f.polynomial()
            Traceback (most recent call last):
            ...
            ValueError: not a polynomial

            sage: L.zero().polynomial()
            0
        """
        S = self.parent()

        if degree is None:
            if isinstance(self._coeff_stream, Stream_zero):
                from sage.rings.all import PolynomialRing
                return PolynomialRing(S.base_ring(), name=name).zero()
            elif isinstance(self._coeff_stream, Stream_exact) and not self._coeff_stream._constant:
                m = self._coeff_stream._degree
            else:
                raise ValueError("not a polynomial")
        else:
            m = degree + 1

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            R = LaurentPolynomialRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self[i] for i in range(n, m)]).shift(n)
        else:
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(S.base_ring(), name=name)
            return R([self[i] for i in range(m)])

    def shift(self, n):
        r"""
        Return ``self`` multiplied by the power `z^n`, where `z` is the
        variable of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1 / (1 + 2*z)
            sage: f
            1 - 2*z + 4*z^2 - 8*z^3 + 16*z^4 - 32*z^5 + 64*z^6 + O(z^7)
            sage: f.shift(3)
            z^3 - 2*z^4 + 4*z^5 - 8*z^6 + 16*z^7 - 32*z^8 + 64*z^9 + O(z^10)
            sage: f << -3  # shorthand
            z^-3 - 2*z^-2 + 4*z^-1 - 8 + 16*z - 32*z^2 + 64*z^3 + O(z^4)
            sage: g = z^-3 + 3 + z^2
            sage: g.shift(5)
            z^2 + 3*z^5 + z^7
            sage: L([2,0,3], valuation=2, degree=7, constant=1) << -2
            2 + 3*z^2 + z^5 + z^6 + z^7 + O(z^8)

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: zero = L.zero()
            sage: zero.shift(10) is zero
            True

            sage: f = 1 / (1 + 2*z + z^2)
            sage: f.shift(5).shift(-5) - f
            0

        """
        if isinstance(self._coeff_stream, Stream_zero):
            return self
        elif isinstance(self._coeff_stream, Stream_shift):
            n += self._coeff_stream._shift
            if n:
                coeff_stream = Stream_shift(self._coeff_stream._series, n)
            else:
                coeff_stream = self._coeff_stream._series
        elif isinstance(self._coeff_stream, Stream_exact):
            init_coeff = self._coeff_stream._initial_coefficients
            degree = self._coeff_stream._degree + n
            valuation = self._coeff_stream._approximate_order + n
            coeff_stream = Stream_exact(init_coeff, self._coeff_stream._is_sparse,
                                                   constant=self._coeff_stream._constant,
                                                   order=valuation, degree=degree)
        else:
            coeff_stream = Stream_shift(self._coeff_stream, n)
        P = self.parent()
        return P.element_class(P, coeff_stream)

    __lshift__ = shift

    def __rshift__(self, n):
        r"""
        Return ``self`` multiplied by the power `z^-n`, where `z` is the
        variable of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 + 2*z); f
            1 - 2*z + 4*z^2 - 8*z^3 + 16*z^4 - 32*z^5 + 64*z^6 + O(z^7)
            sage: f >> 3
            z^-3 - 2*z^-2 + 4*z^-1 - 8 + 16*z - 32*z^2 + 64*z^3 + O(z^4)
            sage: f >> -3
            z^3 - 2*z^4 + 4*z^5 - 8*z^6 + 16*z^7 - 32*z^8 + 64*z^9 + O(z^10)
        """
        return self.shift(-n)

    def _format_series(self, formatter, format_strings=False):
        """
        Return ``self`` formatted by ``formatter``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = 1 / (2 - z^2)
            sage: f._format_series(ascii_art, True)
            1/2 + 1/4*z^2 + 1/8*z^4 + 1/16*z^6 + O(z^7)
        """
        P = self.parent()
        R = P._laurent_poly_ring
        z = R.gen()
        cs = self._coeff_stream
        v = cs._approximate_order
        if format_strings:
            strformat = formatter
        else:
            strformat = lambda x: x

        if isinstance(cs, Stream_exact):
            poly = cs._polynomial_part(R)
            if not cs._constant:
                return formatter(poly)
            m = cs._degree + P.options.constant_length
            poly += sum([cs._constant * z**k for k in range(cs._degree, m)])
            return formatter(poly) + strformat(" + O({})".format(formatter(z**m)))

        # This is an inexact series
        m = v + P.options.display_length

        # Use the polynomial printing
        poly = R([self._coeff_stream[i] for i in range(v, m)]).shift(v)
        if not poly:
            return strformat("O({})".format(formatter(z**m)))
        return formatter(poly) + strformat(" + O({})".format(formatter(z**m)))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: z^-3 + z - 5
            z^-3 - 5 + z
            sage: -1/(1 + 2*z)
            -1 + 2*z - 4*z^2 + 8*z^3 - 16*z^4 + 32*z^5 - 64*z^6 + O(z^7)
            sage: -z^-7/(1 + 2*z)
            -z^-7 + 2*z^-6 - 4*z^-5 + 8*z^-4 - 16*z^-3 + 32*z^-2 - 64*z^-1 + O(1)
            sage: L([1,5,0,3], valuation=-1, degree=5, constant=2)
            z^-1 + 5 + 3*z^2 + 2*z^5 + 2*z^6 + 2*z^7 + O(z^8)
            sage: L(constant=5, valuation=2)
            5*z^2 + 5*z^3 + 5*z^4 + O(z^5)
            sage: L(constant=5, degree=-2)
            5*z^-2 + 5*z^-1 + 5 + O(z)
            sage: L(lambda x: x if x < 0 else 0, valuation=-2)
            -2*z^-2 - z^-1 + O(z^5)
            sage: L(lambda x: x if x < 0 else 0, valuation=2)
            O(z^9)
            sage: L(lambda x: x if x > 0 else 0, valuation=-2)
            z + 2*z^2 + 3*z^3 + 4*z^4 + O(z^5)
            sage: L(lambda x: x if x > 0 else 0, valuation=-10)
            O(z^-3)

            sage: L(None, valuation=0)
            Uninitialized Lazy Laurent Series
            sage: L(0)
            0

            sage: R.<x,y> = QQ[]
            sage: L.<z> = LazyLaurentSeriesRing(R)
            sage: z^-2 / (1 - (x-y)*z) + x^4*z^-3 + (1-y)*z^-4
            (-y + 1)*z^-4 + x^4*z^-3 + z^-2 + (x - y)*z^-1
             + (x^2 - 2*x*y + y^2) + (x^3 - 3*x^2*y + 3*x*y^2 - y^3)*z
             + (x^4 - 4*x^3*y + 6*x^2*y^2 - 4*x*y^3 + y^4)*z^2 + O(z^3)
        """
        if isinstance(self._coeff_stream, Stream_zero):
            return '0'
        if isinstance(self._coeff_stream, Stream_uninitialized) and self._coeff_stream._target is None:
            return 'Uninitialized Lazy Laurent Series'
        return self._format_series(repr)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: latex(z^-3 + z - 5)
            \frac{1}{z^{3}} - 5 + z
            sage: latex(-1/(1 + 2*z))
            -1 + 2z - 4z^{2} + 8z^{3} - 16z^{4} + 32z^{5} - 64z^{6} + O(z^{7})
            sage: latex(-z^-7/(1 + 2*z))
            \frac{-1}{z^{7}} + \frac{2}{z^{6}} + \frac{-4}{z^{5}} + \frac{8}{z^{4}}
             + \frac{-16}{z^{3}} + \frac{32}{z^{2}} + \frac{-64}{z} + O(1)
            sage: latex(L([1,5,0,3], valuation=-1, degree=5, constant=2))
            \frac{1}{z} + 5 + 3z^{2} + 2z^{5} + 2z^{6} + 2z^{7} + O(z^{8})
            sage: latex(L(constant=5, valuation=2))
            5z^{2} + 5z^{3} + 5z^{4} + O(z^{5})
            sage: latex(L(constant=5, degree=-2))
            \frac{5}{z^{2}} + \frac{5}{z} + 5 + O(z)
            sage: latex(L(lambda x: x if x < 0 else 0, valuation=-2))
            \frac{-2}{z^{2}} + \frac{-1}{z} + O(z^{5})
            sage: latex(L(lambda x: x if x < 0 else 0, valuation=2))
            O(z^{9})
            sage: latex(L(lambda x: x if x > 0 else 0, valuation=-2))
            z + 2z^{2} + 3z^{3} + 4z^{4} + O(z^{5})
            sage: latex(L(lambda x: x if x > 0 else 0, valuation=-10))
            O(\frac{1}{z^{3}})

            sage: latex(L(None, valuation=0))
            \text{\texttt{Undef}}
            sage: latex(L(0))
            0

            sage: R.<x,y> = QQ[]
            sage: L.<z> = LazyLaurentSeriesRing(R)
            sage: latex(z^-2 / (1 - (x-y)*z) + x^4*z^-3 + (1-y)*z^-4)
            \frac{-y + 1}{z^{4}} + \frac{x^{4}}{z^{3}} + \frac{1}{z^{2}}
             + \frac{x - y}{z} + x^{2} - 2 x y + y^{2}
             + \left(x^{3} - 3 x^{2} y + 3 x y^{2} - y^{3}\right)z
             + \left(x^{4} - 4 x^{3} y + 6 x^{2} y^{2} - 4 x y^{3} + y^{4}\right)z^{2}
             + O(z^{3})
        """
        from sage.misc.latex import latex
        if isinstance(self._coeff_stream, Stream_zero):
            return latex('0')
        if isinstance(self._coeff_stream, Stream_uninitialized) and self._coeff_stream._target is None:
            return latex("Undef")
        return self._format_series(latex)

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: L.<z> = LazyLaurentSeriesRing(e)
            sage: L.options.display_length = 3
            sage: ascii_art(1 / (1 - e[1]*z))
            e[] + e[1]*z + e[1, 1]*z^2 + O(e[]*z^3)
            sage: L.options._reset()
        """
        from sage.typeset.ascii_art import ascii_art, AsciiArt
        if isinstance(self._coeff_stream, Stream_zero):
            return AsciiArt('0')
        if isinstance(self._coeff_stream, Stream_uninitialized) and self._coeff_stream._target is None:
            return AsciiArt('Uninitialized Lazy Laurent Series')
        return self._format_series(ascii_art, True)

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: L.<z> = LazyLaurentSeriesRing(e)
            sage: L.options.display_length = 3
            sage: unicode_art(1 / (1 - e[1]*z))
            e[] + e[1]*z + e[1, 1]*z^2 + O(e[]*z^3)
            sage: L.options._reset()
        """
        from sage.typeset.unicode_art import unicode_art, UnicodeArt
        if isinstance(self._coeff_stream, Stream_zero):
            return UnicodeArt('0')
        if isinstance(self._coeff_stream, Stream_uninitialized) and self._coeff_stream._target is None:
            return UnicodeArt('Uninitialized Lazy Laurent Series')
        return self._format_series(unicode_art, True)
