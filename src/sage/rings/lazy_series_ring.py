r"""
Lazy Series Rings

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import parent

from sage.categories.algebras import Algebras
from sage.categories.rings import Rings
from sage.categories.integral_domains import IntegralDomains
from sage.categories.fields import Fields
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.lazy_series import (LazyModuleElement,
                                    LazyLaurentSeries,
                                    LazyDirichletSeries)
from sage.structure.global_options import GlobalOptions
from sage.symbolic.ring import SR

from sage.data_structures.stream import (
    Stream_zero,
    Stream_function,
    Stream_exact,
    Stream_uninitialized
)

class LazySeriesRing(UniqueRepresentation, Parent):
    """
    Abstract base class for lazy series.
    """
    def _element_constructor_(self, x=None, valuation=None, degree=None, constant=None, coefficients=None):
        r"""
        Construct a lazy series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a series
        - ``valuation`` -- integer (optional); integer; a lower bound for
          the valuation of the series
        - ``degree`` -- (optional) the degree when the series is ``constant``
        - ``constant`` -- (optional) the eventual constant of the series
        - ``coefficients`` -- (optional) a callable that defines the
          coefficients of the series; must be ``None`` if ``x`` is provided;
          see note below

        If ``valuation`` is specified and ``x`` is convertible into
        an element of the underlying ring corresponding to series
        with finite support or ``x`` is a lazy series of the same
        parent, then the data is shifted so that the result has the
        specified valuation.

        .. WARNING::

            If ``valuation`` is specified and ``x`` is a lazy series, then
            the valuation will be computed. If the series ``x`` is not
            known to be zero, then this will run forever.

        .. NOTE::

            When working over a base ring that takes callables as valid
            input, then passing a function as ``x`` might be converted to
            the base ring. If instead the input is to be treated as the
            function giving the coefficients of the lazy series being
            constructed, then use the ``coefficients`` argument in this
            case and do not provide ``x``.

        .. WARNING::

            Polynomials, but also :class:`LazyLaurentSeries` and
            :class:`LazyDirichletSeries` are callable.  Therefore, an
            argument ``x`` which is not convertible into an element
            of the underlying ring corresponding to series with
            finite support is interpreted as a function providing the
            coefficients when evaluated at integers.  Examples are
            provided below.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)

            sage: L(lambda i: i, valuation=5, constant=1, degree=10)
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)
            sage: L(lambda i: i, valuation=5, constant=(1, 10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: X = L(constant=5, degree=2); X
            5*z^2 + 5*z^3 + 5*z^4 + O(z^5)
            sage: X.valuation()
            2

            sage: def g(i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return 1 + sum(k for k in range(i+1))
            sage: e = L(g, valuation=-5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + O(z^2)
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + O(z^12)
            sage: f.coefficient(10)
            0
            sage: f[20]
            9
            sage: f[30]
            -219

            sage: L(valuation=2, constant=1)
            z^2 + z^3 + z^4 + O(z^5)
            sage: L(constant=1)
            Traceback (most recent call last):
            ...
            ValueError: you must specify the degree for the polynomial 0

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], valuation=-5)
            sage: f
            z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2
            sage: g = L([1,3,5,7,9], valuation=5, constant=-1)
            sage: g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + O(z^13)

        Finally, ``x`` can be a Laurent polynomial::

            sage: P.<x> = LaurentPolynomialRing(QQ)
            sage: p = x^-2 + 3*x^3
            sage: L.<x> = LazyLaurentSeriesRing(ZZ)
            sage: L(p)
            x^-2 + 3*x^3

            sage: L(p, valuation=0)
            1 + 3*x^5

            sage: L(p, valuation=1)
            x + 3*x^6

        We construct a lazy Laurent series over another lazy Laurent series::

            sage: R.<s> = LazyLaurentSeriesRing(QQ)
            sage: L.<z> = LazyLaurentSeriesRing(R)
            sage: e = L(lambda n: 1/factorial(n), 0); e
            1 + z + 1/2*z^2 + 1/6*z^3 + 1/24*z^4 + 1/120*z^5 + 1/720*z^6 + O(z^7)
            sage: L(lambda n: 1/(1 + s^n), 0)
            1/2 + (1 - s + s^2 - s^3 + s^4 - s^5 + s^6 + O(s^7))*z
             + (1 - s^2 + s^4 - s^6 + O(s^7))*z^2
             + (1 - s^3 + s^6 + O(s^7))*z^3 + (1 - s^4 + O(s^7))*z^4
             + (1 - s^5 + O(s^7))*z^5 + (1 - s^6 + O(s^7))*z^6 + O(z^7)

        We note that ``e`` builds correctly because ``R`` additionally
        requires the valuation to be specified.

        In the next example the argument is interpreted as a constant
        polynomial, which happens to be a Dirichlet series::

            sage: D = LazyDirichletSeriesRing(QQ, "s")
            sage: L.<z> = LazyLaurentSeriesRing(D)
            sage: L(lambda n: 1/factorial(n), valuation=0)
            (1 + 1/2/2^s + 1/6/3^s + 1/24/4^s + 1/120/5^s + 1/720/6^s + 1/5040/7^s + O(1/(8^s)))

        We can also specify that the given function should be
        interpreted as the coefficients of the Laurent series::

            sage: L(coefficients=lambda n: 1/factorial(n), valuation=0)
            1 + z + 1/2*z^2 + 1/6*z^3 + 1/24*z^4 + 1/120*z^5 + 1/720*z^6 + O(z^7)

        When the argument ``x`` is callable and not convertible into
        an element of the underlying ring of series of finite
        support, it is evaluated at integers to compute the
        coefficients::

            sage: R.<q> = QQ[]
            sage: D = LazyDirichletSeriesRing(ZZ, 't')
            sage: D(1+2*q)
            3 + 5/2^t + 7/3^t + 9/4^t + 11/5^t + 13/6^t + 15/7^t + O(1/(8^t))

        In this example, the Dirichlet series ``m`` is considered as an
        element in the base ring::

            sage: m = D(moebius)
            sage: s = L(m, valuation=0)
            sage: s[0]
            1 - 1/(2^s) - 1/(3^s) - 1/(5^s) + 1/(6^s) - 1/(7^s) + O(1/(8^s))
            sage: s[1]
            0

        TESTS:

        Checking the valuation is consistent::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L([0,0,2,3], valuation=-4)
            2*z^-4 + 3*z^-3
            sage: L(range(5), valuation=-4)
            z^-4 + 2*z^-3 + 3*z^-2 + 4*z^-1
            sage: P.<x> = ZZ[]
            sage: L(x^2 + x^5, valuation=-4)
            z^-4 + z^-1
            sage: L(1, valuation=-4)
            z^-4
            sage: L(L(1), valuation=-4)
            z^-4
            sage: L(1/(1-z), valuation=-4)
            z^-4 + z^-3 + z^-2 + z^-1 + 1 + z + z^2 + O(z^3)
            sage: L(z^-3/(1-z), valuation=-4)
            z^-4 + z^-3 + z^-2 + z^-1 + 1 + z + z^2 + O(z^3)
            sage: L(z^3/(1-z), valuation=-4)
            z^-4 + z^-3 + z^-2 + z^-1 + 1 + z + z^2 + O(z^3)

            sage: L(z^3/(1-z), valuation=0)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L(lambda n: 1/(n+1), degree=3)
            Traceback (most recent call last):
            ...
            ValueError: the valuation must be specified

            sage: L(5, valuation=3.1)
            Traceback (most recent call last):
            ...
            ValueError: the valuation must be an integer

            sage: L(5, valuation=6/2)
            5*z^3

        Checking that series are not interpreted as coefficients when
        they can be interpreted as series::

            sage: P.<s> = ZZ[]
            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: M.<w> = LazyLaurentSeriesRing(QQ)
            sage: L(M(s^2 + s^5), valuation=-4)
            z^-4 + z^-1

            sage: D = LazyDirichletSeriesRing(ZZ, "s")
            sage: E = LazyDirichletSeriesRing(QQ, "t")
            sage: D(E([1,2,3]))
            1 + 2/2^s + 3/3^s

        This gives zero::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L(lambda n: 0, degree=3, valuation=0)
            0
            sage: L(L.zero(), degree=3)
            0
            sage: L(L.zero(), degree=3, valuation=2)
            0
            sage: L(L.zero(), degree=3, constant=0)
            0
            sage: L(L.zero(), degree=3, valuation=2, constant=0)
            0

        This does not::

            sage: L(lambda n: 0, degree=3, constant=1, valuation=0)
            z^3 + z^4 + z^5 + O(z^6)
            sage: L(L.zero(), degree=-3, constant=1)
            z^-3 + z^-2 + z^-1 + O(1)
            sage: L(L.zero(), valuation=2, constant=1)
            z^2 + z^3 + z^4 + O(z^5)

        This raises an error::

            sage: L(lambda n: 0, valuation=3, constant=1)
            Traceback (most recent call last):
            ...
            ValueError: constant may only be specified if the degree is specified

        We support the old input format for ``constant``::

            sage: f = L(lambda i: i, valuation=-3, constant=-1, degree=3)
            sage: g = L(lambda i: i, valuation=-3, constant=(-1,3))
            sage: f == g
            True
            sage: g = L(lambda i: i, -3, (-1,3))
            sage: f == g
            True

        .. TODO::

            Add a method to change the sparse/dense implementation.

        """
        if valuation is not None and valuation not in ZZ:
            raise ValueError("the valuation must be an integer")

        if x is None and coefficients is None:
            if valuation is None:
                raise ValueError("the valuation must be specified")
            return self.element_class(self, Stream_uninitialized(self._sparse, valuation))

        if coefficients is not None and (x is not None and (not isinstance(x, int) or x)):
            raise ValueError("coefficients must be None if x is provided")

        BR = self.base_ring()
        if isinstance(constant, (tuple, list)):
            constant, degree = constant
        if isinstance(degree, (tuple, list)):
            constant, degree = degree
        if constant is not None:
            constant = BR(constant)

        if coefficients is None:
            # Try to build stuff using the internal polynomial ring constructor
            R = self._internal_poly_ring
            try:
                x = R(x)
            except (TypeError, ValueError):
                pass

            # If x has been converted to the internal polynomial ring
            if parent(x) is R:
                if not x and not constant:
                    return self.zero()
                if x and valuation is not None:
                    x = x.shift(valuation - x.valuation())
                if degree is None and not x:
                    if valuation is None:
                        raise ValueError("you must specify the degree for the polynomial 0")
                    degree = valuation
                if x == R.zero():
                    coeff_stream = Stream_exact([], self._sparse, order=degree, constant=constant)
                    return self.element_class(self, coeff_stream)
                initial_coefficients = [x[i] for i in range(x.valuation(), x.degree() + 1)]
                coeff_stream = Stream_exact(initial_coefficients, self._sparse,
                                            order=x.valuation(), constant=constant, degree=degree)
                return self.element_class(self, coeff_stream)

            # Handle when it is a lazy series
            if isinstance(x, self.Element):
                if x._coeff_stream._is_sparse is not self._sparse:
                    # TODO: Implement a way to make a self._sparse copy
                    raise NotImplementedError("cannot convert between sparse and dense")

                # If x is known to be 0
                if isinstance(x._coeff_stream, Stream_zero):
                    if not constant:
                        if self is parent(x):
                            return x
                        return self.element_class(self, x._coeff_stream)
                    if degree is None:
                        if valuation is None:
                            raise ValueError("you must specify the degree for the polynomial 0")
                        degree = valuation
                    coeff_stream = Stream_exact([], self._sparse, order=degree,
                                                constant=constant)
                    return self.element_class(self, coeff_stream)

                # Make the result exact
                if degree is not None:
                    # truncate the series and then possibly make constant
                    x_val = x._coeff_stream.order()
                    if not valuation:
                        valuation = x_val
                    initial_coefficients = [x[x_val+i] for i in range(degree-valuation)]
                    if not any(initial_coefficients):
                        if not constant:
                            return self.zero()
                        # We learned some stuff about x; pass it along
                        x._coeff_stream._approximate_order += len(initial_coefficients)
                        initial_coefficients = []
                    coeff_stream = Stream_exact(initial_coefficients, self._sparse,
                                                order=valuation, constant=constant, degree=degree)
                    return self.element_class(self, coeff_stream)

                # We are just possibly shifting the result
                ret = self.element_class(self, x._coeff_stream)
                if valuation is None:
                    return ret
                return ret.shift(valuation - x._coeff_stream.order())

        else:
            x = coefficients

        if callable(x):
            if valuation is None:
                raise ValueError("the valuation must be specified")
            if degree is None:
                if constant is not None:
                    raise ValueError("constant may only be specified if the degree is specified")
                coeff_stream = Stream_function(x, self.base_ring(), self._sparse, valuation)
                return self.element_class(self, coeff_stream)

            # degree is not None
            if constant is None:
                constant = BR.zero()
            p = [BR(x(i)) for i in range(valuation, degree)]
            if not any(p) and not constant:
                return self.zero()
            coeff_stream = Stream_exact(p, self._sparse, order=valuation,
                                        constant=constant, degree=degree)
            return self.element_class(self, coeff_stream)

        raise ValueError(f"unable to convert {x} into {self}")


class LazyLaurentSeriesRing(LazySeriesRing):
    """
    The ring of lazy Laurent series.

    The ring of Laurent series over a ring with the usual arithmetic
    where the coefficients are computed lazily.

    INPUT:

    - ``base_ring`` -- base ring
    - ``names`` -- name of the generator
    - ``sparse`` -- (default: ``True``) whether the implementation of
      the series is sparse or not

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(QQ)
        sage: 1 / (1 - z)
        1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
        sage: 1 / (1 - z) == 1 / (1 - z)
        True
        sage: L in Fields
        True

    Lazy Laurent series ring over a finite field::

        sage: L.<z> = LazyLaurentSeriesRing(GF(3)); L
        Lazy Laurent Series Ring in z over Finite Field of size 3
        sage: e = 1 / (1 + z)
        sage: e.coefficient(100)
        1
        sage: e.coefficient(100).parent()
        Finite Field of size 3

    Series can be defined by specifying a coefficient function
    along with a valuation or a degree where after the series
    is evenutally constant::

        sage: R.<x,y> = QQ[]
        sage: L.<z> = LazyLaurentSeriesRing(R)
        sage: def coeff(n):
        ....:     if n < 0:
        ....:         return -2 + n
        ....:     if n == 0:
        ....:         return 6
        ....:     return x + y^n
        sage: f = L(coeff, valuation=-5)
        sage: f
        -7*z^-5 - 6*z^-4 - 5*z^-3 - 4*z^-2 - 3*z^-1 + 6 + (x + y)*z + O(z^2)
        sage: 1 / (1 - f)
        1/7*z^5 - 6/49*z^6 + 1/343*z^7 + 8/2401*z^8 + 64/16807*z^9
         + 17319/117649*z^10 + (1/49*x + 1/49*y - 180781/823543)*z^11 + O(z^12)
        sage: L(coeff, valuation=-3, degree=3, constant=x)
        -5*z^-3 - 4*z^-2 - 3*z^-1 + 6 + (x + y)*z + (y^2 + x)*z^2
         + x*z^3 + x*z^4 + x*z^5 + O(z^6)

    Similarly, we can specify a polynomial or the initial
    coefficients with anything that converts into the
    corresponding Laurent polynomial ring::

        sage: L([1, x, y, 0, x+y])
        1 + x*z + y*z^2 + (x + y)*z^4
        sage: L([1, x, y, 0, x+y], constant=2)
        1 + x*z + y*z^2 + (x + y)*z^4 + 2*z^5 + 2*z^6 + 2*z^7 + O(z^8)
        sage: L([1, x, y, 0, x+y], degree=7, constant=2)
        1 + x*z + y*z^2 + (x + y)*z^4 + 2*z^7 + 2*z^8 + 2*z^9 + O(z^10)
        sage: L([1, x, y, 0, x+y], valuation=-2)
        z^-2 + x*z^-1 + y + (x + y)*z^2
        sage: L([1, x, y, 0, x+y], valuation=-2, constant=3)
        z^-2 + x*z^-1 + y + (x + y)*z^2 + 3*z^3 + 3*z^4 + 3*z^5 + O(z^6)
        sage: L([1, x, y, 0, x+y], valuation=-2, degree=4, constant=3)
        z^-2 + x*z^-1 + y + (x + y)*z^2 + 3*z^4 + 3*z^5 + 3*z^6 + O(z^7)

    Some additional examples over the integer ring::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: L in Fields
        False
        sage: 1 / (1 - 2*z)^3
        1 + 6*z + 24*z^2 + 80*z^3 + 240*z^4 + 672*z^5 + 1792*z^6 + O(z^7)

        sage: R.<x> = LaurentPolynomialRing(ZZ)
        sage: L(x^-2 + 3 + x)
        z^-2 + 3 + z
        sage: L(x^-2 + 3 + x, valuation=-5, constant=2)
        z^-5 + 3*z^-3 + z^-2 + 2*z^-1 + 2 + 2*z + O(z^2)
        sage: L(x^-2 + 3 + x, valuation=-5, degree=0, constant=2)
        z^-5 + 3*z^-3 + z^-2 + 2 + 2*z + 2*z^2 + O(z^3)

    We can also truncate, shift, and make eventually constant any
    Laurent series::

        sage: f = 1 / (z + z^2)
        sage: f
        z^-1 - 1 + z - z^2 + z^3 - z^4 + z^5 + O(z^6)
        sage: L(f, valuation=2)
        z^2 - z^3 + z^4 - z^5 + z^6 - z^7 + z^8 + O(z^9)
        sage: L(f, degree=3)
        z^-1 - 1 + z - z^2
        sage: L(f, degree=3, constant=2)
        z^-1 - 1 + z - z^2 + 2*z^3 + 2*z^4 + 2*z^5 + O(z^6)
        sage: L(f, valuation=1, degree=4)
        z - z^2 + z^3
        sage: L(f, valuation=1, degree=4, constant=5)
        z - z^2 + z^3 + 5*z^4 + 5*z^5 + 5*z^6 + O(z^7)

    Power series can be defined recursively (see :meth:`define()` for
    more examples)::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: s = L(None, valuation=0)
        sage: s.define(1 + z*s^2)
        sage: s
        1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + O(z^7)

    If we do not explcitly know the exact value of every coefficient,
    then equality checking will depend on the computed coefficients.
    If at a certain point we cannot prove two series are different
    (which involves the coefficients we have computed), then we will
    raise an error::

        sage: f = 1 / (z + z^2); f
        z^-1 - 1 + z - z^2 + z^3 - z^4 + z^5 + O(z^6)
        sage: f2 = f * 2  # currently no coefficients computed
        sage: f3 = f * 3  # currently no coefficients computed
        sage: f2 == f3
        Traceback (most recent call last):
        ...
        ValueError: undecidable
        sage: f2  # computes some of the coefficients of f2
        2*z^-1 - 2 + 2*z - 2*z^2 + 2*z^3 - 2*z^4 + 2*z^5 + O(z^6)
        sage: f3  # computes some of the coefficients of f3
        3*z^-1 - 3 + 3*z - 3*z^2 + 3*z^3 - 3*z^4 + 3*z^5 + O(z^6)
        sage: f2 == f3
        False

    The implementation of the ring can be either be a sparse or a dense one.
    The default is a sparse implementation::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: L.is_sparse()
        True
        sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
        sage: L.is_sparse()
        False
    """
    Element = LazyLaurentSeries

    def __init__(self, base_ring, names, sparse=True, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 't')
            sage: elts = L.some_elements()[:-2]  # skip the non-exact elements
            sage: TestSuite(L).run(elements=elts, skip=['_test_elements', '_test_associativity', '_test_distributivity', '_test_zero'])
            sage: L.category()
            Category of infinite commutative no zero divisors algebras over
             (euclidean domains and infinite enumerated sets and metric spaces)

            sage: L = LazyLaurentSeriesRing(QQ, 't')
            sage: L.category()
            Join of Category of complete discrete valuation fields
             and Category of commutative algebras over (number fields and quotient fields and metric spaces)
             and Category of infinite sets
            sage: L = LazyLaurentSeriesRing(ZZ['x,y'], 't')
            sage: L.category()
            Category of infinite commutative no zero divisors algebras over
             (unique factorization domains and commutative algebras over
              (euclidean domains and infinite enumerated sets and metric spaces)
              and infinite sets)
            sage: E.<x,y> = ExteriorAlgebra(QQ)
            sage: L = LazyLaurentSeriesRing(E, 't')  # not tested
        """
        self._sparse = sparse
        self._coeff_ring = base_ring
        # We always use the dense because our CS_exact is implemented densely
        self._laurent_poly_ring = LaurentPolynomialRing(base_ring, names)
        self._internal_poly_ring = self._laurent_poly_ring

        category = Algebras(base_ring.category())
        if base_ring in Fields():
            category &= CompleteDiscreteValuationFields()
        else:
            if "Commutative" in base_ring.category().axioms():
                category = category.Commutative()
            if base_ring in IntegralDomains():
                category &= IntegralDomains()

        if base_ring.is_zero():
            category = category.Finite()
        else:
            category = category.Infinite()

        Parent.__init__(self, base=base_ring, names=names, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LazyLaurentSeriesRing(GF(2), 'z')
            Lazy Laurent Series Ring in z over Finite Field of size 2
        """
        return "Lazy Laurent Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: latex(L)
            \Bold{F}_{2} (\!(z)\!)
        """
        from sage.misc.latex import latex
        return latex(self.base_ring()) + r"(\!({})\!)".format(self.variable_name())

    def characteristic(self):
        """
        Return the characteristic of this lazy power series ring, which
        is the same as the characteristic of its base ring.

        EXAMPLES::

            sage: L.<t> = LazyLaurentSeriesRing(ZZ)
            sage: L.characteristic()
            0
            sage: R.<w> = LazyLaurentSeriesRing(GF(11)); R
            Lazy Laurent Series Ring in w over Finite Field of size 11
            sage: R.characteristic()
            11

        """
        return self.base_ring().characteristic()

    def is_sparse(self):
        """
        Return whether ``self`` is sparse or not.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z', sparse=False)
            sage: L.is_sparse()
            False

            sage: L = LazyLaurentSeriesRing(ZZ, 'z', sparse=True)
            sage: L.is_sparse()
            True
        """
        return self._sparse

    @cached_method
    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.gen()
            z
            sage: L.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        if n != 0:
            raise IndexError("there is only one generator")
        R = self.base_ring()
        coeff_stream = Stream_exact([R.one()], self._sparse,
                                    constant=R.zero(), order=1)
        return self.element_class(self, coeff_stream)

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

        This is always 1.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.ngens()
            1
        """
        return 1

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.gens()
            (z,)
            sage: 1/(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
        """
        return tuple([self.gen(n) for n in range(self.ngens())])

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(GF(2))
            True
        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        R = self._laurent_poly_ring
        return R.has_coerce_map_from(S)

    def _coerce_map_from_base_ring(self):
        """
        Return a coercion map from the base ring of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: phi = L._coerce_map_from_base_ring()
            sage: phi(2)
            2
            sage: phi(2, valuation=-2)
            2*z^-2
            sage: phi(2, valuation=-2, constant=3, degree=1)
            2*z^-2 + 3*z + 3*z^2 + 3*z^3 + O(z^4)
        """
        # Return a DefaultConvertMap_unique; this can pass additional
        # arguments to _element_constructor_, unlike the map returned
        # by UnitalAlgebras.ParentMethods._coerce_map_from_base_ring.
        return self._generic_coerce_map(self.base_ring())

    def _an_element_(self):
        """
        Return a Laurent series in ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.an_element()
            z^-2 + 3*z^-1 + 2*z + z^2 + z^3 + z^4 + z^5 + O(z^6)
        """
        R = self.base_ring()
        coeff_stream = Stream_exact([R.an_element(), 3, 0, 2*R.an_element(), 1],
                                    self._sparse, order=-2, constant=R.one())
        return self.element_class(self, coeff_stream)

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.some_elements()
            [0, 1, z,
             -3*z^-4 + z^-3 - 12*z^-2 - 2*z^-1 - 10 - 8*z + z^2 + z^3,
             z^-2 + 3*z^-1 + 2*z + z^2 + z^3 + z^4 + z^5 + O(z^6),
             -2*z^-3 - 2*z^-2 + 4*z^-1 + 11 - z - 34*z^2 - 31*z^3 + O(z^4),
             4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + O(z^5)]

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L.some_elements()
            [0, 1, z,
             z^-4 + z^-3 + z^2 + z^3,
             z^-1 + z^2 + z^3 + z^4 + z^5 + O(z^6),
             1 + z + z^3 + z^4 + z^6 + O(z^7),
             z^-1 + z + z^3 + O(z^5)]

            sage: L = LazyLaurentSeriesRing(GF(3), 'z')
            sage: L.some_elements()
            [0, 1, z,
             z^-3 + z^-1 + 2 + z + z^2 + z^3,
             z^2 + z^3 + z^4 + z^5 + O(z^6),
             z^-3 + z^-2 + z^-1 + 2 + 2*z + 2*z^2 + 2*z^3 + O(z^4),
             z^-2 + z^-1 + z + z^2 + z^4 + O(z^5)]
        """
        z = self.gen()
        elts = [self.zero(), self.one(), z, (z-3)*(z**-2+2+z)**2, self.an_element(),
                (1 - 2*z**-3)/(1 - z + 3*z**2), self(lambda n: n**2, valuation=-2)]
        return elts

    @cached_method
    def one(self):
        r"""
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.one()
            1
        """
        R = self.base_ring()
        coeff_stream = Stream_exact([R.one()], self._sparse, constant=R.zero(), degree=1)
        return self.element_class(self, coeff_stream)

    @cached_method
    def zero(self):
        r"""
        Return the zero series.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, Stream_zero(self._sparse))

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for Lazy Laurent series.

        If no parameters are set, then the function returns a copy of
        the options dictionary.

        The ``options`` to Lazy Laurent series can be accessed as using
        :class:`LazyLaurentSeriesRing.options` of :class:`LazyLaurentSeriesRing`.

        @OPTIONS@

        EXAMPLES::

            sage: LLS.<z> = LazyLaurentSeriesRing(QQ)
            sage: LLS.options.display_length
            7
            sage: f = 1/(1-z)
            sage: f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: LLS.options.display_length = 10
            sage: f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + z^8 + z^9 + O(z^10)
            sage: g = LLS(lambda n: n^2, valuation=-2, degree=5, constant=42)
            sage: g
            4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + 42*z^5 + 42*z^6 + 42*z^7 + O(z^8)
            sage: LLS.options.constant_length = 1
            sage: g
            4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + 42*z^5 + O(z^6)
            sage: LazyLaurentSeriesRing.options._reset()
            sage: LazyLaurentSeriesRing.options.display_length
            7
        """
        NAME = 'LazyLaurentSeriesRing'
        module = 'sage.rings.lazy_series_ring'
        display_length = dict(default=7,
                              description='the number of coefficients to display from the valuation',
                              checker=lambda x: x in ZZ and x > 0)
        constant_length = dict(default=3,
                               description='the number of coefficients to display for nonzero constant series',
                               checker=lambda x: x in ZZ and x > 0)

    def series(self, coefficient, valuation, degree=None, constant=None):
        r"""
        Return a lazy Laurent series.

        INPUT:

        - ``coefficient`` -- Python function that computes coefficients or a list
        - ``valuation`` -- integer; approximate valuation of the series
        - ``degree`` -- (optional) integer
        - ``constant`` -- (optional) an element of the base ring

        Let the coefficient of index `i` mean the coefficient of the term
        of the series with exponent `i`.

        Python function ``coefficient`` returns the value of the coefficient
        of index `i` from input `s` and `i` where `s` is the series itself.

        Let ``valuation`` be `n`. All coefficients of index below `n` are zero.
        If ``constant`` is not specified, then the ``coefficient`` function is
        responsible to compute the values of all coefficients of index `\ge n`.
        If ``degree`` or ``constant`` is a pair `(c,m)`, then the ``coefficient``
        function is responsible to compute the values of all coefficients of
        index `\ge n` and `< m` and all the coefficients of index `\ge m`
        is the constant `c`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.series(lambda s, i: i, 5, (1,10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: def g(s, i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return s.coefficient(i - 1) + i
            sage: e = L.series(g, -5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + O(z^2)
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + O(z^12)
            sage: f.coefficient(10)
            0
            sage: f.coefficient(20)
            9
            sage: f.coefficient(30)
            -219

        Alternatively, the ``coefficient`` can be a list of elements of the
        base ring. Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero. ::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5); f
            z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2
            sage: g = L.series([1,3,5,7,9], 5, constant=-1); g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + O(z^13)
        """
        if valuation is not None and valuation not in ZZ:
            raise ValueError("the valuation must be an integer")

        if isinstance(constant, (list, tuple)):
            constant, degree = constant
        if isinstance(degree, (list, tuple)):
            constant, degree = degree

        if constant is not None:
            constant = self.base_ring()(constant)

        if isinstance(coefficient, (tuple, list)):
            if constant is None:
                constant = self.base_ring().zero()
            if degree is None:
                degree = valuation + len(coefficient)
            coeff_stream = Stream_exact(coefficient, self._sparse, order=valuation,
                                        constant=constant, degree=degree)
            return self.element_class(self, coeff_stream)

        if degree is not None and valuation > degree and constant:
            raise ValueError('inappropriate valuation')

        t = None
        t = self(lambda n: coefficient(t, n), valuation=valuation,
                 constant=constant, degree=degree)
        return t

    def _monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L._monomial(1, 3)
            z^3
            sage: L._monomial(2, -4)
            2*z^-4
        """
        return self._laurent_poly_ring(c).shift(n)

class LazyDirichletSeriesRing(LazySeriesRing):
    """
    Lazy Dirichlet series ring.

    INPUT:

    - ``base_ring`` -- base ring of this Dirichlet series ring
    - ``names`` -- name of the generator of this Dirichlet series ring
    - ``sparse`` -- (default: ``True``) whether this series is sparse or not

    EXAMPLES::

        sage: LazyDirichletSeriesRing(ZZ, 't')
        Lazy Dirichlet Series Ring in t over Integer Ring
    """
    Element = LazyDirichletSeries

    def __init__(self, base_ring, names, sparse=True, category=None):
        """
        Initialize the ring.

        TESTS::

            sage: L = LazyDirichletSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip=['_test_elements', '_test_associativity', '_test_distributivity', '_test_zero'])
        """
        if base_ring.characteristic() > 0:
            raise ValueError("positive characteristic not allowed for Dirichlet series")

        self._sparse = sparse
        self._coeff_ring = base_ring
        # TODO: it would be good to have something better than the symbolic ring
        self._laurent_poly_ring = SR
        self._internal_poly_ring = PolynomialRing(base_ring, names, sparse=True)

        category = Algebras(base_ring.category())
        if base_ring in IntegralDomains():
            category &= IntegralDomains()
        elif base_ring in Rings().Commutative():
            category = category.Commutative()
        category = category.Infinite()
        Parent.__init__(self, base=base_ring, names=names,
                        category=category)

    def _repr_(self):
        """
        String representation of this Dirichlet series ring.

        EXAMPLES::

            sage: LazyDirichletSeriesRing(QQbar, 'z')
            Lazy Dirichlet Series Ring in z over Algebraic Field
        """
        return "Lazy Dirichlet Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    def characteristic(self):
        """
        Return the characteristic of this lazy power series ring, which
        is the same as the characteristic of its base ring.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, "s")
            sage: L.characteristic()
            0
        """
        return self.base_ring().characteristic()

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(QQ)
            False
        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        return False

    def _coerce_map_from_base_ring(self):
        r"""
        Return a coercion map from the base ring of ``self``.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(QQ, 'z')
            sage: phi = L._coerce_map_from_base_ring()
            sage: phi(2)
            2
            sage: phi(2, valuation=2)
            2/2^z
            sage: phi(2, valuation=2, constant=4)
            2/2^z + 4/3^z + 4/4^z + 4/5^z + O(1/(6^z))
        """
        # Return a DefaultConvertMap_unique; this can pass additional
        # arguments to _element_constructor_, unlike the map returned
        # by UnitalAlgebras.ParentMethods._coerce_map_from_base_ring.
        return self._generic_coerce_map(self.base_ring())

    def _element_constructor_(self, x=None, valuation=None, degree=None, constant=None, coefficients=None):
        r"""
        Construct a Dirichlet series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a Dirichlet series
        - ``valuation`` -- integer (optional); integer; a lower bound for
          the exp of the valuation of the series
        - ``degree`` -- (optional) the degree when the series is ``constant``
        - ``constant`` -- (optional) the eventual constant of the series

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L(3)
            3
            sage: L(lambda i: i, constant=1, degree=6)
            1 + 2/2^z + 3/3^z + 4/4^z + 5/5^z + 1/(6^z) + 1/(7^z) + 1/(8^z) + O(1/(9^z))

            sage: X = L(constant=5, degree=3); X
            5/3^z + 5/4^z + 5/5^z + O(1/(6^z))
            sage: X.valuation()
            log(3)
            sage: e = L(moebius); e
            1 - 1/(2^z) - 1/(3^z) - 1/(5^z) + 1/(6^z) - 1/(7^z) + O(1/(8^z))

            sage: L([0], constant=1)
            1/(2^z) + 1/(3^z) + 1/(4^z) + O(1/(5^z))

            sage: L(constant=1)
            1 + 1/(2^z) + 1/(3^z) + O(1/(4^z))

            sage: L(lambda i: i, valuation=3)
            3/3^z + 4/4^z + 5/5^z + 6/6^z + 7/7^z + 8/8^z + 9/9^z + O(1/(10^z))

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], 4); f
            1/(4^z) + 2/5^z + 3/6^z + 4/7^z
            sage: g = L([1,3,5,7,9], 6, constant=-1); g
            1/(6^z) + 3/7^z + 5/8^z + 7/9^z + 9/10^z - 1/(11^z) - 1/(12^z) - 1/(13^z) + O(1/(14^z))

        TESTS::

            sage: L = LazyDirichletSeriesRing(GF(2), 'z')
            Traceback (most recent call last):
            ...
            ValueError: positive characteristic not allowed for Dirichlet series

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: D = LazyDirichletSeriesRing(QQ, 't')
            sage: D(L.one())
            1 + 1/(2^t) + 1/(3^t) + 1/(4^t) + 1/(5^t) + 1/(6^t) + 1/(7^t) + O(1/(8^t))

            sage: R.<z> = LaurentPolynomialRing(QQ)
            sage: D = LazyDirichletSeriesRing(QQ, 't')
            sage: D(coefficients=z+z^2)
            2 + 6/2^t + 12/3^t + 20/4^t + 30/5^t + 42/6^t + 56/7^t + O(1/(8^t))

            sage: s = D(lambda n: n)
            sage: D(s, valuation=2)
            1/(2^t) + 2/3^t + 3/4^t + 4/5^t + 5/6^t + 6/7^t + 7/8^t + O(1/(9^t))

            sage: Ds = LazyDirichletSeriesRing(ZZ, 's')
            sage: m = Ds(moebius, valuation=2); m
            -1/(2^s) - 1/(3^s) - 1/(5^s) + 1/(6^s) - 1/(7^s) + O(1/(9^s))
            sage: D = LazyDirichletSeriesRing(QQ, 't')
            sage: D(m)
            -1/(2^t) - 1/(3^t) - 1/(5^t) + 1/(6^t) - 1/(7^t) + O(1/(9^t))

        .. TODO::

            Add a method to make a copy of ``self._sparse``.
        """
        if isinstance(x, (list, tuple)):
            p = self._internal_poly_ring(x)
            if valuation is None:
                if not p:
                    valuation = 1 + len(x)
                    x = p
                else:
                    x = p.shift(1)
        else:
            if coefficients is not None:
                if valuation is None:
                    valuation = 1
                return super()._element_constructor_(x, valuation, degree, constant, coefficients)

            BR = self.base_ring()
            if x in BR:
                if valuation is None:
                    valuation = 1
                x = BR(x)

            elif not isinstance(x, LazyDirichletSeries):
                if valuation is None:
                    valuation = 1

                if isinstance(x, LazyModuleElement) or callable(x):
                    if coefficients is not None:
                        raise ValueError("coefficients must be None if x is provided")
                    coefficients = x
                    x = None

        if valuation is not None and (valuation not in ZZ or valuation <= 0):
            raise ValueError("the valuation must be a positive integer")

        return super()._element_constructor_(x, valuation, degree, constant, coefficients)

    def _an_element_(self):
        """
        Return a Dirichlet series in this ring.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.an_element()
            1/(4^z) + 1/(5^z) + 1/(6^z) + O(1/(7^z))
        """
        c = self.base_ring().an_element()
        return self.element_class(self, Stream_exact([], self._sparse, constant=c, order=4))

    @cached_method
    def one(self):
        """
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.one()
            1
            sage: ~L.one()
            1 + O(1/(8^z))
        """
        R = self.base_ring()
        return self.element_class(self, Stream_exact([R.one()], self._sparse, order=1))

    @cached_method
    def zero(self):
        """
        Return the zero series.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, Stream_zero(self._sparse))

    def _monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L._monomial(5, 3)
            5/3^z
        """
        try:
            L = self._laurent_poly_ring
            return L(c) * L(n) ** -L(self.variable_name())
        except (ValueError, TypeError):
            return '({})/{}^{}'.format(self.base_ring()(c), n, self.variable_name())

    options = LazyLaurentSeriesRing.options
