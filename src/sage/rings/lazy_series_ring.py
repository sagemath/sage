r"""
Lazy Series Rings

We provide lazy implementations for various `\NN`-graded rings.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :class:`LazyLaurentSeriesRing` | The ring of lazy Laurent series.
    :class:`LazyPowerSeriesRing` | The ring of (possibly multivariate) lazy Taylor series.
    :class:`LazyCompletionGradedAlgebra` | The completion of a graded algebra consisting of formal series.
    :class:`LazySymmetricFunctions` | The ring of (possibly multivariate) lazy symmetric functions.
    :class:`LazyDirichletSeriesRing` | The ring of lazy Dirichlet series.

.. SEEALSO::

    :class:`sage.rings.padics.generic_nodes.pAdicRelaxedGeneric`,
    :func:`sage.rings.padics.factory.ZpER`

.. WARNING::

    When the halting precision is infinite, the default for ``bool(f)``
    is ``True`` for any lazy series ``f`` that is not known to be zero.
    This could end up resulting in infinite loops::

        sage: L.<x> = LazyPowerSeriesRing(ZZ)
        sage: f = L(lambda n: 0, valuation=0)
        sage: 1 / f  # not tested - infinite loop

.. SEEALSO::

    The examples of :class:`LazyLaurentSeriesRing` contain a discussion
    about the different methods of comparisons the lazy series can use.

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version
- Tejasvi Chebrolu, Martin Rubey, Travis Scrimshaw (2021-08):
  refactored and expanded functionality
"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <ekwankyu@gmail.com>
#                     2022 Martin Rubey <martin.rubey at tuwien.ac.at>
#                     2022 Travis Scrimshaw <tcscrims at gmail.com>
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
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.rings import Rings
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains
from sage.categories.integral_domains import IntegralDomains
from sage.categories.fields import Fields
from sage.categories.complete_discrete_valuation import (CompleteDiscreteValuationFields,
                                                         CompleteDiscreteValuationRings)

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.lazy_series import (LazyModuleElement,
                                    LazyLaurentSeries,
                                    LazyPowerSeries,
                                    LazyPowerSeries_gcd_mixin,
                                    LazyCompletionGradedAlgebraElement,
                                    LazySymmetricFunction,
                                    LazyDirichletSeries)
from sage.structure.global_options import GlobalOptions

from sage.data_structures.stream import (
    Stream_zero,
    Stream_function,
    Stream_iterator,
    Stream_exact,
    Stream_uninitialized,
    Stream_taylor
)

from types import GeneratorType


class LazySeriesRing(UniqueRepresentation, Parent):
    """
    Abstract base class for lazy series.
    """
    # This will never be called directly (as it is an ABC), but we copy it
    #   for use in other subclasses.
    @staticmethod
    def __classcall_private__(cls, base_ring, names, sparse=True, *args, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: Lp = LazyLaurentSeriesRing(QQ, 'z')
            sage: L is Lp
            True
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        return super().__classcall__(cls, base_ring, names, sparse, *args, **kwds)

    def _element_constructor_(self, x=None, valuation=None, degree=None, constant=None, coefficients=None):
        r"""
        Construct a lazy series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a series
        - ``valuation`` -- integer (optional); a lower bound for the valuation
          of the series
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

        .. WARNING::

            If ``x`` is provided as a list, any trailing zeros are
            ignored, because ``x`` is immediately converted into a
            polynomial.

        EXAMPLES:

        If ``x`` can be converted into an element of the underlying
        Laurent polynomial ring, we do this::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1

        In particular, ``x`` can be a Laurent polynomial::

            sage: P.<x> = LaurentPolynomialRing(QQ)
            sage: p = x^-2 + 3*x^3
            sage: L.<x> = LazyLaurentSeriesRing(ZZ)
            sage: L(p)
            x^-2 + 3*x^3

            sage: L(p, valuation=0)
            1 + 3*x^5

            sage: L(p, valuation=1)
            x + 3*x^6

        If ``x`` is callable, its evaluation at the integers,
        beginning at ``valuation``, defines the coefficients of the series::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L(lambda i: i, valuation=5, constant=1, degree=10)
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)
            sage: L(lambda i: i, valuation=5, constant=(1, 10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: def g(i):
            ....:     if i < 0:
            ....:         return 1
            ....:     return 1 + sum(range(i + 1))
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

        We can omit ``x``, when defining a series with constant coefficients::

            sage: X = L(constant=5, degree=2); X
            5*z^2 + 5*z^3 + 5*z^4 + O(z^5)
            sage: X.valuation()
            2

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

        If ``x`` is explicitly passed as ``None``, the resulting
        series is undefined.  This can be used to define it
        implicitly, see
        :meth:`sage.rings.lazy_series.LazyModuleElement.define`::

            sage: f = L(None, valuation=-1)
            sage: f.define(z^-1 + z^2*f^2)
            sage: f
            z^-1 + 1 + 2*z + 5*z^2 + 14*z^3 + 42*z^4 + 132*z^5 + O(z^6)

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
            sage: L(lambda n: 1/factorial(n), valuation=0)                              # needs sage.symbolic
            (1 + 1/2/2^s + 1/6/3^s + 1/24/4^s + 1/120/5^s + 1/720/6^s + 1/5040/7^s + O(1/(8^s)))

        We can also specify that the given function should be
        interpreted as the coefficients of the Laurent series::

            sage: L(coefficients=lambda n: 1/factorial(n), valuation=0)                 # needs sage.symbolic
            1 + z + 1/2*z^2 + 1/6*z^3 + 1/24*z^4 + 1/120*z^5 + 1/720*z^6 + O(z^7)

        When the argument ``x`` is callable and not convertible into
        an element of the underlying ring of series of finite
        support, it is evaluated at integers to compute the
        coefficients::

            sage: R.<q> = QQ[]
            sage: D = LazyDirichletSeriesRing(ZZ, 't')
            sage: D(1+2*q)                                                              # needs sage.symbolic
            3 + 5/2^t + 7/3^t + 9/4^t + 11/5^t + 13/6^t + 15/7^t + O(1/(8^t))

        In this example, the Dirichlet series ``m`` is considered as an
        element in the base ring::

            sage: m = D(moebius)
            sage: s = L(m, valuation=0)                                                 # needs sage.symbolic
            sage: s[0]
            1 - 1/(2^s) - 1/(3^s) - 1/(5^s) + 1/(6^s) - 1/(7^s) + O(1/(8^s))
            sage: s[1]
            0

        Converting various series from a univariate power series::

            sage: # needs sage.rings.finite_rings
            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: R = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.has_coerce_map_from(R)
            True
            sage: L(R(lambda n: n))
            z + z^3 + z^5 + z^7 + O(z^8)
            sage: L(R([2,4,6])) == L.zero()
            True
            sage: L(R([2,4,6], valuation=2, constant=4)) == L.zero()
            True
            sage: L(R([2,4,6], valuation=2, constant=5))
            z^5 + z^6 + z^7 + O(z^8)
            sage: L(R([2,3,4], valuation=2, constant=4))
            z^3
            sage: L(R([2,3,4], valuation=2, constant=5))
            z^3 + z^5 + z^6 + z^7 + O(z^8)

        Can only convert from known to be constant multivariate power series::

            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: R.<x,y> = LazyPowerSeriesRing(QQ)
            sage: L(R(2))
            2
            sage: L(R.zero())
            0
            sage: L(x)
            Traceback (most recent call last):
            ...
            ValueError: unable to convert ...
            sage: L(1 / (1 - x - y))
            Traceback (most recent call last):
            ...
            ValueError: unable to convert ...
            sage: P.<x,y> = QQ[]
            sage: f = R(lambda n: (x+y)^n if n == 0 else P.zero()); f
            1 + O(x,y)^7
            sage: L(f)
            Traceback (most recent call last):
            ...
            ValueError: unable to convert ...

        Converting from the corresponding rational functions::

            sage: L = LazyLaurentSeriesRing(QQ, 't')
            sage: tt = L.gen()
            sage: R.<t> = LaurentPolynomialRing(QQ)
            sage: f = (1 + t) / (1 + t + t^2); f
            (t + 1)/(t^2 + t + 1)
            sage: f.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: L(f)
            1 - t^2 + t^3 - t^5 + t^6 + O(t^7)
            sage: L(f) == (1 + tt) / (1 + tt + tt^2)
            True
            sage: f = (3 + t) / (t^3 - t^5); f
            (-t - 3)/(t^5 - t^3)
            sage: f.parent()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
            sage: L(f)
            3*t^-3 + t^-2 + 3*t^-1 + 1 + 3*t + t^2 + 3*t^3 + O(t^4)
            sage: L(f) - (3 + tt) / (tt^3 - tt^5)
            O(t^4)

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
            z^-4 + z^-3 + z^-2 + O(z^-1)
            sage: L(z^-3/(1-z), valuation=-4)
            z^-4 + z^-3 + z^-2 + O(z^-1)
            sage: L(z^3/(1-z), valuation=-4)
            z^-4 + z^-3 + z^-2 + O(z^-1)

            sage: L(z^3/(1-z), valuation=0)
            1 + z + z^2 + O(z^3)

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
            sage: D(E([1,2,3]))                                                         # needs sage.symbolic
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

        We support passing a generator::

            sage: L(filter(is_odd, NN), -3)
            z^-3 + 3*z^-2 + 5*z^-1 + 7 + 9*z + 11*z^2 + 13*z^3 + O(z^4)
        """
        if valuation is not None and valuation not in ZZ:
            raise ValueError("the valuation must be an integer")

        if x is None and coefficients is None:
            if valuation is None:
                raise ValueError("the valuation must be specified")
            return self.element_class(self, Stream_uninitialized(valuation))

        # WARNING: if x is not explicitly specified as None, it is
        # set to 0 by Parent.__call__
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
                if not x:
                    coeff_stream = Stream_exact([], order=degree, constant=constant)
                    return self.element_class(self, coeff_stream)
                initial_coefficients = [x[i] for i in range(x.valuation(), x.degree() + 1)]
                coeff_stream = Stream_exact(initial_coefficients,
                                            order=x.valuation(), degree=degree, constant=constant)
                return self.element_class(self, coeff_stream)

            # Handle when it is a lazy series
            if isinstance(x, self.Element):
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
                    coeff_stream = Stream_exact([], order=degree, constant=constant)
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
                    coeff_stream = Stream_exact(initial_coefficients,
                                                order=valuation, degree=degree, constant=constant)
                    return self.element_class(self, coeff_stream)

                # We are just possibly shifting the result
                ret = self.element_class(self, x._coeff_stream)
                if valuation is None:
                    return ret
                return ret.shift(valuation - x._coeff_stream.order())

            # Handle when it is a power series
            if isinstance(x, LazyPowerSeries):
                stream = x._coeff_stream
                if isinstance(stream, Stream_zero):
                    return self.zero()
                elif isinstance(stream, Stream_exact):
                    BR = self.base_ring()
                    if x.parent()._arity != 1:
                        # Special case for constant series
                        if stream._degree == 1:
                            return self(BR(stream[0]))
                    else:
                        coeffs = [BR(val) for val in stream._initial_coefficients]
                        valuation = stream._approximate_order
                        for i, c in enumerate(coeffs):
                            if c:
                                valuation += i
                                coeffs = coeffs[i:]
                                break
                        else:
                            valuation += len(coeffs)
                            coeffs = []
                        return self(coeffs,
                                    degree=stream._degree,
                                    constant=BR(stream._constant),
                                    valuation=valuation)
                elif x.parent()._arity == 1:
                    return self.element_class(self, stream)
                raise ValueError(f"unable to convert {x} into {self}")

            # Check if we can realize the input as a rational function
            try:
                FF = self._laurent_poly_ring.fraction_field()
                x = FF(x)
            except (TypeError, ValueError, AttributeError):
                pass
            else:
                return self(x.numerator()) / self(x.denominator())

        else:
            x = coefficients

        if callable(x) or isinstance(x, (GeneratorType, map, filter)):
            if valuation is None:
                raise ValueError("the valuation must be specified")
            if degree is None:
                if constant is not None:
                    raise ValueError("constant may only be specified if the degree is specified")
                if callable(x):
                    coeff_stream = Stream_function(lambda i: BR(x(i)), self._sparse, valuation)
                else:
                    coeff_stream = Stream_iterator(map(BR, _skip_leading_zeros(x)), valuation)
                return self.element_class(self, coeff_stream)

            # degree is not None
            if constant is None:
                constant = BR.zero()
            if callable(x):
                p = [BR(x(i)) for i in range(valuation, degree)]
            else:
                p = [BR(c) for c, _ in zip(_skip_leading_zeros(x), range(valuation, degree))]
            if not any(p) and not constant:
                return self.zero()
            coeff_stream = Stream_exact(p, order=valuation, constant=constant, degree=degree)
            return self.element_class(self, coeff_stream)

        raise ValueError(f"unable to convert {x} into {self}")

    def undefined(self, valuation=None):
        r"""
        Return an uninitialized series.

        INPUT:

        - ``valuation`` -- integer; a lower bound for the valuation of the series

        Power series can be defined recursively (see
        :meth:`sage.rings.lazy_series.LazyModuleElement.define` for
        more examples).

        .. SEEALSO::

            :meth:`sage.rings.padics.generic_nodes.pAdicRelaxedGeneric.unknown`

        EXAMPLES::

            sage: L.<z> = LazyPowerSeriesRing(QQ)
            sage: s = L.undefined(1)
            sage: s.define(z + (s^2+s(z^2))/2)
            sage: s
            z + z^2 + z^3 + 2*z^4 + 3*z^5 + 6*z^6 + 11*z^7 + O(z^8)

        Alternatively::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = L(None, valuation=-1)
            sage: f.define(z^-1 + z^2*f^2)
            sage: f
            z^-1 + 1 + 2*z + 5*z^2 + 14*z^3 + 42*z^4 + 132*z^5 + O(z^6)
        """
        if valuation is None:
            valuation = self._minimal_valuation
        coeff_stream = Stream_uninitialized(valuation)
        return self.element_class(self, coeff_stream)

    unknown = undefined

    class options(GlobalOptions):
        r"""
        Set and display the options for lazy series.

        If no parameters are set, then the function returns a copy of
        the options dictionary.

        The ``options`` to lazy series can be accessed as using
        :class:`LazySeriesRing.options`.

        @OPTIONS@

        EXAMPLES::

            sage: LLS.<z> = LazyLaurentSeriesRing(QQ)
            sage: LLS.options
            Current options for lazy series rings
              - constant_length:   3
              - display_length:    7
              - halting_precision: None
              - secure:            False

            sage: LLS.options.display_length
            7
            sage: f = 1 / (1 + z)
            sage: f
            1 - z + z^2 - z^3 + z^4 - z^5 + z^6 + O(z^7)
            sage: LLS.options.display_length = 10
            sage: f
            1 - z + z^2 - z^3 + z^4 - z^5 + z^6 - z^7 + z^8 - z^9 + O(z^10)
            sage: g = LLS(lambda n: n^2, valuation=-2, degree=5, constant=42)
            sage: g
            4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + 42*z^5 + 42*z^6 + 42*z^7 + O(z^8)
            sage: h = 1 / (1 - z)  # This is exact
            sage: h
            1 + z + z^2 + O(z^3)
            sage: LLS.options.constant_length = 1
            sage: g
            4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + 42*z^5 + O(z^6)
            sage: h
            1 + O(z)
            sage: LazyLaurentSeriesRing.options._reset()
            sage: LazyLaurentSeriesRing.options.display_length
            7
        """
        NAME = 'lazy series rings'
        module = 'sage.rings.lazy_series_ring'
        display_length = dict(default=7,
                              description='the number of coefficients to display from the valuation',
                              checker=lambda x: x in ZZ and x > 0)
        constant_length = dict(default=3,
                               description='the number of coefficients to display for nonzero constant series',
                               checker=lambda x: x in ZZ and x > 0)
        halting_precision = dict(default=None,
                                 description='the number of coefficients, beginning with the approximate valuation, to check in equality tests',
                                 checker=lambda x: x is None or x in ZZ and x > 0)
        secure = dict(default=False,
                      description='whether to raise an error when a comparison is unknown',
                      checker=lambda x: x is True or x is False)

    @cached_method
    def one(self):
        r"""
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.one()
            1

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.one()
            1

            sage: m = SymmetricFunctions(ZZ).m()                                        # needs sage.modules
            sage: L = LazySymmetricFunctions(m)                                         # needs sage.modules
            sage: L.one()                                                               # needs sage.modules
            m[]
        """
        R = self.base_ring()
        coeff_stream = Stream_exact([R.one()], constant=R.zero(), order=0)
        return self.element_class(self, coeff_stream)

    @cached_method
    def zero(self):
        r"""
        Return the zero series.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.zero()
            0

            sage: s = SymmetricFunctions(ZZ).s()                                        # needs sage.modules
            sage: L = LazySymmetricFunctions(s)                                         # needs sage.modules
            sage: L.zero()                                                              # needs sage.modules
            0

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.zero()
            0

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, Stream_zero())

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

            sage: R.<x, y> = LazyPowerSeriesRing(GF(7)); R
            Multivariate Lazy Taylor Series Ring in x, y over Finite Field of size 7
            sage: R.characteristic()
            7

            sage: L = LazyDirichletSeriesRing(ZZ, "s")
            sage: L.characteristic()
            0
        """
        return self.base_ring().characteristic()

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(GF(2))
            True
            sage: R = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.has_coerce_map_from(R)
            True

            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: R = LazyPowerSeriesRing(QQ, 'z')
            sage: L.has_coerce_map_from(R)
            True
            sage: R = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.has_coerce_map_from(R)
            True
            sage: R = LazyPowerSeriesRing(ZZ['t'], 'z')
            sage: L.has_coerce_map_from(R)
            False

            sage: L = LazyPowerSeriesRing(GF(2), 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(GF(2))
            True

            sage: # needs sage.modules sage.rings.finite_rings
            sage: s = SymmetricFunctions(GF(2)).s()
            sage: L = LazySymmetricFunctions(s)
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(GF(2))
            True
        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        R = self._laurent_poly_ring
        if R.has_coerce_map_from(S):
            return True

        if (isinstance(S, LazySeriesRing)
            and self._laurent_poly_ring.has_coerce_map_from(S._laurent_poly_ring)):
            return True

        return None

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

            sage: L = LazyDirichletSeriesRing(QQ, 'z')
            sage: phi = L._coerce_map_from_base_ring()
            sage: m = phi(2)
            sage: m                                                                     # needs sage.symbolic
            2
            sage: m = phi(2, valuation=2)
            sage: m                                                                     # needs sage.symbolic
            2/2^z
            sage: m = phi(2, valuation=2, constant=4)
            sage: m                                                                     # needs sage.symbolic
            2/2^z + 4/3^z + 4/4^z + 4/5^z + O(1/(6^z))
        """
        # Return a DefaultConvertMap_unique; this can pass additional
        # arguments to _element_constructor_, unlike the map returned
        # by UnitalAlgebras.ParentMethods._coerce_map_from_base_ring.
        return self._generic_coerce_map(self.base_ring())

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

    def is_exact(self):
        """
        Return if ``self`` is exact or not.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.is_exact()
            True
            sage: L = LazyLaurentSeriesRing(RR, 'z')
            sage: L.is_exact()
            False
        """
        return self.base_ring().is_exact()

    def prod(self, f, a=None, b=infinity, add_one=False):
        r"""
        The product of elements of ``self``.

        INPUT:

        - ``f`` -- list (or iterable) of elements of ``self``
        - ``a``, ``b`` -- optional arguments
        - ``add_one`` -- (default: ``False``) if ``True``, then converts a
          lazy series `p_i` from ``args`` into `1 + p_i` for the product

        If ``a`` and ``b`` are both integers, then this returns the product
        `\prod_{i=a}^b f(i)`, where `f(i) = p_i` if ``add_one=False`` or
        `f(i) = 1 + p_i` otherwise. If ``b`` is not specified, then we consider
        `b = \infty`. Note this corresponds to the Python ``range(a, b+1)``.

        If `a` is any other iterable, then this returns the product
        `\prod_{i \in a} f(i)`, where `f(i) = p_i` if ``add_one=False`` or
        `f(i) = 1 + p_i`.

        .. NOTE::

            For infinite products, it is faster to use ``add_one=True`` since
            the implementation is based on `p_i` in `\prod_i (1 + p_i)`.

        .. WARNING::

            When ``f`` is an infinite generator, then the first argument
            ``a`` must be ``True``. Otherwise this will loop forever.

        .. WARNING::

            For an *infinite* product of the form `\prod_i (1 + p_i)`,
            if `p_i = 0`, then this will loop forever.

        EXAMPLES::

            sage: L.<t> = LazyLaurentSeriesRing(QQ)
            sage: euler = L.prod(lambda n: 1 - t^n, PositiveIntegers())
            sage: euler
            1 - t - t^2 + t^5 + O(t^7)
            sage: 1 / euler
            1 + t + 2*t^2 + 3*t^3 + 5*t^4 + 7*t^5 + 11*t^6 + O(t^7)
            sage: euler - L.euler()
            O(t^7)
            sage: L.prod(lambda n: -t^n, 1, add_one=True)
            1 - t - t^2 + t^5 + O(t^7)

            sage: L.prod((1 - t^n for n in PositiveIntegers()), True)
            1 - t - t^2 + t^5 + O(t^7)
            sage: L.prod((-t^n for n in PositiveIntegers()), True, add_one=True)
            1 - t - t^2 + t^5 + O(t^7)

            sage: L.prod((1 + t^(n-3) for n in PositiveIntegers()), True)
            2*t^-3 + 4*t^-2 + 4*t^-1 + 4 + 6*t + 10*t^2 + 16*t^3 + O(t^4)

            sage: L.prod(lambda n: 2 + t^n, -3, 5)
            96*t^-6 + 240*t^-5 + 336*t^-4 + 840*t^-3 + 984*t^-2 + 1248*t^-1
             + 1980 + 1668*t + 1824*t^2 + 1872*t^3 + 1782*t^4 + 1710*t^5
             + 1314*t^6 + 1122*t^7 + 858*t^8 + 711*t^9 + 438*t^10 + 282*t^11
             + 210*t^12 + 84*t^13 + 60*t^14 + 24*t^15
            sage: L.prod(lambda n: t^n / (1 + abs(n)), -2, 2, add_one=True)
            1/3*t^-3 + 5/6*t^-2 + 13/9*t^-1 + 25/9 + 13/9*t + 5/6*t^2 + 1/3*t^3
            sage: L.prod(lambda n: t^-2 + t^n / n, -4, -2)
            1/24*t^-9 - 1/8*t^-8 - 1/6*t^-7 + 1/2*t^-6

            sage: D = LazyDirichletSeriesRing(QQ, "s")
            sage: D.prod(lambda p: (1+D(1, valuation=p)).inverse(), Primes())
            1 - 1/(2^s) - 1/(3^s) + 1/(4^s) - 1/(5^s) + 1/(6^s) - 1/(7^s) + O(1/(8^s))

            sage: D.prod(lambda p: D(1, valuation=p), Primes(), add_one=True)
            1 + 1/(2^s) + 1/(3^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + O(1/(8^s))
        """
        if a is None:
            if add_one:
                return super().prod(self.one() + g for g in f)
            return super().prod(f)

        if a is True:
            it = f
        elif a in ZZ:
            if b != infinity:
                if add_one:
                    return super().prod(self.one() + f(i) for i in range(a, b+1))
                return super().prod(f(i) for i in range(a, b+1))
            from sage.sets.non_negative_integers import NonNegativeIntegers
            it = (f(i+a) for i in NonNegativeIntegers())
        else:
            it = (f(i) for i in a)

        # NOTE: We must have a new variable name for each new iterator
        if not add_one:
            data = (g - self.one() for g in it)
        else:
            data = it

        from sage.data_structures.stream import Stream_infinite_product
        coeff_stream = Stream_infinite_product(data)
        return self.element_class(self, coeff_stream)

    def sum(self, f, a=None, b=infinity):
        r"""
        The sum of elements of ``self``.

        INPUT:

        - ``f`` -- list (or iterable or function) of elements of ``self``
        - ``a``, ``b`` -- optional arguments

        If ``a`` and ``b`` are both integers, then this returns the sum
        `\sum_{i=a}^b f(i)`. If ``b`` is not specified, then we consider
        `b = \infty`. Note this corresponds to the Python ``range(a, b+1)``.

        If `a` is any other iterable, then this returns the sum
        `\sum{i \in a} f(i)`.

        .. WARNING::

            When ``f`` is an infinite generator, then the first argument
            ``a`` must be ``True``. Otherwise this will loop forever.

        .. WARNING::

            For an *infinite* sum of the form `\sum_i s_i`,
            if `s_i = 0`, then this will loop forever.

        EXAMPLES::

            sage: L.<t> = LazyLaurentSeriesRing(QQ)
            sage: L.sum(lambda n: t^n / (n+1), PositiveIntegers())
            1/2*t + 1/3*t^2 + 1/4*t^3 + 1/5*t^4 + 1/6*t^5 + 1/7*t^6 + 1/8*t^7 + O(t^8)

            sage: L.<z> = LazyPowerSeriesRing(QQ)
            sage: T = L.undefined(1)
            sage: D = L.undefined(0)
            sage: H = L.sum(lambda k: T(z^k)/k, 2)
            sage: T.define(z*exp(T)*D)
            sage: D.define(exp(H))
            sage: T
            z + z^2 + 2*z^3 + 4*z^4 + 9*z^5 + 20*z^6 + 48*z^7 + O(z^8)
            sage: D
            1 + 1/2*z^2 + 1/3*z^3 + 7/8*z^4 + 11/30*z^5 + 281/144*z^6 + O(z^7)

        We verify the Rogers-Ramanujan identities up to degree 100::

            sage: L.<q> = LazyPowerSeriesRing(QQ)
            sage: Gpi = L.prod(lambda k: -q^(1+5*k), 0, oo, add_one=True)
            sage: Gpi *= L.prod(lambda k: -q^(4+5*k), 0, oo, add_one=True)
            sage: Gp = 1 / Gpi
            sage: G = L.sum(lambda n: q^(n^2) / prod(1 - q^(k+1) for k in range(n)), 0, oo)
            sage: G - Gp
            O(q^7)
            sage: all(G[k] == Gp[k] for k in range(100))
            True

            sage: Hpi = L.prod(lambda k: -q^(2+5*k), 0, oo, add_one=True)
            sage: Hpi *= L.prod(lambda k: -q^(3+5*k), 0, oo, add_one=True)
            sage: Hp = 1 / Hpi
            sage: H = L.sum(lambda n: q^(n^2+n) / prod(1 - q^(k+1) for k in range(n)), 0, oo)
            sage: H - Hp
            O(q^7)
            sage: all(H[k] == Hp[k] for k in range(100))
            True

        ::

            sage: D = LazyDirichletSeriesRing(QQ, "s")
            sage: D.sum(lambda p: D(1, valuation=p), Primes())
            1/(2^s) + 1/(3^s) + 1/(5^s) + 1/(7^s) + O(1/(9^s))
        """
        if a is None:
            return super().sum(f)

        if a is True:
            it = f
        elif a in ZZ:
            if b != infinity:
                return super().sum(f(i) for i in range(a, b+1))
            from sage.sets.non_negative_integers import NonNegativeIntegers
            it = (f(i+a) for i in NonNegativeIntegers())
        else:
            it = (f(i) for i in a)

        from sage.data_structures.stream import Stream_infinite_sum
        coeff_stream = Stream_infinite_sum(it)
        return self.element_class(self, coeff_stream)

    def _test_invert(self, **options):
        """
        Test multiplicative inversion of elements of ``self``.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES::

            sage: LazyLaurentSeriesRing.options.halting_precision(5)
            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: L._test_invert()
            sage: LazyLaurentSeriesRing.options._reset()  # reset the options

        .. SEEALSO::

            :class:`TestSuite`
        """
        tester = self._tester(**options)

        elements = tester.some_elements()
        for x in elements:
            # because of laziness, creating the inverse of x should
            # always succeed except if the series is 'exact'
            if not x.is_unit():
                continue
            y = ~x
            e = y * x
            tester.assertFalse(x.is_zero(), "zero should not be invertible")
            tester.assertTrue(e.is_one(), "an element (%s) times its inverse should be 1" % x)
            tester.assertEqual(y.valuation(), -x.valuation(), "the valuation of the inverse should be the negative of the valuation of the element (%s)" % x)

    def _test_div(self, **options):
        r"""
        Test division of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES::

            sage: LazyLaurentSeriesRing.options.halting_precision(5)
            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: L._test_div()
            sage: LazyLaurentSeriesRing.options._reset()  # reset the options

        .. SEEALSO::

            :class:`TestSuite`
        """
        from sage.misc.misc import some_tuples
        tester = self._tester(**options)

        elements = list(tester.some_elements())
        for x, y in some_tuples(elements, 2, tester._max_runs):
            # because of laziness, creating the inverse of x should
            # always succeed except if the series is 'exact'
            if not y.is_unit():
                continue
            z = x / y
            xx = z * y
            try:
                v_z = z.valuation()
            except Exception as error:
                raise ValueError("could not compute the valuation of the quotient (%s)/(%s): %s" % (x, y, error))
            else:
                v_x = x.valuation()
                v_y = y.valuation()
                tester.assertEqual(v_z, v_x - v_y, "the valuation of the quotient should be the difference of the valuations of the elements (%s and %s)" % (x, y))
                tester.assertEqual(xx, x, "the element (%s) should be the quotient times the divisor (%s)" % (x, y))

    def _test_revert(self, **options):
        """
        Test compositional inverse of elements of this ring.

        INPUT:

        - ``options`` -- any keyword arguments accepted by :meth:`_tester`

        EXAMPLES::

            sage: LazyLaurentSeriesRing.options.halting_precision(5)
            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: L._test_revert()
            sage: LazyLaurentSeriesRing.options._reset()

        .. SEEALSO::

            :class:`TestSuite`
        """
        if not hasattr(self.element_class, "revert") or self._arity != 1:
            return
        tester = self._tester(**options)

        elements = tester.some_elements()
        count = 0
        for x in elements:
            # because of laziness, creating the compositional inverse
            # of x should always succeed, except if the series is
            # 'exact' or if it has negative valuation
            vx = x.valuation()
            if (vx != 1
                and not (isinstance(x._coeff_stream, Stream_exact)
                         and ((vx == 0
                               and x._coeff_stream._degree == 2
                               and not x._coeff_stream._constant)
                              or (vx == -1
                                  and x._coeff_stream._degree == 0
                                  and not x._coeff_stream._constant)))):
                continue
            try:
                y = x.revert()
            except Exception as error:
                raise AssertionError("compositional inverse of %s should exist: %s" % (x, error))
            try:
                vy = y.valuation()
                _ = y[vy]
            except NotImplementedError:
                pass
            except (ValueError, TypeError):
                tester.assertFalse(vx == 1 and x[vx].is_unit(),
                                   ("the series %s should be reversible "
                                    "- its valuation is one and its leading coefficient is a unit") % x)
            else:
                count += 1
                e1 = y(x)
                e2 = x(y)
                tester.assertEqual(e1, e2, "y(x) and x(y) differ for x = %s and y = %s" % (x, y))
                # tester.assertEqual(e1, self.gen())
        # we want to test at least 2 elements
        tester.assertGreater(count, 1, msg="only %s elements in %s.some_elements() have a compositional inverse" % (count, self))


class LazyLaurentSeriesRing(LazySeriesRing):
    r"""
    The ring of lazy Laurent series.

    The ring of Laurent series over a ring with the usual arithmetic
    where the coefficients are computed lazily.

    INPUT:

    - ``base_ring`` -- base ring
    - ``names`` -- name of the generator
    - ``sparse`` -- boolean (default: ``True``); whether the implementation of
      the series is sparse or not

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(QQ)
        sage: 1 / (1 - z)
        1 + z + z^2 + O(z^3)
        sage: 1 / (1 - z) == 1 / (1 - z)
        True
        sage: L in Fields
        True

    Lazy Laurent series ring over a finite field::

        sage: # needs sage.rings.finite_rings
        sage: L.<z> = LazyLaurentSeriesRing(GF(3)); L
        Lazy Laurent Series Ring in z over Finite Field of size 3
        sage: e = 1 / (1 + z)
        sage: e.coefficient(100)
        1
        sage: e.coefficient(100).parent()
        Finite Field of size 3

    Series can be defined by specifying a coefficient function
    and a valuation::

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

    We can also specify a polynomial or the initial coefficients.
    Additionally, we may specify that all coefficients are equal to a
    given constant, beginning at a given degree::

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

    We can truncate a series, shift its coefficients, or replace all
    coefficients beginning with a given degree by a constant::

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

    Power series can be defined recursively (see
    :meth:`sage.rings.lazy_series.LazyModuleElement.define` for
    more examples)::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: s = L.undefined(valuation=0)
        sage: s.define(1 + z*s^2)
        sage: s
        1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + O(z^7)

    By default, any two series ``f`` and ``g`` that are not known to
    be equal are considered to be different::

        sage: f = L(lambda n: 0, valuation=0)
        sage: f == 0
        False

        sage: f = L(constant=1, valuation=0).derivative(); f
        1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + O(z^7)
        sage: g = L(lambda n: (n+1), valuation=0); g
        1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + O(z^7)
        sage: f == g
        False

    .. WARNING::

        We have imposed that ``(f == g) == not (f != g)``, and so
        ``f != g`` returning ``True`` might not mean that the two
        series are actually different::

            sage: f = L(lambda n: 0, valuation=0)
            sage: g = L.zero()
            sage: f != g
            True

        This can be verified by :meth:`~sage.rings.lazy_series.is_nonzero()`,
        which only returns ``True`` if the series is known to be nonzero::

            sage: (f - g).is_nonzero()
            False

    The implementation of the ring can be either be a sparse or a dense one.
    The default is a sparse implementation::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: L.is_sparse()
        True
        sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=False)
        sage: L.is_sparse()
        False

    We additionally provide two other methods of performing comparisons.
    The first is raising a :exc:`ValueError` and the second uses a check
    up to a (user set) finite precision. These behaviors are set using the
    options ``secure`` and ``halting_precision``. In particular,
    this applies to series that are not specified by a finite number
    of initial coefficients and a constant for the remaining coefficients.
    Equality checking will depend on the coefficients which have
    already been computed. If this information is not enough to
    check that two series are different, then if ``L.options.secure``
    is set to ``True``, then we raise a :exc:`ValueError`::

        sage: L.options.secure = True
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
        sage: f2a = f + f
        sage: f2 == f2a
        Traceback (most recent call last):
        ...
        ValueError: undecidable
        sage: zf = L(lambda n: 0, valuation=0)
        sage: zf == 0
        Traceback (most recent call last):
        ...
        ValueError: undecidable

    For boolean checks, an error is raised when it is not known to be nonzero::

        sage: bool(zf)
        Traceback (most recent call last):
        ...
        ValueError: undecidable

    If the halting precision is set to a finite number `p` (for unlimited
    precision, it is set to ``None``), then it will check up to `p` values
    from the current position::

        sage: L.options.halting_precision = 20
        sage: f2 = f * 2  # currently no coefficients computed
        sage: f3 = f * 3  # currently no coefficients computed
        sage: f2 == f3
        False
        sage: f2a = f + f
        sage: f2 == f2a
        True
        sage: zf = L(lambda n: 0, valuation=0)
        sage: zf == 0
        True

    TESTS:

    We reset the options::

        sage: L.options._reset()
    """
    Element = LazyLaurentSeries

    # Follow the "generic" normalization
    __classcall_private__ = LazySeriesRing.__classcall_private__

    def __init__(self, base_ring, names, sparse=True, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: LazyLaurentSeriesRing.options.halting_precision(12)

            sage: L = LazyLaurentSeriesRing(ZZ, 't')
            sage: TestSuite(L).run()
            sage: L.category()
            Category of infinite commutative no zero divisors algebras over
             (Dedekind domains and euclidean domains
              and noetherian rings
              and infinite enumerated sets and metric spaces)

            sage: L = LazyLaurentSeriesRing(QQ, 't')
            sage: TestSuite(L).run()
            sage: L.category()
            Join of Category of complete discrete valuation fields
             and Category of commutative algebras over (number fields and quotient fields and metric spaces)
             and Category of infinite sets

            sage: L = LazyLaurentSeriesRing(ZZ['x, y'], 't')
            sage: TestSuite(L).run()                                                    # needs sage.libs.singular
            sage: L.category()
            Category of infinite commutative no zero divisors algebras over
             (unique factorization domains and algebras with basis over
              (Dedekind domains and euclidean domains
               and noetherian rings
               and infinite enumerated sets and metric spaces)
              and commutative algebras over
               (Dedekind domains and euclidean domains
                and noetherian rings
                and infinite enumerated sets and metric spaces)
              and infinite sets)

            sage: L = LazyLaurentSeriesRing(GF(5), 't')
            sage: TestSuite(L).run()

            sage: L = LazyLaurentSeriesRing(GF(5)['x'], 't')
            sage: TestSuite(L).run()

            sage: L = LazyLaurentSeriesRing(GF(5)['x, y'], 't')
            sage: TestSuite(L).run()

            sage: L = LazyLaurentSeriesRing(Zmod(6), 't')
            sage: TestSuite(L).run(skip=['_test_revert'])
            sage: L.category()
            Category of infinite commutative algebras over
             (finite commutative rings and subquotients of monoids
              and quotients of semigroups and finite enumerated sets)

            sage: E.<x,y> = ExteriorAlgebra(QQ)                                         # needs sage.modules
            sage: L = LazyLaurentSeriesRing(E, 't')     # not tested                    # needs sage.modules

            sage: LazyLaurentSeriesRing.options._reset()  # reset the options
        """
        self._sparse = sparse
        if len(names) != 1:
            raise ValueError("only univariate lazy Laurent series are implemented")
        self._arity = 1
        self._minimal_valuation = None
        self._laurent_poly_ring = LaurentPolynomialRing(base_ring, names, sparse=sparse)
        self._internal_poly_ring = self._laurent_poly_ring

        category = Algebras(base_ring.category())
        if base_ring in Fields():
            category &= CompleteDiscreteValuationFields()
        elif base_ring in IntegralDomains():
            category &= IntegralDomains()
        elif "Commutative" in base_ring.category().axioms():
            category = category.Commutative()

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
        coeff_stream = Stream_exact([R.one()], constant=R.zero(), order=1)
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
    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.gens()
            (z,)
            sage: 1/(1 - z)
            1 + z + z^2 + O(z^3)
        """
        return tuple([self.gen(n) for n in range(self.ngens())])

    def _an_element_(self):
        """
        Return a Laurent series in ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.an_element()
            z^-2 + z^3 + z^4 + z^5 + O(z^6)
        """
        return self(self._laurent_poly_ring.an_element(),
                    valuation=-2,
                    degree=3,
                    constant=self.base_ring().an_element())

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.some_elements()[:7]
            [0, 1, z,
             -3*z^-4 + z^-3 - 12*z^-2 - 2*z^-1 - 10 - 8*z + z^2 + z^3,
             z^-2 + z^3 + z^4 + z^5 + O(z^6),
             -2*z^-3 - 2*z^-2 + 4*z^-1 + 11 - z - 34*z^2 - 31*z^3 + O(z^4),
             4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + O(z^5)]

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L.some_elements()[:7]
            [0, 1, z,
             z^-4 + z^-3 + z^2 + z^3,
             z^-2,
             1 + z + z^3 + z^4 + z^6 + O(z^7),
             z^-1 + z + z^3 + O(z^5)]

            sage: L = LazyLaurentSeriesRing(GF(3), 'z')
            sage: L.some_elements()[:7]
            [0, 1, z,
             z^-3 + z^-1 + 2 + z + z^2 + z^3,
             z^-2,
             z^-3 + z^-2 + z^-1 + 2 + 2*z + 2*z^2 + O(z^3),
             z^-2 + z^-1 + z + z^2 + z^4 + O(z^5)]
        """
        z = self.gen()
        elts = [self.zero(), self.one(), z, (z-3)*(z**-2+2+z)**2, self.an_element(),
                (1 - 2*z**-3)/(1 - z + 3*z**2),
                self(lambda n: n**2, valuation=-2),
                self(lambda n: n**2, valuation=1),
                self([3, 2, 1], valuation=1, constant=1)]
        return elts

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
            coeff_stream = Stream_exact(coefficient, order=valuation,
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

    def uniformizer(self):
        """
        Return a uniformizer of ``self``..

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: L.uniformizer()
            z
        """
        R = self.base_ring()
        if R not in Fields():
            raise TypeError("the base ring is not a field")
        return self.gen()

    def residue_field(self):
        """
        Return the residue field of the ring of integers of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(QQ, 'z')
            sage: L.residue_field()
            Rational Field
        """
        R = self.base_ring()
        if R not in Fields():
            raise TypeError("the base ring is not a field")
        return R

    def taylor(self, f):
        r"""
        Return the Taylor expansion around `0` of the function ``f``.

        INPUT:

        - ``f`` -- a function such that one of the following works:

          * the substitution `f(z)`, where `z` is a generator of ``self``
          * `f` is a function of a single variable with no poles at `0`
            and has a ``derivative`` method

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: x = SR.var('x')
            sage: f(x) = (1 + x) / (1 - x^2)
            sage: L.taylor(f)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)

        For inputs as symbolic functions/expressions, the function must
        not have any poles at `0`::

            sage: f(x) = (1 + x^2) / sin(x^2)
            sage: L.taylor(f)
            <repr(...) failed: ValueError: power::eval(): division by zero>
            sage: def g(a): return (1 + a^2) / sin(a^2)
            sage: L.taylor(g)
            z^-2 + 1 + 1/6*z^2 + 1/6*z^4 + O(z^5)
        """
        try:
            return f(self.gen())
        except (ValueError, TypeError):
            pass
        stream = Stream_taylor(f, self.is_sparse())
        return self.element_class(self, stream)

    # === special functions ===

    def q_pochhammer(self, q=None):
        r"""
        Return the infinite ``q``-Pochhammer symbol `(a; q)_{\infty}`,
        where `a` is the variable of ``self``.

        This is also one version of the quantum dilogarithm or
        the `q`-Exponential function.

        INPUT:

        - ``q`` -- (default: `q \in \QQ(q)`) the parameter `q`

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: L.<z> = LazyLaurentSeriesRing(q.parent())
            sage: qpoch = L.q_pochhammer(q)
            sage: qpoch
            1
             + (-1/(-q + 1))*z
             + (q/(q^3 - q^2 - q + 1))*z^2
             + (-q^3/(-q^6 + q^5 + q^4 - q^2 - q + 1))*z^3
             + (q^6/(q^10 - q^9 - q^8 + 2*q^5 - q^2 - q + 1))*z^4
             + (-q^10/(-q^15 + q^14 + q^13 - q^10 - q^9 - q^8 + q^7 + q^6 + q^5 - q^2 - q + 1))*z^5
             + (q^15/(q^21 - q^20 - q^19 + q^16 + 2*q^14 - q^12 - q^11 - q^10 - q^9 + 2*q^7 + q^5 - q^2 - q + 1))*z^6
             + O(z^7)

        We show that `(z; q)_n = \frac{(z; q)_{\infty}}{(q^n z; q)_{\infty}}`::

            sage: qpoch / qpoch(q*z)
            1 - z + O(z^7)
            sage: qpoch / qpoch(q^2*z)
            1 + (-q - 1)*z + q*z^2 + O(z^7)
            sage: qpoch / qpoch(q^3*z)
            1 + (-q^2 - q - 1)*z + (q^3 + q^2 + q)*z^2 - q^3*z^3 + O(z^7)
            sage: qpoch / qpoch(q^4*z)
            1 + (-q^3 - q^2 - q - 1)*z + (q^5 + q^4 + 2*q^3 + q^2 + q)*z^2
             + (-q^6 - q^5 - q^4 - q^3)*z^3 + q^6*z^4 + O(z^7)

        We can also construct part of Euler's function::

            sage: M.<a> = LazyLaurentSeriesRing(QQ)
            sage: phi = sum(qpoch[i](q=a)*a^i for i in range(10))
            sage: phi[:20] == M.euler()[:20]
            True

        TESTS::

            sage: R = ZZ['q'].fraction_field()
            sage: q = R.gen()
            sage: L.<z> = LazyLaurentSeriesRing(LazyDirichletSeriesRing(R, "s"))
            sage: z.q_pochhammer(q)                                                     # needs sage.symbolic
            1 + ((1/(q-1)))*z + ((q/(q^3-q^2-q+1)))*z^2 + ... + O(z^7)

        REFERENCES:

        - :wikipedia:`Q-Pochhammer_symbol`
        - :wikipedia:`Quantum_dilogarithm`
        - :wikipedia:`Q-exponential`
        """
        if q is None:
            q = ZZ['q'].fraction_field().gen()
        if q not in self.base_ring():
            raise ValueError("q must be in the base ring")
        from sage.arith.misc import binomial
        qP = q.parent()
        one = qP.one()

        def coeff(n):
            return (-1)**n * q**binomial(n, 2) / qP.prod(one - q**i for i in range(1, n+1))
        return self(coefficients=coeff, valuation=0)

    def euler(self):
        r"""
        Return the Euler function as an element of ``self``.

        The *Euler function* is defined as

        .. MATH::

            \phi(z) = (z; z)_{\infty}
            = \sum_{n=0}^{\infty} (-1)^n q^{(3n^2-n)/2}.

        EXAMPLES::

            sage: L.<q> = LazyLaurentSeriesRing(ZZ)
            sage: phi = q.euler()
            sage: phi
            1 - q - q^2 + q^5 + O(q^7)

        We verify that `1 / phi` gives the generating function
        for all partitions::

            sage: P = 1 / phi; P
            1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + O(q^7)
            sage: P[:20] == [Partitions(n).cardinality() for n in range(20)]            # needs sage.libs.flint
            True

        TESTS::

            sage: L.<q> = LazyLaurentSeriesRing(LazyDirichletSeriesRing(QQ, "s"))
            sage: q.euler()                                                             # needs sage.symbolic
            1 - q - q^2 + q^5 + O(q^7)

        REFERENCES:

        - :wikipedia:`Euler_function`
        """
        def coeff(n):
            k = ZZ(24 * n + 1)
            m, rem = k.sqrtrem()
            if rem:
                return ZZ.zero()
            return (-1) ** ((m + 1) // 6)
        return self(coefficients=coeff, valuation=0)

######################################################################


class LazyPowerSeriesRing(LazySeriesRing):
    """
    The ring of (possibly multivariate) lazy Taylor series.

    INPUT:

    - ``base_ring`` -- base ring of this Taylor series ring
    - ``names`` -- name(s) of the generator of this Taylor series ring
    - ``sparse`` -- boolean (default: ``True``); whether this series is sparse or not

    EXAMPLES::

        sage: LazyPowerSeriesRing(ZZ, 't')
        Lazy Taylor Series Ring in t over Integer Ring

        sage: L.<x, y> = LazyPowerSeriesRing(QQ); L
        Multivariate Lazy Taylor Series Ring in x, y over Rational Field
    """
    Element = LazyPowerSeries

    # Follow the "generic" normalization
    __classcall_private__ = LazySeriesRing.__classcall_private__

    def __init__(self, base_ring, names, sparse=True, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: LazyPowerSeriesRing.options.halting_precision(12)

            sage: L = LazyPowerSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip='_test_fraction_field')
            sage: L = LazyPowerSeriesRing(ZZ, 's, t')
            sage: TestSuite(L).run(skip='_test_fraction_field')

            sage: L = LazyPowerSeriesRing(QQ, 't')
            sage: TestSuite(L).run(skip='_test_fraction_field')
            sage: L = LazyPowerSeriesRing(QQ, 's, t')
            sage: TestSuite(L).run(skip='_test_fraction_field')

            sage: L = LazyPowerSeriesRing(GF(5), 't')
            sage: TestSuite(L).run()

            sage: L = LazyPowerSeriesRing(GF(5), 's, t')
            sage: TestSuite(L).run(skip=['_test_fraction_field'])

            sage: L = LazyPowerSeriesRing(Zmod(6), 't')
            sage: TestSuite(L).run(skip=['_test_revert'])
            sage: L = LazyPowerSeriesRing(Zmod(6), 's, t')
            sage: TestSuite(L).run(skip=['_test_revert'])

            sage: L = LazyPowerSeriesRing(QQ['q'], 't')
            sage: TestSuite(L).run(skip='_test_fraction_field')
            sage: L = LazyPowerSeriesRing(QQ['q'], 's, t')
            sage: TestSuite(L).run(skip='_test_fraction_field')  # long time

            sage: L = LazyPowerSeriesRing(ZZ['q'], 't')
            sage: TestSuite(L).run(skip='_test_fraction_field')
            sage: L = LazyPowerSeriesRing(ZZ['q'], 's, t')
            sage: TestSuite(L).run(skip='_test_fraction_field')  # long time

            sage: LazyPowerSeriesRing.options._reset()  # reset the options

        Check that :issue:`34470` is fixed::

            sage: L.<t> = LazyPowerSeriesRing(QQ)
            sage: L in CompleteDiscreteValuationRings
            True
            sage: L.uniformizer()
            t
            sage: lcm(1/(1 - t^2) - 1, t)
            t^2

            sage: L.<t> = LazyPowerSeriesRing(ZZ)
            sage: L in PrincipalIdealDomains
            False

        The ideal generated by `s` and `t` is not principal::

            sage: L = LazyPowerSeriesRing(QQ, 's, t')
            sage: L in PrincipalIdealDomains
            False
        """
        self._sparse = sparse
        self._minimal_valuation = 0
        self._arity = len(names)
        if self._arity == 1:
            self._laurent_poly_ring = PolynomialRing(base_ring, names, sparse=sparse)
            self._internal_poly_ring = self._laurent_poly_ring
        else:
            self._laurent_poly_ring = PolynomialRing(base_ring, names)
            self._internal_poly_ring = PolynomialRing(self._laurent_poly_ring, "DUMMY_VARIABLE", sparse=sparse)
        category = Algebras(base_ring.category())
        mixin_gcd = False
        if self._arity == 1:
            if base_ring in Fields():
                category &= CompleteDiscreteValuationRings()
                mixin_gcd = True
        elif base_ring in Fields():
            category &= UniqueFactorizationDomains()
            mixin_gcd = True
        if base_ring in IntegralDomains():
            category &= IntegralDomains()
        elif base_ring in Rings().Commutative():
            category = category.Commutative()

        if mixin_gcd:
            from sage.structure.dynamic_class import dynamic_class
            self.Element = dynamic_class(
                f"{self.Element.__name__}_gcd",
                (self.Element, LazyPowerSeries_gcd_mixin),
                doccls=self.Element)

        if base_ring.is_zero():
            category = category.Finite()
        else:
            category = category.Infinite()
        Parent.__init__(self, base=base_ring, names=names,
                        category=category)

    def _repr_(self):
        """
        String representation of this Taylor series ring.

        EXAMPLES::

            sage: LazyPowerSeriesRing(GF(2), 'z')
            Lazy Taylor Series Ring in z over Finite Field of size 2
        """
        BR = self.base_ring()
        if len(self.variable_names()) == 1:
            return "Lazy Taylor Series Ring in {} over {}".format(self.variable_name(), BR)
        generators_rep = ", ".join(self.variable_names())
        return "Multivariate Lazy Taylor Series Ring in {} over {}".format(generators_rep, BR)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(GF(2), 'z')
            sage: latex(L)
            \Bold{F}_{2} [\![z]\!]
        """
        from sage.misc.latex import latex
        generators_rep = ", ".join(self.variable_names())
        return latex(self.base_ring()) + r"[\![{}]\!]".format(generators_rep)

    def _monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L._monomial(2, 3)
            2*z^3
        """
        m = len(self.variable_names())
        L = self._laurent_poly_ring
        if m == 1:
            return L(c) * L.gen() ** n
        return L(c)

    @cached_method
    def gen(self, n=0):
        """
        Return the ``n``-th generator of ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.gen()
            z
            sage: L.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        m = len(self.variable_names())
        if n > m:
            if m == 1:
                raise IndexError("there is only one generator")
            raise IndexError("there are only %s generators" % m)

        R = self._laurent_poly_ring
        BR = self.base_ring()
        if len(self.variable_names()) == 1:
            coeff_stream = Stream_exact([BR.one()], constant=BR.zero(), order=1)
        else:
            coeff_stream = Stream_exact([R.gen(n)], constant=BR.zero(), order=1)
        return self.element_class(self, coeff_stream)

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyPowerSeriesRing(ZZ)
            sage: L.ngens()
            1
        """
        return len(self.variable_names())

    @cached_method
    def gens(self) -> tuple:
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(ZZ, 'x,y')
            sage: L.gens()
            (x, y)
        """
        return tuple([self.gen(n) for n in range(self.ngens())])

    def _element_constructor_(self, x=None, valuation=None, constant=None, degree=None, coefficients=None, check=True):
        """
        Construct a Taylor series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a Taylor series
        - ``valuation`` -- integer (optional); a lower bound for the valuation
          of the series
        - ``constant`` -- (optional) the eventual constant of the series
        - ``degree`` -- (optional) the degree when the series is ``constant``
        - ``check`` -- (optional) check that coefficients are homogeneous of
          the correct degree when they are retrieved

        .. WARNING::

            The behaviour of ``LazyPowerSeries(c)`` for a list ``c``
            with nonzero last element `e` changed with
            :issue:`32367`.  To obtain the old behaviour, use
            ``LazyPowerSeries(c, constant=e)``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L(lambda i: i, 5, 1, 10)
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)
            sage: L(lambda i: i, 5, (1, 10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: X = L(constant=5, degree=2); X
            5*z^2 + 5*z^3 + 5*z^4 + O(z^5)
            sage: X.valuation()
            2

            sage: e = L(lambda n: n + 1); e
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + O(z^7)
            sage: f = e^-1; f
            1 - 2*z + z^2 + O(z^7)
            sage: f.coefficient(10)
            0
            sage: f[20]
            0

            sage: L(valuation=2, constant=1)
            z^2 + z^3 + z^4 + O(z^5)
            sage: L(constant=1)
            1 + z + z^2 + O(z^3)

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], 1); f
            z + 2*z^2 + 3*z^3 + 4*z^4

            sage: g = L([1,3,5,7,9], 5, -1); g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + O(z^13)

        Additionally, ``x`` can be a polynomial::

            sage: P.<x> = QQ[]
            sage: p = x + 3*x^2 + x^5
            sage: L.<x> = LazyPowerSeriesRing(ZZ)
            sage: L(p)
            x + 3*x^2 + x^5

            sage: L(p, valuation=0)
            1 + 3*x + x^4

            sage: P.<x, y> = QQ[]
            sage: p = x + y^2 + x*y
            sage: L.<x,y> = LazyPowerSeriesRing(ZZ)
            sage: L(p)
            x + (x*y+y^2)

        Finally ``x`` can be in the corresponding fraction field::

            sage: R.<a,b,c> = PolynomialRing(ZZ)
            sage: L = LazyPowerSeriesRing(ZZ, 'a,b,c')
            sage: aa, bb, cc = L.gens()
            sage: f = (1 + a + b) / (1 + a*b + c^3); f
            (a + b + 1)/(c^3 + a*b + 1)
            sage: f.parent()
            Fraction Field of Multivariate Polynomial Ring in a, b, c over Integer Ring
            sage: L(f)                                                                  # needs sage.libs.singular
            1 + (a+b) - a*b - (a^2*b+a*b^2+c^3) + (a^2*b^2-a*c^3-b*c^3)
            + (a^3*b^2+a^2*b^3+2*a*b*c^3) - (a^3*b^3-2*a^2*b*c^3-2*a*b^2*c^3-c^6)
            + O(a,b,c)^7
            sage: L(f) == (1 + aa + bb) / (1 + aa*bb + cc^3)                            # needs sage.libs.singular
            True

        TESTS::

            sage: L.<x,y> = LazyPowerSeriesRing(ZZ)
            sage: L(constant=1)
            Traceback (most recent call last):
            ...
            ValueError: constant must be zero for multivariate Taylor series

            sage: L(lambda n: 0)
            O(x,y)^7

            sage: L(lambda n: n)[3];
            Traceback (most recent call last):
            ...
            ValueError: coefficient 3 at degree 3 is not a homogeneous polynomial

            sage: L([1, 2, 3]);
            Traceback (most recent call last):
            ...
            ValueError: unable to convert [1, 2, 3] into a lazy Taylor series

            sage: L(lambda n: n, degree=3);
            Traceback (most recent call last):
            ...
            ValueError: coefficients must be homogeneous polynomials of the correct degree
        """
        if valuation is not None:
            if valuation < 0:
                raise ValueError("the valuation of a Taylor series must be nonnegative")
            # TODO: the following is nonsense, think of an iterator
            if self._arity > 1:
                raise ValueError("valuation must not be specified for multivariate Taylor series")
        if self._arity > 1:
            valuation = 0

        R = self._laurent_poly_ring
        BR = self.base_ring()
        if x is None:
            assert degree is None
            coeff_stream = Stream_uninitialized(valuation)
            return self.element_class(self, coeff_stream)

        try:
            # Try to build stuff using the polynomial ring constructor
            x = R(x)
        except (TypeError, ValueError):
            pass
        if isinstance(constant, (tuple, list)):
            constant, degree = constant
        if constant is not None:
            if self._arity > 1 and constant:
                raise ValueError("constant must be zero for multivariate Taylor series")
            constant = BR(constant)

        if x in R:
            if not x and not constant:
                coeff_stream = Stream_zero()
            else:
                if not x:
                    coeff_stream = Stream_exact([],
                                                order=valuation,
                                                degree=degree,
                                                constant=constant)
                    return self.element_class(self, coeff_stream)

                if self._arity == 1:
                    v = x.valuation()
                    d = x.degree()
                    p_list = [x[i] for i in range(v, d + 1)]
                    if valuation is not None:
                        v = valuation
                else:
                    p_dict = x.homogeneous_components()
                    v = min(p_dict.keys())
                    d = max(p_dict.keys())
                    p_list = [p_dict.get(i, 0) for i in range(v, d + 1)]

                coeff_stream = Stream_exact(p_list,
                                            order=v,
                                            constant=constant,
                                            degree=degree)
            return self.element_class(self, coeff_stream)

        if isinstance(x, LazyPowerSeries):
            stream = x._coeff_stream
            if isinstance(stream, Stream_exact):
                if self._arity == 1:
                    BR = self.base_ring()
                else:
                    BR = self._laurent_poly_ring
                coeffs = [BR(val) for val in stream._initial_coefficients]
                valuation = stream._approximate_order
                for i, c in enumerate(coeffs):
                    if c:
                        valuation += i
                        coeffs = coeffs[i:]
                        break
                else:
                    valuation += len(coeffs)
                    coeffs = []
                return self(coeffs,
                            degree=stream._degree,
                            constant=self.base_ring()(stream._constant),
                            valuation=valuation)
            return self.element_class(self, stream)

        # Check if we can realize the input as a rational function
        try:
            FF = self._laurent_poly_ring.fraction_field()
            x = FF(x)
        except (TypeError, ValueError, AttributeError):
            pass
        else:
            return self(x.numerator()) / self(x.denominator())

        if callable(x) or isinstance(x, (GeneratorType, map, filter)):
            if valuation is None:
                valuation = 0
            if degree is not None:
                if constant is None:
                    constant = ZZ.zero()
                if callable(x):
                    p = [x(i) for i in range(valuation, degree)]
                else:
                    p = [c for c, _ in zip(_skip_leading_zeros(x), range(valuation, degree))]
                if self._arity == 1:
                    p = [BR(c) for c in p]
                else:
                    p = [R(c) for c in p]
                    if not all(e.is_homogeneous() and e.degree() == i
                               for i, e in enumerate(p, valuation)):
                        raise ValueError("coefficients must be homogeneous polynomials of the correct degree")
                coeff_stream = Stream_exact(p,
                                            order=valuation,
                                            constant=constant,
                                            degree=degree)
                return self.element_class(self, coeff_stream)
            if check and self._arity > 1:
                if callable(x):
                    def y(n):
                        e = R(x(n))
                        if not e or e.is_homogeneous() and e.degree() == n:
                            return e
                        raise ValueError("coefficient %s at degree %s is not a homogeneous polynomial" % (e, n))
                    coeff_stream = Stream_function(y, self._sparse, valuation)
                else:
                    coeff_stream = Stream_iterator(map(R, _skip_leading_zeros(x)), valuation)
            else:
                if callable(x):
                    coeff_stream = Stream_function(lambda i: BR(x(i)), self._sparse, valuation)
                else:
                    coeff_stream = Stream_iterator(map(BR, _skip_leading_zeros(x)), valuation)
            return self.element_class(self, coeff_stream)
        raise ValueError(f"unable to convert {x} into a lazy Taylor series")

    def _an_element_(self):
        """
        Return a Taylor series in ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.an_element()
            z + z^2 + z^3 + O(z^4)

            sage: L = LazyPowerSeriesRing(ZZ, 'x, y')
            sage: L.an_element()
            x
        """
        if self._arity == 1:
            return self(self._laurent_poly_ring.an_element(),
                        constant=self.base_ring().an_element())
        return self(self._laurent_poly_ring.an_element())

    def uniformizer(self):
        """
        Return a uniformizer of ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ, 'x')
            sage: L.uniformizer()
            x
        """
        R = self.base_ring()
        if R not in Fields():
            raise TypeError("the base ring is not a field")
        if self._arity != 1:
            raise TypeError("the arity must be one")
        return self.gen()

    def residue_field(self):
        """
        Return the residue field of the ring of integers of ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ, 'x')
            sage: L.residue_field()
            Rational Field
        """
        R = self.base_ring()
        if R not in Fields():
            raise TypeError("the base ring is not a field")
        if self._arity != 1:
            raise TypeError("the arity must be one")
        return R

    def fraction_field(self):
        """
        Return the fraction field of ``self``.

        If this is with a single variable over a field, then the fraction
        field is the field of (lazy) formal Laurent series.

        .. TODO::

            Implement other fraction fields.

        EXAMPLES::

            sage: L.<x> = LazyPowerSeriesRing(QQ)
            sage: L.fraction_field()
            Lazy Laurent Series Ring in x over Rational Field
        """
        if self not in IntegralDomains():
            raise TypeError("must be an integral domain")
        R = self.base_ring()
        if self._arity == 1 and R in Fields():
            return LazyLaurentSeriesRing(R, names=self.variable_names())
        raise NotImplementedError("the fraction field is not yet implemented")

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(ZZ, 'z')
            sage: L.some_elements()[:6]
            [0, 1, z + z^2 + z^3 + O(z^4),
             -12 - 8*z + z^2 + z^3,
             1 + z - 2*z^2 - 7*z^3 - z^4 + 20*z^5 + 23*z^6 + O(z^7),
             z + 4*z^2 + 9*z^3 + 16*z^4 + 25*z^5 + 36*z^6 + O(z^7)]

            sage: L = LazyPowerSeriesRing(GF(3)["q"], 'z')
            sage: L.some_elements()[:6]
            [0, 1, z + q*z^2 + q*z^3 + q*z^4 + O(z^5),
             z + z^2 + z^3,
             1 + z + z^2 + 2*z^3 + 2*z^4 + 2*z^5 + O(z^6),
             z + z^2 + z^4 + z^5 + O(z^7)]

            sage: L = LazyPowerSeriesRing(GF(3), 'q, t')
            sage: L.some_elements()[:6]
            [0, 1, q,
             q + q^2 + q^3,
             1 + q + q^2 - q^3 - q^4 - q^5 - q^6 + O(q,t)^7,
             1 + (q+t) + (q^2-q*t+t^2) + (q^3+t^3) + (q^4+q^3*t+q*t^3+t^4)
             + (q^5-q^4*t+q^3*t^2+q^2*t^3-q*t^4+t^5) + (q^6-q^3*t^3+t^6) + O(q,t)^7]
        """
        z = self.gen(0)
        elts = [self.zero(), self.one(), self.an_element()]
        if self._arity == 1:
            elts.extend([(z-3)*(2+z)**2, (1 - 2*z**3)/(1 - z + 3*z**2), self(lambda n: n**2)])
        else:
            PR = self._laurent_poly_ring
            sum_gens = PR.sum(PR.gens())
            elts.extend([(z-3)*(2+z)**2, (1 - 2*z**3)/(1 - z + 3*z**2), self(lambda n: sum_gens**n)])
        return elts

    def taylor(self, f):
        r"""
        Return the Taylor expansion around `0` of the function ``f``.

        INPUT:

        - ``f`` -- a function such that one of the following works:

          * the substitution `f(z_1, \ldots, z_n)`, where `(z_1, \ldots, z_n)`
            are the generators of ``self``
          * `f` is a function with no poles at `0` and has a ``derivative``
            method

        .. WARNING::

            For inputs as symbolic functions/expressions, this does not check
            that the function does not have poles at `0`.

        EXAMPLES::

            sage: L.<z> = LazyPowerSeriesRing(QQ)
            sage: x = SR.var('x')
            sage: f(x) = (1 + x) / (1 - x^3)
            sage: L.taylor(f)
            1 + z + z^3 + z^4 + z^6 + O(z^7)
            sage: (1 + z) / (1 - z^3)
            1 + z + z^3 + z^4 + z^6 + O(z^7)
            sage: f(x) = cos(x + pi/2)
            sage: L.taylor(f)
            -z + 1/6*z^3 - 1/120*z^5 + O(z^7)

        For inputs as symbolic functions/expressions, the function must
        not have any poles at `0`::

            sage: L.<z> = LazyPowerSeriesRing(QQ, sparse=True)
            sage: f = 1 / sin(x)
            sage: L.taylor(f)
            <repr(...) failed: ValueError: power::eval(): division by zero>

        Different multivariate inputs::

            sage: L.<a,b> = LazyPowerSeriesRing(QQ)
            sage: def f(x, y): return (1 + x) / (1 + y)
            sage: L.taylor(f)
            1 + (a-b) - (a*b-b^2) + (a*b^2-b^3) - (a*b^3-b^4) + (a*b^4-b^5) - (a*b^5-b^6) + O(a,b)^7
            sage: g(w, z) = (1 + w) / (1 + z)
            sage: L.taylor(g)
            1 + (a-b) - (a*b-b^2) + (a*b^2-b^3) - (a*b^3-b^4) + (a*b^4-b^5) - (a*b^5-b^6) + O(a,b)^7
            sage: y = SR.var('y')
            sage: h = (1 + x) / (1 + y)
            sage: L.taylor(h)
            1 + (a-b) - (a*b-b^2) + (a*b^2-b^3) - (a*b^3-b^4) + (a*b^4-b^5) - (a*b^5-b^6) + O(a,b)^7
        """
        try:
            return f(*self.gens())
        except (ValueError, TypeError):
            pass

        if self._arity != 1:
            R = self._laurent_poly_ring
            BR = R.base_ring()
            args = f.arguments()
            subs = {str(va): ZZ.zero() for va in args}
            ell = len(subs)
            from sage.combinat.integer_vector import integer_vectors_nk_fast_iter
            from sage.arith.misc import factorial

            def taylor_expand(deg):
                if deg == 0:
                    return BR(f(**subs))
                return R.sum(BR(f.diff(*sum(([g] * e for g, e in zip(args, al)), []))(**subs)
                                / ZZ.prod(factorial(a) for a in al))
                             * R.monomial(*al) for al in integer_vectors_nk_fast_iter(deg, ell))

            coeff_stream = Stream_function(taylor_expand, self._sparse, self._minimal_valuation)
        else:
            coeff_stream = Stream_taylor(f, self._sparse)
        return self.element_class(self, coeff_stream)


######################################################################


class LazyCompletionGradedAlgebra(LazySeriesRing):
    r"""
    The completion of a graded algebra consisting of formal series.

    For a graded algebra `A`, we can form a completion of `A` consisting of
    all formal series of `A` such that each homogeneous component is
    a finite linear combination of basis elements of `A`.

    INPUT:

    - ``basis`` -- a graded algebra
    - ``names`` -- name(s) of the alphabets
    - ``sparse`` -- boolean (default: ``True``); whether we use a sparse or
      a dense representation

    EXAMPLES::

        sage: # needs sage.modules
        sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
        sage: S = NCSF.Complete()
        sage: L = S.formal_series_ring(); L
        Lazy completion of Non-Commutative Symmetric Functions
         over the Rational Field in the Complete basis
        sage: f = 1 / (1 - L(S[1])); f
        S[] + S[1] + (S[1,1]) + (S[1,1,1]) + (S[1,1,1,1]) + (S[1,1,1,1,1])
         + (S[1,1,1,1,1,1]) + O^7
        sage: g = 1 / (1 - L(S[2])); g
        S[] + S[2] + (S[2,2]) + (S[2,2,2]) + O^7
        sage: f * g
        S[] + S[1] + (S[1,1]+S[2]) + (S[1,1,1]+S[1,2])
         + (S[1,1,1,1]+S[1,1,2]+S[2,2]) + (S[1,1,1,1,1]+S[1,1,1,2]+S[1,2,2])
         + (S[1,1,1,1,1,1]+S[1,1,1,1,2]+S[1,1,2,2]+S[2,2,2]) + O^7
        sage: g * f
        S[] + S[1] + (S[1,1]+S[2]) + (S[1,1,1]+S[2,1])
         + (S[1,1,1,1]+S[2,1,1]+S[2,2]) + (S[1,1,1,1,1]+S[2,1,1,1]+S[2,2,1])
         + (S[1,1,1,1,1,1]+S[2,1,1,1,1]+S[2,2,1,1]+S[2,2,2]) + O^7
        sage: f * g - g * f
        (S[1,2]-S[2,1]) + (S[1,1,2]-S[2,1,1])
         + (S[1,1,1,2]+S[1,2,2]-S[2,1,1,1]-S[2,2,1])
         + (S[1,1,1,1,2]+S[1,1,2,2]-S[2,1,1,1,1]-S[2,2,1,1]) + O^7
    """
    Element = LazyCompletionGradedAlgebraElement

    def __init__(self, basis, sparse=True, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: LazySymmetricFunctions.options.halting_precision(6)

            sage: # needs sage.modules
            sage: s = SymmetricFunctions(QQ).s()
            sage: L = LazySymmetricFunctions(s)
            sage: TestSuite(L).run()                                                    # needs lrcalc_python
            sage: p = SymmetricFunctions(GF(5)).p()
            sage: L = LazySymmetricFunctions(p)
            sage: TestSuite(L).run()

        Reversion will only work when the base ring is a field::

            sage: # needs sage.modules
            sage: s = SymmetricFunctions(ZZ).s()
            sage: L = LazySymmetricFunctions(s)
            sage: TestSuite(L).run(skip=['_test_revert'])                               # needs lrcalc_python
            sage: s = SymmetricFunctions(QQ["q"]).s()
            sage: L = LazySymmetricFunctions(s)
            sage: TestSuite(L).run(skip=['_test_revert'])                               # needs lrcalc_python

        Options are remembered across doctests::

            sage: LazySymmetricFunctions.options._reset()

        Check that :issue:`34470` is fixed.  The ideal generated by
        `p[1]` and `p[2]` is not principal::

            sage: p = SymmetricFunctions(QQ).p()                                        # needs sage.modules
            sage: L = LazySymmetricFunctions(s)                                         # needs sage.modules
            sage: L in PrincipalIdealDomains                                            # needs sage.modules
            False

        Check that a basis which is not graded is not enough::

            sage: ht = SymmetricFunctions(ZZ).ht()                                      # needs sage.modules
            sage: L = LazySymmetricFunctions(ht)                                        # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: basis should be in GradedAlgebrasWithBasis

        Check that :issue:`37625` is fixed::

            sage: R = algebras.Free(QQ, ('a', 'b'), degrees=(1, 2))
            sage: L = R.completion()
            sage: a, b = R.gens()
            sage: (L(a) + L(b))^2
            a^2 + (a*b+b*a) + b^2
        """
        base_ring = basis.base_ring()
        self._minimal_valuation = 0
        if basis in Algebras.TensorProducts:
            self._arity = len(basis._sets)
        else:
            if basis not in GradedAlgebrasWithBasis:
                raise ValueError("basis should be in GradedAlgebrasWithBasis")
            self._arity = 1
        category = Algebras(basis.category())
        if basis in IntegralDomains():
            category &= IntegralDomains()
        elif basis in Rings().Commutative():
            category = category.Commutative()

        if base_ring.is_zero():
            category = category.Finite()
        else:
            category = category.Infinite()
        Parent.__init__(self, base=base_ring, category=category)
        self._sparse = sparse
        self._laurent_poly_ring = basis
        self._internal_poly_ring = PolynomialRing(self._laurent_poly_ring, "DUMMY_VARIABLE", sparse=sparse)

    def _repr_(self):
        """
        String representation of the lazy symmetric functions ring.

        EXAMPLES::

            sage: s = SymmetricFunctions(GF(2)).s()                                     # needs sage.modules
            sage: LazySymmetricFunctions(s)                                             # needs sage.modules
            Lazy completion of Symmetric Functions over Finite Field of size 2 in the Schur basis
        """
        return "Lazy completion of {}".format(self._laurent_poly_ring)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: s = SymmetricFunctions(GF(2)).s()                                     # needs sage.modules
            sage: L = LazySymmetricFunctions(s)                                         # needs sage.modules
            sage: latex(L)                                                              # needs sage.modules
            \text{\texttt{Symmetric{ }Functions{ }over{ }Finite{ }Field{ }of{ }size{ }2{ }in{ }the{ }Schur{ }basis}}
        """
        from sage.misc.latex import latex
        return latex(self._laurent_poly_ring)

    def _monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: m = SymmetricFunctions(ZZ).m()
            sage: s = SymmetricFunctions(ZZ).s()
            sage: L = LazySymmetricFunctions(m)
            sage: L._monomial(s[2,1], 3)
            2*m[1, 1, 1] + m[2, 1]
        """
        L = self._laurent_poly_ring
        return L(c)

    def _element_constructor_(self, x=None, valuation=None, degree=None, constant=None, check=True):
        r"""
        Construct a lazy element in ``self`` from ``x``.

        INPUT:

        - ``x`` -- data used to the define a lazy element
        - ``valuation`` -- integer (optional); a lower bound for the valuation
          of the series
        - ``degree`` -- (optional) the degree when the lazy element
          has finite support
        - ``check`` -- (optional) check that coefficients are homogeneous of
          the correct degree when they are retrieved

        EXAMPLES::

            sage: # needs sage.modules
            sage: m = SymmetricFunctions(GF(2)).m()
            sage: L = LazySymmetricFunctions(m)
            sage: L(2)
            0
            sage: L(3)
            m[]

            sage: # needs sage.modules
            sage: m = SymmetricFunctions(ZZ).m()
            sage: L = LazySymmetricFunctions(m)
            sage: f = L(lambda i: m([i]), valuation=5, degree=10); f
            m[5] + m[6] + m[7] + m[8] + m[9]
            sage: f.coefficient(6)
            m[6]
            sage: f[20]
            0

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``::

            sage: f = L([m[1],m[2],m[3]], valuation=1); f                               # needs sage.modules
            m[1] + m[2] + m[3]

        Finally, ``x`` can be a symmetric function::

            sage: # needs sage.modules
            sage: m = SymmetricFunctions(ZZ).m()
            sage: s = SymmetricFunctions(ZZ).s()
            sage: L = LazySymmetricFunctions(m)
            sage: L(s.an_element())
            2*m[] + 2*m[1] + (3*m[1,1]+3*m[2])

        TESTS::

            sage: # needs sage.modules
            sage: e = SymmetricFunctions(ZZ).e()
            sage: h = SymmetricFunctions(ZZ).h()
            sage: L = LazySymmetricFunctions(tensor([h, e]))
            sage: L(lambda n: 0)
            O^7

            sage: # needs sage.modules
            sage: L(lambda n: tensor([h[n], e([])]) + tensor([h([]), e[n]]), degree=3)
            (2*h[]#e[]) + (h[]#e[1]+h[1]#e[]) + (h[]#e[2]+h[2]#e[])
            sage: L(lambda n: n)[3];
            Traceback (most recent call last):
            ...
            ValueError: coefficient 3*h[] # e[] should be an element
            of homogeneous degree 3 but has degree 0
            sage: L([1, 2, 3]);
            Traceback (most recent call last):
            ...
            ValueError: coefficient 2*h[] # e[] should be an element
            of homogeneous degree 1 but has degree 0
            sage: L(lambda n: n, degree=3);
            Traceback (most recent call last):
            ...
            ValueError: coefficient h[] # e[] should be an element
            of homogeneous degree 1 but has degree 0
        """
        if valuation is None:
            valuation = 0
        if valuation < 0:
            raise ValueError("the valuation of a lazy completion element must be nonnegative")

        R = self._laurent_poly_ring
        if x is None:
            assert degree is None
            coeff_stream = Stream_uninitialized(valuation)
            return self.element_class(self, coeff_stream)
        try:
            # Try to build stuff using the polynomial ring constructor
            x = R(x)
        except (TypeError, ValueError, NotImplementedError):
            pass
        if x in R:
            if not x:
                coeff_stream = Stream_zero()
            else:
                p_dict = {}
                if self._arity == 1:
                    for f in x.terms():
                        d = f.degree()
                        p_dict[d] = p_dict.get(d, 0) + f
                else:
                    for f in x.terms():
                        try:
                            d = f.degree()
                        except (TypeError, ValueError, AttributeError):
                            # FIXME: Fallback for symmetric functions in multiple variables
                            d = sum(sum(mu.size() for mu in p) for p in f.support())
                        p_dict[d] = p_dict.get(d, 0) + f
                v = min(p_dict)
                d = max(p_dict)
                p_list = [p_dict.get(i, 0) for i in range(v, d + 1)]

                coeff_stream = Stream_exact(p_list,
                                            order=v,
                                            constant=0,
                                            degree=degree)
            return self.element_class(self, coeff_stream)

        if isinstance(x, self.Element):
            return self.element_class(self, x._coeff_stream)

        if self._arity == 1:
            def check_homogeneous_of_degree(f, d):
                if not f:
                    return
                try:
                    d1 = f.homogeneous_degree()
                    if d1 == d:
                        return
                except ValueError:
                    raise ValueError("coefficient %s should be an element of homogeneous degree %s" % (f, d))
                raise ValueError("coefficient %s should be an element of homogeneous degree %s but has degree %s" % (f, d, d1))
        else:
            def check_homogeneous_of_degree(f, d):
                if not f:
                    return
                for m in f.monomials():
                    try:
                        d1 = m.degree()
                    except AttributeError:
                        # FIXME: Fallback for symmetric functions in multiple variables
                        for t in m.support():
                            d1 = sum(p.size() for p in t)
                            if d1 != d:
                                raise ValueError("coefficient %s should be an element of homogeneous degree %s but has degree %s" % (f, d, d1))
                    except (TypeError, ValueError):
                        raise ValueError("coefficient %s is not homogeneous")
                    if d1 != d:
                        raise ValueError("coefficient %s should be an element of homogeneous degree %s but has degree %s" % (f, d, d1))

        if isinstance(x, (tuple, list)):
            if degree is None:
                degree = valuation + len(x)
            p = [R(e) for e in x]
            for i, e in enumerate(p, valuation):
                check_homogeneous_of_degree(e, i)
            coeff_stream = Stream_exact(p,
                                        order=valuation,
                                        constant=0,
                                        degree=degree)
            return self.element_class(self, coeff_stream)
        if callable(x):
            if degree is not None:
                p = [R(x(i)) for i in range(valuation, degree)]
                for i, e in enumerate(p, valuation):
                    check_homogeneous_of_degree(e, i)
                coeff_stream = Stream_exact(p,
                                            order=valuation,
                                            constant=0,
                                            degree=degree)
                return self.element_class(self, coeff_stream)
            if check:
                def y(n):
                    e = R(x(n))
                    check_homogeneous_of_degree(e, n)
                    return e

                coeff_stream = Stream_function(y, self._sparse, valuation)
            else:
                coeff_stream = Stream_function(x, self._sparse, valuation)
            return self.element_class(self, coeff_stream)
        raise ValueError(f"unable to convert {x} into a lazy completion element")

    def _an_element_(self):
        """
        Return a lazy symmetric function in ``self``.

        EXAMPLES::

            sage: m = SymmetricFunctions(ZZ).m()                                        # needs sage.modules
            sage: L = LazySymmetricFunctions(m)                                         # needs sage.modules
            sage: L.an_element()                                                        # needs sage.modules
            2*m[] + 2*m[1] + 3*m[2]
        """
        return self(self._laurent_poly_ring.an_element())

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: m = SymmetricFunctions(GF(5)).m()                                     # needs sage.modules
            sage: L = LazySymmetricFunctions(m)                                         # needs sage.modules
            sage: L.some_elements()[:5]                                                 # needs sage.modules
            [0, m[], 2*m[] + 2*m[1] + 3*m[2], 2*m[1] + 3*m[2],
             3*m[] + 2*m[1] + (m[1,1]+m[2])
                   + (2*m[1,1,1]+m[3])
                   + (2*m[1,1,1,1]+4*m[2,1,1]+2*m[2,2])
                   + (3*m[2,1,1,1]+3*m[3,1,1]+4*m[3,2]+m[5])
                   + (2*m[2,2,1,1]+m[2,2,2]+2*m[3,2,1]+2*m[3,3]+m[4,1,1]+3*m[4,2]+4*m[5,1]+4*m[6])
                   + O^7]

            sage: # needs sage.modules
            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: S = NCSF.Complete()
            sage: L = S.formal_series_ring()
            sage: L.some_elements()[:4]
            [0, S[], 2*S[] + 2*S[1] + (3*S[1,1]), 2*S[1] + (3*S[1,1])]
        """
        elt = self.an_element()
        elts = [self.zero(), self.one(), elt]
        # an element with no constant term
        elts.append(elt - elt[0])
        # the inverse of an element
        try:
            if elt.is_unit():
                elts.append(~elt)
            else:
                elts.append(~(1 - elt[0] + elt))
        except NotImplementedError:
            pass
        # an element with no constant term and an invertible
        # coefficient of the linear term
        it = iter(self._laurent_poly_ring.basis())
        temp = self.sum(b for _ in range(4) if (b := next(it)).degree())
        if temp:
            elts.append(temp)

        return elts

######################################################################


class LazySymmetricFunctions(LazyCompletionGradedAlgebra):
    """
    The ring of lazy symmetric functions.

    INPUT:

    - ``basis`` -- the ring of symmetric functions
    - ``names`` -- name(s) of the alphabets
    - ``sparse`` -- boolean (default: ``True``); whether we use a sparse or a
      dense representation

    EXAMPLES::

        sage: s = SymmetricFunctions(ZZ).s()                                            # needs sage.modules
        sage: LazySymmetricFunctions(s)                                                 # needs sage.modules
        Lazy completion of Symmetric Functions over Integer Ring in the Schur basis

        sage: m = SymmetricFunctions(ZZ).m()                                            # needs sage.modules
        sage: LazySymmetricFunctions(tensor([s, m]))                                    # needs sage.modules
        Lazy completion of
         Symmetric Functions over Integer Ring in the Schur basis
          # Symmetric Functions over Integer Ring in the monomial basis
    """
    Element = LazySymmetricFunction


######################################################################

class LazyDirichletSeriesRing(LazySeriesRing):
    r"""
    The ring of lazy Dirichlet series.

    INPUT:

    - ``base_ring`` -- base ring of this Dirichlet series ring
    - ``names`` -- name of the generator of this Dirichlet series ring
    - ``sparse`` -- boolean (default: ``True``); whether this series is sparse or not

    Unlike formal univariate Laurent/power series (over a field),
    the ring of formal Dirichlet series is not a
    :wikipedia:`discrete_valuation_ring`.  On the other hand, it
    is a :wikipedia:`local_ring`.  The unique maximal ideal
    consists of all non-invertible series, i.e., series with
    vanishing constant term.

    .. TODO::

        According to the answers in
        https://mathoverflow.net/questions/5522/dirichlet-series-with-integer-coefficients-as-a-ufd,
        (which, in particular, references :arxiv:`math/0105219`)
        the ring of formal Dirichlet series is actually a
        :wikipedia:`Unique_factorization_domain` over `\ZZ`.

    .. NOTE::

        An interesting valuation is described in Emil Daniel
        Schwab; Gheorghe Silberberg *A note on some discrete
        valuation rings of arithmetical functions*, Archivum
        Mathematicum, Vol. 36 (2000), No. 2, 103-109,
        http://dml.cz/dmlcz/107723.  Let `J_k` be the ideal of
        Dirichlet series whose coefficient `f[n]` of `n^s`
        vanishes if `n` has less than `k` prime factors, counting
        multiplicities.  For any Dirichlet series `f`, let `D(f)`
        be the largest integer `k` such that `f` is in `J_k`.
        Then `D` is surjective, `D(f g) = D(f) + D(g)` for
        nonzero `f` and `g`, and `D(f + g) \geq \min(D(f), D(g))`
        provided that `f + g` is nonzero.

        For example, `J_1` are series with no constant term, and
        `J_2` are series such that `f[1]` and `f[p]` for prime
        `p` vanish.

        Since this is a chain of increasing ideals, the ring of
        formal Dirichlet series is not a
        :wikipedia:`Noetherian_ring`.

        Evidently, this valuation cannot be computed for a given
        series.

    EXAMPLES::

        sage: LazyDirichletSeriesRing(ZZ, 't')
        Lazy Dirichlet Series Ring in t over Integer Ring

    The ideal generated by `2^-s` and `3^-s` is not principal::

        sage: L = LazyDirichletSeriesRing(QQ, 's')
        sage: L in PrincipalIdealDomains
        False
    """
    Element = LazyDirichletSeries

    # Follow the "generic" normalization
    __classcall_private__ = LazySeriesRing.__classcall_private__

    @lazy_attribute
    def _laurent_poly_ring(self):
        r"""
        Return the symbolic ring.

        .. TODO::

            It would be good to have something better than the symbolic ring.

        TESTS::

            sage: L = LazyDirichletSeriesRing(ZZ, 't')
            sage: L._laurent_poly_ring is SR                                            # needs sage.symbolic
            True
        """
        from sage.symbolic.ring import SR
        return SR

    def __init__(self, base_ring, names, sparse=True, category=None):
        r"""
        Initialize the ring.

        TESTS::

            sage: LazyDirichletSeriesRing.options.halting_precision(12)

            sage: L = LazyDirichletSeriesRing(ZZ, 't')
            sage: TestSuite(L).run()                                                    # needs sage.symbolic

            sage: L = LazyDirichletSeriesRing(QQ, 't')
            sage: TestSuite(L).run()                                                    # needs sage.symbolic

            sage: LazyDirichletSeriesRing.options._reset()  # reset the options
        """
        if base_ring.characteristic() > 0:
            raise ValueError("positive characteristic not allowed for Dirichlet series")

        self._sparse = sparse
        self._minimal_valuation = 1
        self._arity = 1
        self._internal_poly_ring = PolynomialRing(base_ring, names, sparse=sparse)

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

            sage: LazyDirichletSeriesRing(QQbar, 'z')                                   # needs sage.rings.number_field
            Lazy Dirichlet Series Ring in z over Algebraic Field
        """
        return "Lazy Dirichlet Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    @cached_method
    def one(self):
        r"""
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.one()                                                               # needs sage.symbolic
            1
            sage: ~L.one()                                                              # needs sage.symbolic
            1 + O(1/(8^z))
        """
        R = self.base_ring()
        coeff_stream = Stream_exact([R.one()], constant=R.zero(), order=1)
        return self.element_class(self, coeff_stream)

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

    def _element_constructor_(self, x=None, valuation=None, degree=None, constant=None, coefficients=None):
        r"""
        Construct a Dirichlet series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a Dirichlet series
        - ``valuation`` -- integer (optional); a lower bound for the exp of the
          valuation of the series
        - ``degree`` -- (optional) the degree when the series is ``constant``
        - ``constant`` -- (optional) the eventual constant of the series

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: R = L(3)
            sage: R                                                                     # needs sage.symbolic
            3
            sage: S = L(lambda i: i, constant=1, degree=6)
            sage: S                                                                     # needs sage.symbolic
            1 + 2/2^z + 3/3^z + 4/4^z + 5/5^z + 1/(6^z) + 1/(7^z) + 1/(8^z) + O(1/(9^z))

            sage: X = L(constant=5, degree=3)
            sage: X                                                                     # needs sage.symbolic
            5/3^z + 5/4^z + 5/5^z + O(1/(6^z))
            sage: X.valuation()                                                         # needs sage.symbolic
            log(3)
            sage: e = L(moebius)
            sage: e                                                                     # needs sage.symbolic
            1 - 1/(2^z) - 1/(3^z) - 1/(5^z) + 1/(6^z) - 1/(7^z) + O(1/(8^z))

            sage: T = L([0], constant=1)
            sage: T                                                                     # needs sage.symbolic
            1/(2^z) + 1/(3^z) + 1/(4^z) + O(1/(5^z))

            sage: U = L(constant=1)
            sage: U                                                                     # needs sage.symbolic
            1 + 1/(2^z) + 1/(3^z) + O(1/(4^z))

            sage: V = L(lambda i: i, valuation=3)
            sage: V                                                                     # needs sage.symbolic
            3/3^z + 4/4^z + 5/5^z + 6/6^z + 7/7^z + 8/8^z + 9/9^z + O(1/(10^z))

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], 4)
            sage: f                                                                     # needs sage.symbolic
            1/(4^z) + 2/5^z + 3/6^z + 4/7^z
            sage: g = L([1,3,5,7,9], 6, constant=-1)
            sage: g                                                                     # needs sage.symbolic
            1/(6^z) + 3/7^z + 5/8^z + 7/9^z + 9/10^z - 1/(11^z) - 1/(12^z)
             - 1/(13^z) + O(1/(14^z))

        TESTS::

            sage: L = LazyDirichletSeriesRing(GF(2), 'z')
            Traceback (most recent call last):
            ...
            ValueError: positive characteristic not allowed for Dirichlet series

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: D = LazyDirichletSeriesRing(QQ, 't')
            sage: d = D(L.one())                                                        # needs sage.symbolic
            sage: d                                                                     # needs sage.symbolic
            1 + 1/(2^t) + 1/(3^t) + 1/(4^t) + 1/(5^t) + 1/(6^t) + 1/(7^t) + O(1/(8^t))

            sage: R.<z> = LaurentPolynomialRing(QQ)
            sage: D = LazyDirichletSeriesRing(QQ, 't')
            sage: dd = D(coefficients=z + z^2)
            sage: dd                                                                    # needs sage.symbolic
            2 + 6/2^t + 12/3^t + 20/4^t + 30/5^t + 42/6^t + 56/7^t + O(1/(8^t))

            sage: s = D(lambda n: n)
            sage: d2 = D(s, valuation=2)                                                # needs sage.symbolic
            sage: d2                                                                    # needs sage.symbolic
            1/(2^t) + 2/3^t + 3/4^t + 4/5^t + 5/6^t + 6/7^t + 7/8^t + O(1/(9^t))

            sage: Ds = LazyDirichletSeriesRing(ZZ, 's')
            sage: m = Ds(moebius, valuation=2)
            sage: m                                                                     # needs sage.symbolic
            -1/(2^s) - 1/(3^s) - 1/(5^s) + 1/(6^s) - 1/(7^s) + O(1/(9^s))
            sage: D = LazyDirichletSeriesRing(QQ, 't')
            sage: dm = D(m)                                                             # needs sage.libs.pari
            sage: dm                                                                    # needs sage.libs.pari sage.symbolic
            -1/(2^t) - 1/(3^t) - 1/(5^t) + 1/(6^t) - 1/(7^t) + O(1/(9^t))
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
            sage: m = L.an_element()
            sage: m                                                                     # needs sage.symbolic
            1/(4^z) + 1/(5^z) + 1/(6^z) + O(1/(7^z))
        """
        c = self.base_ring().an_element()
        return self.element_class(self, Stream_exact([], constant=c, order=4))

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: l = L.some_elements()
            sage: l                                                                     # needs sage.symbolic
            [0, 1,
             1/(4^z) + 1/(5^z) + 1/(6^z) + O(1/(7^z)),
             1/(2^z) - 1/(3^z) + 2/4^z - 2/5^z + 3/6^z - 3/7^z + 4/8^z - 4/9^z,
             1/(2^z) - 1/(3^z) + 2/4^z - 2/5^z + 3/6^z - 3/7^z + 4/8^z - 4/9^z + 1/(10^z) + 1/(11^z) + 1/(12^z) + O(1/(13^z)),
             1 + 4/2^z + 9/3^z + 16/4^z + 25/5^z + 36/6^z + 49/7^z + O(1/(8^z))]

            sage: L = LazyDirichletSeriesRing(QQ, 'z')
            sage: l = L.some_elements()
            sage: l                                                                     # needs sage.symbolic
            [0, 1,
             1/2/4^z + 1/2/5^z + 1/2/6^z + O(1/(7^z)),
             1/2 - 1/2/2^z + 2/3^z - 2/4^z + 1/(6^z) - 1/(7^z) + 42/8^z + 2/3/9^z,
             1/2 - 1/2/2^z + 2/3^z - 2/4^z + 1/(6^z) - 1/(7^z) + 42/8^z + 2/3/9^z + 1/2/10^z + 1/2/11^z + 1/2/12^z + O(1/(13^z)),
             1 + 4/2^z + 9/3^z + 16/4^z + 25/5^z + 36/6^z + 49/7^z + O(1/(8^z))]
        """
        R = self.base_ring()
        some_numbers = [c for c, _ in zip(R.some_elements(), range(9))]
        elts = [self.zero(), self.one(), self.an_element(),
                self(some_numbers),
                self(some_numbers, constant=R.an_element()),
                self(lambda n: n**2)]
        return elts

    def _monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: m = L._monomial(5, 3); m                                              # needs sage.symbolic
            5/3^z
        """
        try:
            L = self._laurent_poly_ring
            return L(c) * L(n) ** -L(self.variable_name())
        except (ValueError, TypeError):
            return '({})/{}^{}'.format(self.base_ring()(c), n, self.variable_name())


def _skip_leading_zeros(iterator):
    """
    Return an iterator which discards all leading zeros.

    EXAMPLES::

        sage: from sage.rings.lazy_series_ring import _skip_leading_zeros
        sage: it = map(lambda x: 0 if x < 10 else x, NN)
        sage: [x for x, _ in zip(_skip_leading_zeros(it), range(10))]
        [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

        sage: # needs sage.rings.finite_rings
        sage: it = map(GF(3), NN)
        sage: [x for x, _ in zip(it, range(10))]
        [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
        sage: it = map(GF(3), NN)
        sage: [x for x, _ in zip(_skip_leading_zeros(it), range(10))]
        [1, 2, 0, 1, 2, 0, 1, 2, 0, 1]
    """
    while True:
        c = next(iterator)
        if c:
            yield c
            break
    yield from iterator
