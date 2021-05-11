r"""
Lazy Laurent Series Rings

The ring of lazy Laurent series over a ring has usual arithmetic operations,
but it is actually not a ring in the usual sense since every
arithmetic operation gives a new series.

EXAMPLES:

The definition of Laurent series rings is not initially imported into the
global namespace. You need to import it explicitly to use it::

    sage: L.<z> = LazyLaurentSeriesRing(QQ)
    sage: L.category()
    Category of magmas and additive magmas
    sage: 1/(1 - z)
    1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
    sage: 1/(1 - z) == 1/(1 - z)
    True

Lazy Laurent series ring over a finite field::

    sage: L.<z> = LazyLaurentSeriesRing(GF(3)); L
    Lazy Laurent Series Ring in z over Finite Field of size 3
    sage: e = 1/(1 + z)
    sage: e.coefficient(100)
    1
    sage: e.coefficient(100).parent()
    Finite Field of size 3

Generating functions of integer sequences are Laurent series over the integer
ring::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ); L
    Lazy Laurent Series Ring in z over Integer Ring
    sage: 1/(1 - 2*z)^3
    1 + 6*z + 24*z^2 + 80*z^3 + 240*z^4 + 672*z^5 + 1792*z^6 + ...

Power series can be defined recursively::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: L.series(lambda s,n: (1 + z*s^2)[n], valuation=0)
    1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...

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
import random

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets

from sage.misc.cachefunc import cached_method

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_generic

from .lazy_laurent_series import LazyLaurentSeries
from .lazy_laurent_series_operator import (
    LazyLaurentSeriesOperator_gen,
    LazyLaurentSeriesOperator_constant,
    LazyLaurentSeriesOperator_list,
    LazyLaurentSeriesOperator_polynomial
)


class LazyLaurentSeriesRing(UniqueRepresentation, Parent):
    """
    Lazy Laurent series ring.

    INPUT:

    - ``base_ring`` -- base ring of this Laurent series ring

    - ``names`` -- name of the generator of this Laurent series ring

    EXAMPLES::

        sage: LazyLaurentSeriesRing(ZZ, 't')
        Lazy Laurent Series Ring in t over Integer Ring
    """
    Element = LazyLaurentSeries

    def __init__(self, base_ring, names, category=None):
        """
        Initialize.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip='_test_elements')
        """
        Parent.__init__(self, base=base_ring, names=names,
                        category=MagmasAndAdditiveMagmas().or_subcategory(category))

    def _repr_(self):
        """
        String representation of this Laurent series ring.

        EXAMPLES::

            sage: LazyLaurentSeriesRing(GF(2), 'z')
            Lazy Laurent Series Ring in z over Finite Field of size 2
        """
        return "Lazy Laurent Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    @cached_method
    def gen(self, n=0):
        """
        Return the generator of this Laurent series ring.

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

        op = LazyLaurentSeriesOperator_gen(self)
        c = (self.base_ring().zero(), 2)
        return self.element_class(self, coefficient=op, valuation=1, constant=c)

    def ngens(self):
        """
        Return the number of generators of this Laurent series ring.

        This is always 1.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.ngens()
            1
        """
        return 1

    def gens(self):
        """
        Return the tuple of the generator.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: 1/(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        return (self.gen(),)

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

        if isinstance(S, (PolynomialRing_general, LaurentPolynomialRing_generic)) and S.ngens() == 1:
            def make_series_from(poly):
                op = LazyLaurentSeriesOperator_polynomial(self, poly)
                a = poly.valuation()
                c = (self.base_ring().zero(), poly.degree() + 1)
                return self.element_class(self, coefficient=op, valuation=a, constant=c)

            return SetMorphism(Hom(S, self, Sets()), make_series_from)

        return False

    def _element_constructor_(self, x):
        """
        Construct a Laurent series from ``x``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1
        """
        R = self.base_ring()

        op = LazyLaurentSeriesOperator_constant(self, R(x))

        return self.element_class(self, coefficient=op, valuation=0, constant=(R.zero(), 1))

    def _an_element_(self):
        """
        Return a Laurent series in this ring.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.an_element()  # random
            z^-10 + z^-9 + z^-8 + z^-7 + z^-6 + z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + ...
        """
        N = 10

        e = self.base_ring().an_element()

        def r(s, i):
            return self.base_ring().an_element()

        n = random.randint(-N,N)
        m = random.randint(0,N)

        return self.element_class(self, coefficient=r, valuation=n, constant=(e,n+m))

    def one(self):
        """
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.one()
            1
        """
        return self._element_constructor_(1)

    def zero(self):
        """
        Return the zero series.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self._element_constructor_(0)

    def series(self, coefficient, valuation, constant=None):
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

        Let ``valuation`` be `n`. All coefficients of index below `n` are zero.
        If ``constant`` is ``None``, then the ``coefficient`` function is responsible to
        compute the values of all coefficients of index `\ge n`. If
        ``constant`` is a pair `(c,m)`, then the ``coefficient`` function is responsible
        to compute the values of all coefficients of index `\ge n` and `< m`
        and all the coefficients of index `\ge m` is the constant `c`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.series(lambda s, i: i, 5, (1,10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + ...

            sage: def g(s, i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return s.coefficient(i - 1) + i
            sage: e = L.series(g, -5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + ...
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + ...
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
        simply omitted if it is zero.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5)
            sage: f
            z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2
            sage: g = L.series([1,3,5,7,9], 5, -1)
            sage: g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + ...
        """
        if isinstance(coefficient, (tuple, list)):
            if isinstance(constant, tuple):
                constant = constant[0]
            if constant is None:
                constant = self.base_ring().zero()
            elif constant not in self.base_ring():
                raise ValueError("constant is not an element of the base ring")
            constant = (constant, valuation + len(coefficient))
            coefficient = LazyLaurentSeriesOperator_list(self, coefficient, valuation)
        elif constant is not None:
            try:
                c,m = constant
            except TypeError:
                raise TypeError('not a tuple')

            if valuation > m and c: # weird case
                raise ValueError('inappropriate valuation')

            constant = (self.base_ring()(c), m)

        return self.element_class(self, coefficient=coefficient, valuation=valuation, constant=constant)
