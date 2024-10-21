# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.modules
r"""
Ring of Laurent Polynomials (base class)

If `R` is a commutative ring, then the ring of Laurent polynomials in `n`
variables over `R` is `R[x_1^{\pm 1}, x_2^{\pm 1}, \ldots, x_n^{\pm 1}]`.

AUTHORS:

- David Roe (2008-2-23): created
- David Loeffler (2009-07-10): cleaned up docstrings
"""
# ****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>,
#                          William Stein <wstein@gmail.com>,
#                          Mike Hansen <mhansen@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ring import CommutativeRing
from sage.structure.parent import Parent
from sage.combinat.integer_vector import IntegerVectors

class LaurentPolynomialRing_generic(CommutativeRing, Parent):
    """
    Laurent polynomial ring (base class).

    EXAMPLES:

    This base class inherits from :class:`~sage.rings.ring.CommutativeRing`.
    Since :issue:`11900`, it is also initialised as such::

        sage: R.<x1,x2> = LaurentPolynomialRing(QQ)
        sage: R.category()
        Join of Category of unique factorization domains
            and Category of algebras with basis
                over (number fields and quotient fields and metric spaces)
            and Category of commutative algebras
                over (number fields and quotient fields and metric spaces)
            and Category of infinite sets
        sage: TestSuite(R).run()
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ, 2, 'x')
            sage: R == loads(dumps(R))
            True
        """
        self._n = R.ngens()
        self._R = R
        names = R.variable_names()
        self._one_element = self.element_class(self, R.one())
        CommutativeRing.__init__(self, R.base_ring(), names=names,
                                 category=R.category())
        ernames = []
        for n in names:
            ernames.append(n)
            ernames.append(n + "inv")
        ER = PolynomialRing(R.base_ring(), ernames)
        self._extended_ring = ER
        self._extended_ring_ideal = ER.ideal([ER.gen(2*i)*ER.gen(2*i+1)-1 for i in range(self._n)])

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').ngens()
            2
            sage: LaurentPolynomialRing(QQ, 1, 'x').ngens()
            1
        """
        return self._n

    def gen(self, i=0):
        r"""
        Return the `i`-th generator of ``self``.  If `i` is not specified, then
        the first generator will be returned.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').gen()
            x0
            sage: LaurentPolynomialRing(QQ, 2, 'x').gen(0)
            x0
            sage: LaurentPolynomialRing(QQ, 2, 'x').gen(1)
            x1

        TESTS::

            sage: LaurentPolynomialRing(QQ, 2, 'x').gen(3)
            Traceback (most recent call last):
            ...
            ValueError: generator not defined
        """
        if i < 0 or i >= self._n:
            raise ValueError("generator not defined")
        try:
            return self.__generators[i]
        except AttributeError:
            self.__generators = tuple(self(x) for x in self._R.gens())
            return self.__generators[i]

    def variable_names_recursive(self, depth=infinity):
        r"""
        Return the list of variable names of this ring and its base rings,
        as if it were a single multi-variate Laurent polynomial.

        INPUT:

        - ``depth`` -- integer or :mod:`Infinity <sage.rings.infinity>`

        OUTPUT: a tuple of strings

        EXAMPLES::

            sage: T = LaurentPolynomialRing(QQ, 'x')
            sage: S = LaurentPolynomialRing(T, 'y')
            sage: R = LaurentPolynomialRing(S, 'z')
            sage: R.variable_names_recursive()
            ('x', 'y', 'z')
            sage: R.variable_names_recursive(2)
            ('y', 'z')
        """
        if depth <= 0:
            return ()
        if depth == 1:
            return self.variable_names()
        my_vars = self.variable_names()
        try:
            return self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
        except AttributeError:
            return my_vars

    def is_integral_domain(self, proof=True):
        """
        Return ``True`` if ``self`` is an integral domain.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').is_integral_domain()
            True

        The following used to fail; see :issue:`7530`::

            sage: L = LaurentPolynomialRing(ZZ, 'X')
            sage: L['Y']
            Univariate Polynomial Ring in Y over
             Univariate Laurent Polynomial Ring in X over Integer Ring
        """
        return self.base_ring().is_integral_domain(proof)

    def is_noetherian(self):
        """
        Return ``True`` if ``self`` is Noetherian.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').is_noetherian()
            True
        """
        return self.base_ring().is_noetherian()

    def construction(self):
        """
        Return the construction of ``self``.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x,y').construction()
            (LaurentPolynomialFunctor,
             Univariate Laurent Polynomial Ring in x over Rational Field)
        """
        from sage.categories.pushout import LaurentPolynomialFunctor
        from .laurent_polynomial_ring import LaurentPolynomialRing

        vars = self.variable_names()
        if len(vars) == 1:
            return LaurentPolynomialFunctor(vars[0], False), self.base_ring()
        else:
            return LaurentPolynomialFunctor(vars[-1], True), LaurentPolynomialRing(self.base_ring(), vars[:-1])

    def completion(self, p=None, prec=20, extras=None):
        r"""
        Return the completion of ``self``.

        Currently only implemented for the ring of formal Laurent series.
        The ``prec`` variable controls the precision used in the
        Laurent series ring. If ``prec`` is `\infty`, then this
        returns a :class:`LazyLaurentSeriesRing`.

        EXAMPLES::

            sage: P.<x> = LaurentPolynomialRing(QQ); P
            Univariate Laurent Polynomial Ring in x over Rational Field
            sage: PP = P.completion(x); PP
            Laurent Series Ring in x over Rational Field
            sage: f = 1 - 1/x
            sage: PP(f)
            -x^-1 + 1
            sage: g = 1 / PP(f); g
            -x - x^2 - x^3 - x^4 - x^5 - x^6 - x^7 - x^8 - x^9 - x^10 - x^11
             - x^12 - x^13 - x^14 - x^15 - x^16 - x^17 - x^18 - x^19 - x^20 + O(x^21)
            sage: 1 / g
            -x^-1 + 1 + O(x^19)

            sage: # needs sage.combinat
            sage: PP = P.completion(x, prec=oo); PP
            Lazy Laurent Series Ring in x over Rational Field
            sage: g = 1 / PP(f); g
            -x - x^2 - x^3 + O(x^4)
            sage: 1 / g == f
            True

        TESTS:

        Check that the precision is taken into account (:issue:`24431`)::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: L.completion('x', 100).default_prec()
            100
            sage: L.completion('x', 20).default_prec()
            20
        """
        if p is None or str(p) == self._names[0] and self._n == 1:
            if prec == float('inf'):
                from sage.rings.lazy_series_ring import LazyLaurentSeriesRing
                sparse = self.polynomial_ring().is_sparse()
                return LazyLaurentSeriesRing(self.base_ring(), names=(self._names[0],), sparse=sparse)
            from sage.rings.laurent_series_ring import LaurentSeriesRing
            R = self.polynomial_ring().completion(self._names[0], prec)
            return LaurentSeriesRing(R)

        raise TypeError("cannot complete %s with respect to %s" % (self, p))

    def remove_var(self, var):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,'x,y,z')
            sage: R.remove_var('x')
            Multivariate Laurent Polynomial Ring in y, z over Rational Field
            sage: R.remove_var('x').remove_var('y')
            Univariate Laurent Polynomial Ring in z over Rational Field
        """
        from .laurent_polynomial_ring import LaurentPolynomialRing

        vars = list(self.variable_names())
        vars.remove(str(var))
        return LaurentPolynomialRing(self.base_ring(), vars)

    def _coerce_map_from_(self, R):
        """
        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: L.coerce_map_from(QQ)
            Generic morphism:
              From: Rational Field
              To:   Multivariate Laurent Polynomial Ring in x, y over Rational Field

        Let us check that coercion between Laurent Polynomials over
        different base rings works (:issue:`15345`)::

            sage: R = LaurentPolynomialRing(ZZ, 'x')
            sage: T = LaurentPolynomialRing(QQ, 'x')
            sage: R.gen() + 3*T.gen()
            4*x
        """
        if R is self._R:
            return self._generic_coerce_map(R)
        f = self._coerce_map_via([self._R], R)
        if f is not None:
            return f
        if (isinstance(R, LaurentPolynomialRing_generic)
            and self._R.has_coerce_map_from(R._R)):
            return self._generic_coerce_map(R)

    def __eq__(self, right):
        """
        Check whether ``self`` is equal to ``right``.

        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,'x,y,z')
            sage: P = LaurentPolynomialRing(ZZ,'x,y,z')
            sage: Q = LaurentPolynomialRing(QQ,'x,y')

            sage: R == R
            True
            sage: R == Q
            False
            sage: Q == P
            False
            sage: P == R
            False
        """
        if type(self) is not type(right):
            return False
        return self._R == right._R

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,'x,y,z')
            sage: P = LaurentPolynomialRing(ZZ,'x,y,z')
            sage: Q = LaurentPolynomialRing(QQ,'x,y')

            sage: R != R
            False
            sage: R != Q
            True
            sage: Q != P
            True
            sage: P != R
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: h1 = hash(LaurentPolynomialRing(ZZ,'x,y,z'))
            sage: h2 = hash(LaurentPolynomialRing(ZZ,'x,y,z'))
            sage: h3 = hash(LaurentPolynomialRing(QQ,'x,y,z'))
            sage: h4 = hash(LaurentPolynomialRing(ZZ,'x,y'))
            sage: h1 == h2 and h1 != h3 and h1 != h4
            True
        """
        return hash(self._R) ^ 12059065606945654693

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(LaurentPolynomialRing(QQ, 2, 'x'))
            \Bold{Q}[x_{0}^{\pm 1}, x_{1}^{\pm 1}]
        """
        from sage.misc.latex import latex

        vars = ', '.join(a + r'^{\pm 1}' for a in self.latex_variable_names())
        return "%s[%s]" % (latex(self.base_ring()), vars)

    def _ideal_class_(self, n=0):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x')._ideal_class_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # One may eventually want ideal classes in these guys.
        raise NotImplementedError

    def ideal(self, *args, **kwds):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').ideal([1])
            Ideal (1) of Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field

        TESTS:

        check that :issue:`26421` is fixed::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: P.<x> = PolynomialRing(R)
            sage: p = x-t
            sage: p.content_ideal()    # indirect doctest
            Ideal (-t, 1) of Univariate Laurent Polynomial Ring in t over Integer Ring
        """
        from sage.rings.polynomial.laurent_polynomial_ideal import LaurentPolynomialIdeal
        return LaurentPolynomialIdeal(self, *args, **kwds)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: T.<t> = ZZ[]
            sage: K.<i> = NumberField(t^2 + 1)
            sage: L.<x,y> = LaurentPolynomialRing(K)
            sage: L._is_valid_homomorphism_(K, (K(1/2), K(3/2)))
            True
            sage: Q5 = Qp(5); i5 = Q5(-1).sqrt()                                                    # needs sage.rings.padics
            sage: L._is_valid_homomorphism_(Q5, (Q5(1/2), Q5(3/2)))  # no coercion                  # needs sage.rings.padics
            False
            sage: L._is_valid_homomorphism_(Q5, (Q5(1/2), Q5(3/2)), base_map=K.hom([i5]))           # needs sage.rings.padics
            True
        """
        if base_map is None and not codomain.has_coerce_map_from(self.base_ring()):
            # we need that elements of the base ring
            # canonically coerce into codomain.
            return False
        for a in im_gens:
            # in addition, the image of each generator must be invertible.
            if not a.is_unit():
                return False
        return True

    def term_order(self):
        """
        Return the term order of ``self``.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').term_order()
            Degree reverse lexicographic term order
        """
        return self._R.term_order()

    def is_finite(self):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').is_finite()
            False
        """
        return False

    def is_field(self, proof=True):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').is_field()
            False
        """
        return False

    def polynomial_ring(self):
        """
        Return the polynomial ring associated with ``self``.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').polynomial_ring()
            Multivariate Polynomial Ring in x0, x1 over Rational Field
            sage: LaurentPolynomialRing(QQ, 1, 'x').polynomial_ring()
            Multivariate Polynomial Ring in x over Rational Field
        """
        return self._R

    def characteristic(self):
        """
        Return the characteristic of the base ring.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').characteristic()
            0
            sage: LaurentPolynomialRing(GF(3), 2, 'x').characteristic()
            3
        """
        return self.base_ring().characteristic()

    def krull_dimension(self):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').krull_dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def random_element(self, min_valuation=-2, max_degree=2, *args, **kwds):
        """
        Return a random polynomial with degree at most ``max_degree`` and
        lowest valuation at least ``min_valuation``.

        Uses the random sampling from the base polynomial ring then divides out
        by a monomial to ensure correct ``max_degree`` and ``min_valuation``.

        INPUT:

        - ``min_valuation`` -- integer (default: `-2`); the
          minimal allowed valuation of the polynomial

        - ``max_degree`` -- integer (default: `2`); the
          maximal allowed degree of the polynomial

        - ``*args``, ``**kwds`` -- passed to the random element generator of the
          base polynomial ring and base ring itself

        EXAMPLES::

            sage: L.<x> = LaurentPolynomialRing(QQ)
            sage: f = L.random_element()
            sage: f.degree() <= 2
            True
            sage: f.valuation() >= -2
            True
            sage: f.parent() is L
            True

        ::

            sage: L = LaurentPolynomialRing(ZZ, 2, 'x')
            sage: f = L.random_element(10, 20)
            sage: f.degree() <= 20
            True
            sage: f.valuation() >= 10
            True
            sage: f.parent() is L
            True

        ::

            sage: L = LaurentPolynomialRing(GF(13), 3, 'x')
            sage: f = L.random_element(-10, -1)
            sage: f.degree() <= -1
            True
            sage: f.valuation() >= -10
            True
            sage: f.parent() is L
            True

        ::

            sage: L.<x, y> = LaurentPolynomialRing(RR)
            sage: f = L.random_element()
            sage: f.degree() <= 2
            True
            sage: f.valuation() >= -2
            True
            sage: f.parent() is L
            True

        ::

            sage: L = LaurentPolynomialRing(QQbar, 5, 'x')
            sage: f = L.random_element(-1, 1)
            sage: f = L.random_element(-1, 1)
            sage: f.degree() <= 1
            True
            sage: f.valuation() >= -1
            True
            sage: f.parent() is L
            True

        TESTS:

        Ensure everything works for the multivariate case with only
        one generator::

            sage: L = LaurentPolynomialRing(ZZ, 1, 'x')
            sage: f = L.random_element(10, 20)
            sage: f.degree() <= 20
            True
            sage: f.valuation() >= 10
            True
            sage: f.parent() is L
            True

        Test for constructions which use multivariate polynomial rings::

            sage: rings = [RR, QQ, ZZ, GF(13), GF(7^3)]
            sage: for ring in rings:
            ....:     d = randint(1, 6)
            ....:     t = randint(5, 20)
            ....:     L = LaurentPolynomialRing(ring, d, 'x')
            ....:     for _ in range(100):
            ....:         n, m = randint(-10, 10), randint(-10, 10)
            ....:         if n > m:
            ....:             n, m = m, n
            ....:         f = L.random_element(n, m, terms=t)
            ....:         if f.is_zero(): continue # the zero polynomial is defined to have degree -1
            ....:         assert len(list(f)) <= t
            ....:         assert f.degree() <= m
            ....:         assert f.valuation() >= n

        Test for constructions which use univariate polynomial rings::

            sage: rings = [RR, QQ, ZZ, GF(13), GF(7^3)]
            sage: for ring in rings:
            ....:     L.<x> = LaurentPolynomialRing(ring)
            ....:     for _ in range(100):
            ....:         n, m = randint(-10, 10), randint(-10, 10)
            ....:         if n > m:
            ....:             n, m = m, n
            ....:         f = L.random_element(n, m)
            ....:         if f.is_zero(): continue # the zero polynomial is defined to have degree -1
            ....:         for x in L.gens():
            ....:             assert f.degree() <= m
            ....:             assert f.valuation() >= n

        The ``max_degree`` must be greater than or equal to ``min_valuation``::

            sage: L.<x> = LaurentPolynomialRing(QQ)
            sage: f = L.random_element(1, -1)
            Traceback (most recent call last):
            ...
            ValueError: `max_degree` must be greater than or equal to `min_valuation`
        """
        # Ensure the degree parameters are sensible
        if max_degree < min_valuation:
            raise ValueError("`max_degree` must be greater than or equal to `min_valuation`")

        # Sample a polynomial in the base ring of degree `max_degree - min_valuation`
        abs_deg = (max_degree - min_valuation)
        f_rand = self._R.random_element(degree=abs_deg, *args, **kwds)

        # Cast this polynomial back the ``self``
        f = self(f_rand)

        # For the univariate case we simply shift by x**min_valuation
        if self._n == 1:
            # When there is one generator we either have a univariate class
            # and can shift by min_valuation
            try:
                return f.shift(min_valuation)
            # Or we have a multivariate class with one variable, which does not
            # have the shift method
            except AttributeError:
                return f * self.gen() ** min_valuation

        # For the multivariate case, we sample a single monomial of degree
        # exactly min_valuation and then use this to shift the polynomial
        # into the correct bounds
        s = self.monomial(*IntegerVectors(abs(min_valuation), self._n).random_element())
        if min_valuation < 0:
            s = ~s
        return f * s

    def is_exact(self):
        """
        Return ``True`` if the base ring is exact.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ, 2, 'x').is_exact()
            True
            sage: LaurentPolynomialRing(RDF, 2, 'x').is_exact()
            False
        """
        return self.base_ring().is_exact()

    def change_ring(self, base_ring=None, names=None, sparse=False, order=None):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ, 2, 'x')
            sage: R.change_ring(ZZ)
            Multivariate Laurent Polynomial Ring in x0, x1 over Integer Ring

        Check that the distinction between a univariate ring and a multivariate ring with one
        generator is preserved::

            sage: P.<x> = LaurentPolynomialRing(QQ, 1)
            sage: P
            Multivariate Laurent Polynomial Ring in x over Rational Field
            sage: K.<i> = CyclotomicField(4)                                                        # needs sage.rings.number_field
            sage: P.change_ring(K)                                                                  # needs sage.rings.number_field
            Multivariate Laurent Polynomial Ring in x over
             Cyclotomic Field of order 4 and degree 2
        """
        from .laurent_polynomial_ring import LaurentPolynomialRing, LaurentPolynomialRing_univariate

        if base_ring is None:
            base_ring = self.base_ring()
        if names is None:
            names = self.variable_names()
        if isinstance(self, LaurentPolynomialRing_univariate):
            return LaurentPolynomialRing(base_ring, names[0], sparse=sparse)

        if order is None:
            order = self.polynomial_ring().term_order()
        return LaurentPolynomialRing(base_ring, self._n, names, order=order)

    def fraction_field(self):
        """
        The fraction field is the same as the fraction field of the
        polynomial ring.

        EXAMPLES::

            sage: L.<x> = LaurentPolynomialRing(QQ)
            sage: L.fraction_field()
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
            sage: (x^-1 + 2) / (x - 1)
            (2*x + 1)/(x^2 - x)
        """
        return self.polynomial_ring().fraction_field()
