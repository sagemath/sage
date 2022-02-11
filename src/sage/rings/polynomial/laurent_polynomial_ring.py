r"""
Ring of Laurent Polynomials

If `R` is a commutative ring, then the ring of Laurent polynomials in `n`
variables over `R` is `R[x_1^{\pm 1}, x_2^{\pm 1}, \ldots, x_n^{\pm 1}]`.
We implement it as a quotient ring

.. MATH::

    R[x_1, y_1, x_2, y_2, \ldots, x_n, y_n] / (x_1 y_1 - 1, x_2 y_2 - 1, \ldots, x_n y_n - 1).

TESTS::

    sage: P.<q> = LaurentPolynomialRing(QQ)
    sage: qi = q^(-1)
    sage: qi in P
    True
    sage: P(qi)
    q^-1

    sage: A.<Y> = QQ[]
    sage: R.<X> = LaurentPolynomialRing(A)
    sage: matrix(R,2,2,[X,0,0,1])
    [X 0]
    [0 1]

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

from sage.structure.element import parent
from sage.structure.parent import Parent
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.latex import latex
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_mpair, LaurentPolynomial_univariate
from sage.rings.ring import CommutativeRing

import sage.rings.polynomial.laurent_polynomial_ideal as lp_ideal

def is_LaurentPolynomialRing(R):
    """
    Returns True if and only if R is a Laurent polynomial ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import is_LaurentPolynomialRing
        sage: P = PolynomialRing(QQ,2,'x')
        sage: is_LaurentPolynomialRing(P)
        False

        sage: R = LaurentPolynomialRing(QQ,3,'x')
        sage: is_LaurentPolynomialRing(R)
        True
    """
    return isinstance(R, LaurentPolynomialRing_generic)

_cache = {}
def LaurentPolynomialRing(base_ring, *args, **kwds):
    r"""
    Return the globally unique univariate or multivariate Laurent polynomial
    ring with given properties and variable name or names.

    There are four ways to call the Laurent polynomial ring constructor:

    1. ``LaurentPolynomialRing(base_ring, name,    sparse=False)``
    2. ``LaurentPolynomialRing(base_ring, names,   order='degrevlex')``
    3. ``LaurentPolynomialRing(base_ring, name, n, order='degrevlex')``
    4. ``LaurentPolynomialRing(base_ring, n, name, order='degrevlex')``

    The optional arguments sparse and order *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:

    - ``base_ring`` -- a commutative ring
    - ``name`` -- a string
    - ``names`` -- a list or tuple of names, or a comma separated string
    - ``n`` -- a positive integer
    - ``sparse`` -- bool (default: False), whether or not elements are sparse
    - ``order`` -- string or
      :class:`~sage.rings.polynomial.term_order.TermOrder`, e.g.,

        - ``'degrevlex'`` (default) -- degree reverse lexicographic
        - ``'lex'`` -- lexicographic
        - ``'deglex'`` -- degree lexicographic
        - ``TermOrder('deglex',3) + TermOrder('deglex',3)`` -- block ordering

    OUTPUT:

    ``LaurentPolynomialRing(base_ring, name, sparse=False)`` returns a
    univariate Laurent polynomial ring; all other input formats return a
    multivariate Laurent polynomial ring.

    UNIQUENESS and IMMUTABILITY: In Sage there is exactly one
    single-variate Laurent polynomial ring over each base ring in each choice
    of variable and sparseness.  There is also exactly one multivariate
    Laurent polynomial ring over each base ring for each choice of names of
    variables and term order.

    ::

        sage: R.<x,y> = LaurentPolynomialRing(QQ,2); R
        Multivariate Laurent Polynomial Ring in x, y over Rational Field
        sage: f = x^2 - 2*y^-2

    You can't just globally change the names of those variables.
    This is because objects all over Sage could have pointers to
    that polynomial ring.

    ::

        sage: R._assign_names(['z','w'])
        Traceback (most recent call last):
        ...
        ValueError: variable names cannot be changed after object creation.


    EXAMPLES:

    1. ``LaurentPolynomialRing(base_ring, name, sparse=False)``

       ::

           sage: LaurentPolynomialRing(QQ, 'w')
           Univariate Laurent Polynomial Ring in w over Rational Field

       Use the diamond brackets notation to make the variable
       ready for use after you define the ring::

           sage: R.<w> = LaurentPolynomialRing(QQ)
           sage: (1 + w)^3
           1 + 3*w + 3*w^2 + w^3

       You must specify a name::

           sage: LaurentPolynomialRing(QQ)
           Traceback (most recent call last):
           ...
           TypeError: you must specify the names of the variables

           sage: R.<abc> = LaurentPolynomialRing(QQ, sparse=True); R
           Univariate Laurent Polynomial Ring in abc over Rational Field

           sage: R.<w> = LaurentPolynomialRing(PolynomialRing(GF(7),'k')); R
           Univariate Laurent Polynomial Ring in w over Univariate Polynomial Ring in k over Finite Field of size 7

       Rings with different variables are different::

           sage: LaurentPolynomialRing(QQ, 'x') == LaurentPolynomialRing(QQ, 'y')
           False

    2. ``LaurentPolynomialRing(base_ring, names,   order='degrevlex')``

       ::

           sage: R = LaurentPolynomialRing(QQ, 'a,b,c'); R
           Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

           sage: S = LaurentPolynomialRing(QQ, ['a','b','c']); S
           Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

           sage: T = LaurentPolynomialRing(QQ, ('a','b','c')); T
           Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

       All three rings are identical.

       ::

           sage: (R is S) and  (S is T)
           True

       There is a unique Laurent polynomial ring with each term order::

           sage: R = LaurentPolynomialRing(QQ, 'x,y,z', order='degrevlex'); R
           Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
           sage: S = LaurentPolynomialRing(QQ, 'x,y,z', order='invlex'); S
           Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
           sage: S is LaurentPolynomialRing(QQ, 'x,y,z', order='invlex')
           True
           sage: R == S
           False


    3. ``LaurentPolynomialRing(base_ring, name, n, order='degrevlex')``

       If you specify a single name as a string and a number of
       variables, then variables labeled with numbers are created.

       ::

           sage: LaurentPolynomialRing(QQ, 'x', 10)
           Multivariate Laurent Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 over Rational Field

           sage: LaurentPolynomialRing(GF(7), 'y', 5)
           Multivariate Laurent Polynomial Ring in y0, y1, y2, y3, y4 over Finite Field of size 7

           sage: LaurentPolynomialRing(QQ, 'y', 3, sparse=True)
           Multivariate Laurent Polynomial Ring in y0, y1, y2 over Rational Field

       By calling the
       :meth:`~sage.structure.category_object.CategoryObject.inject_variables`
       method, all those variable names are available for interactive use::

           sage: R = LaurentPolynomialRing(GF(7),15,'w'); R
           Multivariate Laurent Polynomial Ring in w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14 over Finite Field of size 7
           sage: R.inject_variables()
           Defining w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14
           sage: (w0 + 2*w8 + w13)^2
           w0^2 + 4*w0*w8 + 4*w8^2 + 2*w0*w13 + 4*w8*w13 + w13^2
    """
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing

    R = PolynomialRing(base_ring, *args, **kwds)
    if R in _cache:
        return _cache[R]   # put () here to re-enable weakrefs

    if is_PolynomialRing(R):
        # univariate case
        P = LaurentPolynomialRing_univariate(R)
    else:
        assert is_MPolynomialRing(R)
        P = LaurentPolynomialRing_mpair(R)

    _cache[R] = P
    return P

def _split_dict_(D, indices, group_by=None):
    r"""
    Split the dictionary ``D`` by ``indices`` and ``group_by``.

    INPUT:

    - ``D`` -- a dictionary.

    - ``indices`` -- a tuple or list of nonnegative integers.

    - ``group_by`` -- a tuple or list of nonnegative integers.
      If this is ``None`` (default), then no grouping is done.

    OUTPUT:

    A dictionary.

    TESTS::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import _split_dict_
        sage: D = {(0,0,0,0): 'a', (1,0,0,0): 'b',
        ....:      (1,0,0,2): 'c', (1,2,0,3): 'd'}
        sage: _split_dict_(D, [1,0,3])
        {(0, 0, 0): 'a', (0, 1, 0): 'b', (0, 1, 2): 'c', (2, 1, 3): 'd'}
        sage: _split_dict_(D, [2,3], [0,1])
        {(0, 0): {(0, 0): 'a'},
         (1, 0): {(0, 0): 'b', (0, 2): 'c'},
         (1, 2): {(0, 3): 'd'}}
        sage: _split_dict_(D, [3,1], [0])
        {(0,): {(0, 0): 'a'}, (1,): {(0, 0): 'b', (2, 0): 'c', (3, 2): 'd'}}

        sage: _split_dict_(D, [0,None,1,3])
        {(0, 0, 0, 0): 'a', (1, 0, 0, 0): 'b',
         (1, 0, 0, 2): 'c', (1, 0, 2, 3): 'd'}
        sage: _split_dict_(D, [0,1], [None,3,None])
        {(0, 0, 0): {(0, 0): 'a', (1, 0): 'b'},
         (0, 2, 0): {(1, 0): 'c'},
         (0, 3, 0): {(1, 2): 'd'}}
        sage: _split_dict_(D, [None,3,1], [0,None])
        {(0, 0): {(0, 0, 0): 'a'},
         (1, 0): {(0, 0, 0): 'b', (0, 2, 0): 'c',
                     (0, 3, 2): 'd'}}

        sage: _split_dict_(D, [0,1])
        Traceback (most recent call last):
        ...
        SplitDictError: split not possible
        sage: _split_dict_(D, [0], [1])
        Traceback (most recent call last):
        ...
        SplitDictError: split not possible
        sage: _split_dict_({}, [])
        {}
    """
    if not D:
        return {}
    if group_by is None:
        group_by = tuple()

    class SplitDictError(ValueError):
        pass
    def get(T, i):
        return T[i] if i is not None else 0
    def extract(T, indices):
        return tuple(get(T, i) for i in indices)

    remaining = sorted(set(range(len(next(iter(D)))))
                       - set(indices) - set(group_by))
    result = {}
    for K, V in D.items():
        if not all(r == 0 for r in extract(K, remaining)):
            raise SplitDictError('split not possible')
        G = extract(K, group_by)
        I = extract(K, indices)
        result.setdefault(G, dict()).update({I: V})
    if not group_by:
        return result.popitem()[1]
    else:
        return result

def _split_laurent_polynomial_dict_(P, M, d):
    r"""
    Helper function for splitting a multivariate Laurent polynomial
    during conversion.

    INPUT:

    - ``P`` -- the parent to which we want to convert.

    - ``M`` -- the parent from which we want to convert.

    - ``d`` -- a dictionary mapping tuples (representing the exponents)
      to their coefficients. This is the dictionary corresponding to
      an element of ``M``.

    OUTPUT:

    A dictionary corresponding to an element of ``P``.

    TESTS::

        sage: L.<a, b, c, d> = LaurentPolynomialRing(ZZ)
        sage: M = LaurentPolynomialRing(ZZ, 'c, d')
        sage: N = LaurentPolynomialRing(M, 'a, b')
        sage: M(c/d + 1/c)  # indirect doctest
        c*d^-1 + c^-1
        sage: N(a + b/c/d + 1/b)  # indirect doctest
        a + (c^-1*d^-1)*b + b^-1
    """
    vars_P = P.variable_names()
    vars_M = M.variable_names()
    if not set(vars_M) & set(vars_P):
        raise TypeError('no common variables')

    def index(T, value):
        try:
            return T.index(value)
        except ValueError:
            return None

    def value(d, R):
        assert d
        if len(d) == 1:
            k, v = next(iter(d.items()))
            if all(i == 0 for i in k):
                return R(v)
        return R(M(d))

    group_by = tuple(index(vars_M, var) for var in vars_P)
    indices = list(range(len(vars_M)))
    for g in group_by:
        if g is not None:
            indices[g] = None
    D = _split_dict_(d, indices, group_by)
    try:
        return {k: value(v, P.base_ring()) for k, v in D.items()}
    except (ValueError, TypeError):
        pass
    return sum(P({k: 1}) * value(v, P) for k, v in D.items()).dict()


class LaurentPolynomialRing_generic(CommutativeRing, Parent):
    """
    Laurent polynomial ring (base class).

    EXAMPLES:

    This base class inherits from :class:`~sage.rings.ring.CommutativeRing`.
    Since :trac:`11900`, it is also initialised as such::

        sage: R.<x1,x2> = LaurentPolynomialRing(QQ)
        sage: R.category()
        Join of Category of unique factorization domains and Category of commutative algebras over (number fields and quotient fields and metric spaces) and Category of infinite sets
        sage: TestSuite(R).run()

    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,2,'x')
            sage: R == loads(dumps(R))
            True
        """
        self._n = R.ngens()
        self._R = R
        names = R.variable_names()
        self._one_element = self.element_class(self, R.one())
        CommutativeRing.__init__(self, R.base_ring(), names=names,
                                 category=R.category())

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').ngens()
            2
            sage: LaurentPolynomialRing(QQ,1,'x').ngens()
            1
        """
        return self._n

    def gen(self, i=0):
        r"""
        Returns the `i^{th}` generator of self.  If i is not specified, then
        the first generator will be returned.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').gen()
            x0
            sage: LaurentPolynomialRing(QQ,2,'x').gen(0)
            x0
            sage: LaurentPolynomialRing(QQ,2,'x').gen(1)
            x1

        TESTS::

            sage: LaurentPolynomialRing(QQ,2,'x').gen(3)
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

        - ``depth`` -- an integer or :mod:`Infinity <sage.rings.infinity>`.

        OUTPUT:

        A tuple of strings.

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
        elif depth == 1:
            return self.variable_names()
        else:
            my_vars = self.variable_names()
            try:
               return self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
            except AttributeError:
                return my_vars


    def is_integral_domain(self, proof = True):
        """
        Returns True if self is an integral domain.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').is_integral_domain()
            True

        The following used to fail; see :trac:`7530`::

            sage: L = LaurentPolynomialRing(ZZ, 'X')
            sage: L['Y']
            Univariate Polynomial Ring in Y over Univariate Laurent Polynomial Ring in X over Integer Ring
        """
        return self.base_ring().is_integral_domain(proof)

    def is_noetherian(self):
        """
        Returns True if self is Noetherian.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').is_noetherian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def construction(self):
        """
        Return the construction of ``self``.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x,y').construction()
            (LaurentPolynomialFunctor,
            Univariate Laurent Polynomial Ring in x over Rational Field)

        """
        from sage.categories.pushout import LaurentPolynomialFunctor
        vars = self.variable_names()
        if len(vars) == 1:
            return LaurentPolynomialFunctor(vars[0], False), self.base_ring()
        else:
            return LaurentPolynomialFunctor(vars[-1], True), LaurentPolynomialRing(self.base_ring(), vars[:-1])

    def completion(self, p, prec=20, extras=None):
        """
        EXAMPLES::

            sage: P.<x>=LaurentPolynomialRing(QQ)
            sage: P
            Univariate Laurent Polynomial Ring in x over Rational Field
            sage: PP=P.completion(x)
            sage: PP
            Laurent Series Ring in x over Rational Field
            sage: f=1-1/x
            sage: PP(f)
            -x^-1 + 1
            sage: 1/PP(f)
            -x - x^2 - x^3 - x^4 - x^5 - x^6 - x^7 - x^8 - x^9 - x^10 - x^11 - x^12 - x^13 - x^14 - x^15 - x^16 - x^17 - x^18 - x^19 - x^20 + O(x^21)

        TESTS:

        Check that the precision is taken into account (:trac:`24431`)::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: L.completion('x', 100).default_prec()
            100
            sage: L.completion('x', 20).default_prec()
            20
        """
        if str(p) == self._names[0] and self._n == 1:
            from sage.rings.laurent_series_ring import LaurentSeriesRing
            R = self.polynomial_ring().completion(self._names[0], prec)
            return LaurentSeriesRing(R)
        else:
            raise TypeError("Cannot complete %s with respect to %s" % (self, p))

    def remove_var(self, var):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,'x,y,z')
            sage: R.remove_var('x')
            Multivariate Laurent Polynomial Ring in y, z over Rational Field
            sage: R.remove_var('x').remove_var('y')
            Univariate Laurent Polynomial Ring in z over Rational Field
        """
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
        different base rings works (:trac:`15345`)::

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
        if type(self) != type(right):
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

            sage: latex(LaurentPolynomialRing(QQ,2,'x'))
            \Bold{Q}[x_{0}^{\pm 1}, x_{1}^{\pm 1}]
        """
        vars = ', '.join(a + r'^{\pm 1}' for a in self.latex_variable_names())
        return "%s[%s]" % (latex(self.base_ring()), vars)

    def _ideal_class_(self, n=0):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x')._ideal_class_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # One may eventually want ideal classes in these guys.
        raise NotImplementedError

    def ideal(self, *args, **kwds):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').ideal([1])
            Ideal (1) of Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field

        TESTS:
 
        check that :trac:`26421` is fixed:

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: P.<x> = PolynomialRing(R)
            sage: p = x-t
            sage: p.content_ideal()    # indirect doctest
            Ideal (-t, 1) of Univariate Laurent Polynomial Ring in t over Integer Ring
        """
        return lp_ideal.LaurentPolynomialIdeal(self, *args, **kwds)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        EXAMPLES::

            sage: T.<t> = ZZ[]
            sage: K.<i> = NumberField(t^2 + 1)
            sage: L.<x,y> = LaurentPolynomialRing(K)
            sage: L._is_valid_homomorphism_(K, (K(1/2), K(3/2)))
            True
            sage: Q5 = Qp(5); i5 = Q5(-1).sqrt()
            sage: L._is_valid_homomorphism_(Q5, (Q5(1/2), Q5(3/2))) # no coercion
            False
            sage: L._is_valid_homomorphism_(Q5, (Q5(1/2), Q5(3/2)), base_map=K.hom([i5]))
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
        Returns the term order of self.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').term_order()
            Degree reverse lexicographic term order

        """
        return self._R.term_order()

    def is_finite(self):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').is_finite()
            False

        """
        return False

    def is_field(self, proof = True):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').is_field()
            False
        """
        return False

    def polynomial_ring(self):
        """
        Returns the polynomial ring associated with self.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').polynomial_ring()
            Multivariate Polynomial Ring in x0, x1 over Rational Field
            sage: LaurentPolynomialRing(QQ,1,'x').polynomial_ring()
            Multivariate Polynomial Ring in x over Rational Field
        """
        return self._R

    def characteristic(self):
        """
        Returns the characteristic of the base ring.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').characteristic()
            0
            sage: LaurentPolynomialRing(GF(3),2,'x').characteristic()
            3

        """
        return self.base_ring().characteristic()

    def krull_dimension(self):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').krull_dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def random_element(self, low_degree = -2, high_degree = 2, terms = 5, choose_degree=False,*args, **kwds):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_exact(self):
        """
        Returns True if the base ring is exact.

        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').is_exact()
            True
            sage: LaurentPolynomialRing(RDF,2,'x').is_exact()
            False
        """
        return self.base_ring().is_exact()

    def change_ring(self, base_ring=None, names=None, sparse=False, order=None):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,2,'x')
            sage: R.change_ring(ZZ)
            Multivariate Laurent Polynomial Ring in x0, x1 over Integer Ring

        Check that the distinction between a univariate ring and a multivariate ring with one
        generator is preserved::

            sage: P.<x> = LaurentPolynomialRing(QQ, 1)
            sage: P
            Multivariate Laurent Polynomial Ring in x over Rational Field
            sage: K.<i> = CyclotomicField(4)
            sage: P.change_ring(K)
            Multivariate Laurent Polynomial Ring in x over Cyclotomic Field of order 4 and degree 2
        """
        if base_ring is None:
            base_ring = self.base_ring()
        if names is None:
            names = self.variable_names()
        if isinstance(self, LaurentPolynomialRing_univariate):
            return LaurentPolynomialRing(base_ring, names[0], sparse = sparse)

        if order is None:
            order = self.polynomial_ring().term_order()
        return LaurentPolynomialRing(base_ring, self._n, names, order = order)

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

class LaurentPolynomialRing_univariate(LaurentPolynomialRing_generic):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,'x')
            sage: type(L)
            <class 'sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_univariate_with_category'>
            sage: L == loads(dumps(L))
            True


        TESTS::

            sage: TestSuite(LaurentPolynomialRing(Zmod(4), 'y')).run()
            sage: TestSuite(LaurentPolynomialRing(ZZ, 'u')).run()
            sage: TestSuite(LaurentPolynomialRing(Zmod(4)['T'], 'u')).run()
        """
        if R.ngens() != 1:
            raise ValueError("must be 1 generator")
        LaurentPolynomialRing_generic.__init__(self, R)

    Element = LaurentPolynomial_univariate

    def _repr_(self):
        """
        TESTS::

            sage: LaurentPolynomialRing(QQ,'x')  # indirect doctest
            Univariate Laurent Polynomial Ring in x over Rational Field
        """
        return "Univariate Laurent Polynomial Ring in %s over %s"%(self._R.variable_name(), self._R.base_ring())

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: L(1/2)
            1/2

            sage: L(x + 3/x)
            3*x^-1 + x

        ::

            sage: L(exp(x))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^x to a rational

        ::

            sage: U = LaurentPolynomialRing(QQ, 'a')
            sage: V = LaurentPolynomialRing(QQ, 'c')
            sage: L.<a, b, c, d> = LaurentPolynomialRing(QQ)
            sage: M = LaurentPolynomialRing(QQ, 'c, d')
            sage: Mc, Md = M.gens()
            sage: N = LaurentPolynomialRing(M, 'a, b')
            sage: Na, Nb = N.gens()
            sage: U(Na)
            a
            sage: V(Mc)
            c

            sage: M(L(0))
            0
            sage: N(L(0))
            0
            sage: L(M(0))
            0
            sage: L(N(0))
            0

        ::

            sage: A.<a> = LaurentPolynomialRing(QQ)
            sage: B.<b> = LaurentPolynomialRing(A)
            sage: B(a)
            a
            sage: C.<c> = LaurentPolynomialRing(B)
            sage: B(C(b))
            b
            sage: D.<d, e> = LaurentPolynomialRing(B)
            sage: B(D(b))
            b

        TESTS:

        Check that conversion back from fraction field does work (:trac:`26425`)::

            sage: R.<t> = LaurentPolynomialRing(ZZ)
            sage: F = FractionField(R)
            sage: R(F(25/(5*t**2)))
            5*t^-2
            sage: R(F(1/(1+t**2)))
            Traceback (most recent call last):
            ...
            TypeError: fraction must have unit denominator
        """
        from sage.structure.element import Expression
        from sage.rings.fraction_field_element import FractionFieldElement
        if isinstance(x, Expression):
            return x.laurent_polynomial(ring=self)

        elif isinstance(x, (LaurentPolynomial_univariate, LaurentPolynomial_mpair)):
            P = x.parent()
            if set(self.variable_names()) & set(P.variable_names()):
                if isinstance(x, LaurentPolynomial_univariate):
                    d = {(k,): v for k, v in x.dict().items()}
                else:
                    d = x.dict()
                x = _split_laurent_polynomial_dict_(self, P, d)
                x = {k[0]: v for k, v in x.items()}
            elif P is self.base_ring():
                x = {0: x}
            elif x.is_constant() and self.has_coerce_map_from(x.parent().base_ring()):
                return self(x.constant_coefficient())
            elif len(self.variable_names()) == len(P.variable_names()):
                x = x.dict()

        elif isinstance(x, FractionFieldElement):
            # since the field of fraction of self is defined corresponding to the polynomial ring of self
            # the conversion of its elements back must be treated separately (:trac:`26425`).
            P = x.parent()
            d = self(x.denominator())
            if not d.is_unit():
                raise TypeError("fraction must have unit denominator")
            return self(x.numerator()) * d.inverse_of_unit()

        return self.element_class(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: loads(dumps(L)) == L
            True
        """
        return LaurentPolynomialRing_univariate, (self._R,)

class LaurentPolynomialRing_mpair(LaurentPolynomialRing_generic):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,2,'x')
            sage: type(L)
            <class
            'sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_mpair_with_category'>
            sage: L == loads(dumps(L))
            True
        """
        if R.ngens() <= 0:
            raise ValueError("n must be positive")
        if not R.base_ring().is_integral_domain():
            raise ValueError("base ring must be an integral domain")
        LaurentPolynomialRing_generic.__init__(self, R)

    Element = LaurentPolynomial_mpair

    def _repr_(self):
        """
        TESTS::

            sage: LaurentPolynomialRing(QQ,2,'x').__repr__()
            'Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field'
            sage: LaurentPolynomialRing(QQ,1,'x').__repr__()
            'Multivariate Laurent Polynomial Ring in x over Rational Field'
        """
        return "Multivariate Laurent Polynomial Ring in %s over %s"%(", ".join(self._R.variable_names()), self._R.base_ring())

    def monomial(self, *args):
        r"""
        Return the monomial whose exponents are given in argument.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x', 2)
            sage: L.monomial(-3, 5)
            x0^-3*x1^5
            sage: L.monomial(1, 1)
            x0*x1
            sage: L.monomial(0, 0)
            1
            sage: L.monomial(-2, -3)
            x0^-2*x1^-3

            sage: x0, x1 = L.gens()
            sage: L.monomial(-1, 2) == x0^-1 * x1^2
            True

            sage: L.monomial(1, 2, 3)
            Traceback (most recent call last):
            ...
            TypeError: tuple key must have same length as ngens
        """
        if len(args) != self.ngens():
            raise TypeError("tuple key must have same length as ngens")

        from sage.rings.polynomial.polydict import ETuple
        m = ETuple(args, int(self.ngens()))
        return self.element_class(self, self.polynomial_ring().one(), m)

    def _element_constructor_(self, x, mon=None):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,2,'x')
            sage: L(1/2)
            1/2

            sage: M = LaurentPolynomialRing(QQ, 'x, y')
            sage: var('x, y')
            (x, y)
            sage: M(x/y + 3/x)
            x*y^-1 + 3*x^-1

        ::

            sage: M(exp(x))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^x to a rational

        ::

            sage: L.<a, b, c, d> = LaurentPolynomialRing(QQ)
            sage: M = LaurentPolynomialRing(QQ, 'c, d')
            sage: Mc, Md = M.gens()
            sage: N = LaurentPolynomialRing(M, 'a, b')
            sage: Na, Nb = N.gens()
            sage: M(c/d)
            c*d^-1
            sage: N(a*b/c/d)
            (c^-1*d^-1)*a*b
            sage: N(c/d)
            c*d^-1
            sage: L(Mc)
            c
            sage: L(Nb)
            b

            sage: M(L(0))
            0
            sage: N(L(0))
            0
            sage: L(M(0))
            0
            sage: L(N(0))
            0

            sage: U = LaurentPolynomialRing(QQ, 'a')
            sage: Ua = U.gen()
            sage: V = LaurentPolynomialRing(QQ, 'c')
            sage: Vc = V.gen()
            sage: L(Ua)
            a
            sage: L(Vc)
            c
            sage: N(Ua)
            a
            sage: M(Vc)
            c

            sage: P = LaurentPolynomialRing(QQ, 'a, b')
            sage: Q = LaurentPolynomialRing(P, 'c, d')
            sage: Q(P.0)
            a

        ::

            sage: A.<a> = LaurentPolynomialRing(QQ)
            sage: B.<b> = LaurentPolynomialRing(A)
            sage: C = LaurentPolynomialRing(QQ, 'a, b')
            sage: C(B({1: a}))
            a*b
            sage: D.<d, e> = LaurentPolynomialRing(B)
            sage: F.<f, g> = LaurentPolynomialRing(D)
            sage: D(F(d*e))
            d*e

        ::

            sage: from sage.rings.polynomial.polydict import ETuple
            sage: R.<x,y,z> = LaurentPolynomialRing(QQ)
            sage: mon = ETuple({}, int(3))
            sage: P = R.polynomial_ring()
            sage: R(sum(P.gens()), mon)
            x + y + z
            sage: R(sum(P.gens()), (-1,-1,-1))
            y^-1*z^-1 + x^-1*z^-1 + x^-1*y^-1
        """
        from sage.structure.element import Expression

        if mon is not None:
            return self.element_class(self, x, mon)

        P = parent(x)
        if P is self.polynomial_ring():
            from sage.rings.polynomial.polydict import ETuple
            return self.element_class( self, x, mon=ETuple({}, int(self.ngens())) )

        elif isinstance(x, Expression):
            return x.laurent_polynomial(ring=self)

        elif isinstance(x, (LaurentPolynomial_univariate, LaurentPolynomial_mpair)):
            if self.variable_names() == P.variable_names():
                # No special processing needed here;
                #   handled by LaurentPolynomial_mpair.__init__
                pass
            elif set(self.variable_names()) & set(P.variable_names()):
                if isinstance(x, LaurentPolynomial_univariate):
                    d = {(k,): v for k, v in x.dict().items()}
                else:
                    d = x.dict()
                x = _split_laurent_polynomial_dict_(self, P, d)
            elif P is self.base_ring():
                from sage.rings.polynomial.polydict import ETuple
                mz = ETuple({}, int(self.ngens()))
                return self.element_class(self, {mz: x}, mz)
            elif x.is_constant() and self.has_coerce_map_from(P.base_ring()):
                return self(x.constant_coefficient())
            elif len(self.variable_names()) == len(P.variable_names()):
                x = x.dict()

        return self.element_class(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 2, 'x')
            sage: loads(dumps(L)) == L
            True
        """
        return LaurentPolynomialRing_mpair, (self._R,)


