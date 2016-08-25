"""
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
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>,
#                          William Stein <wstein@gmail.com>,
#                          Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.category_object import normalize_names
from sage.structure.element import is_Element
from sage.rings.ring import is_Ring
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import _single_variate as _single_variate_poly
from sage.rings.polynomial.polynomial_ring_constructor import _multi_variate as _multi_variate_poly
from sage.misc.latex import latex
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_mpair, LaurentPolynomial_univariate
from sage.rings.ring import CommutativeRing
from sage.structure.parent_gens import ParentWithGens

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

def LaurentPolynomialRing(base_ring, arg1=None, arg2=None, sparse = False, order='degrevlex', names = None, name=None):
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
           TypeError: You must specify the names of the variables.

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
    if is_Element(arg1) and not isinstance(arg1, (int, long, Integer)):
        arg1 = repr(arg1)
    if is_Element(arg2) and not isinstance(arg2, (int, long, Integer)):
        arg2 = repr(arg2)

    if isinstance(arg1, (int, long, Integer)):
        arg1, arg2 = arg2, arg1

    if not names is None:
        arg1 = names
    elif not name is None:
        arg1 = name

    if not is_Ring(base_ring):
        raise TypeError('base_ring must be a ring')

    if arg1 is None:
        raise TypeError("You must specify the names of the variables.")

    R = None
    if isinstance(arg1, (list, tuple)):
        arg1 = [str(x) for x in arg1]
    if isinstance(arg2, (list, tuple)):
        arg2 = [str(x) for x in arg2]
    if isinstance(arg2, (int, long, Integer)):
        # 3. LaurentPolynomialRing(base_ring, names, n, order='degrevlex'):
        if not isinstance(arg1, (list, tuple, str)):
            raise TypeError("You *must* specify the names of the variables.")
        n = int(arg2)
        names = arg1
        R = _multi_variate(base_ring, names, n, sparse, order)

    elif isinstance(arg1, str) or (isinstance(arg1, (list,tuple)) and len(arg1) == 1) and isinstance(arg1[0], str):
        if isinstance(arg1, (list,tuple)):
            arg1 = arg1[0]
        if not ',' in arg1:
            # 1. LaurentPolynomialRing(base_ring, name, sparse=False):
            if not arg2 is None:
                raise TypeError("if second arguments is a string with no commas, then there must be no other non-optional arguments")
            name = arg1
            R = _single_variate(base_ring, name, sparse)
        else:
            # 2-4. LaurentPolynomialRing(base_ring, names, order='degrevlex'):
            if not arg2 is None:
                raise TypeError("invalid input to LaurentPolynomialRing function; please see the docstring for that function")
            names = arg1.split(',')
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)
    elif isinstance(arg1, (list, tuple)):
        # LaurentPolynomialRing(base_ring, names (list or tuple), order='degrevlex'):
        names = arg1
        n = len(names)
        R = _multi_variate(base_ring, names, n, sparse, order)

    if arg1 is None and arg2 is None:
        raise TypeError("you *must* specify the indeterminates (as not None).")
    if R is None:
        raise TypeError("invalid input (%s, %s, %s) to PolynomialRing function; please see the docstring for that function"%(base_ring, arg1, arg2))

    return R

_cache = {}
def _get_from_cache(key):
    """
    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import _get_from_cache
        sage: L = LaurentPolynomialRing(QQ,2,'x')
        sage: L2 = _get_from_cache( (QQ,('x0','x1'),2,False,TermOrder('degrevlex')) ); L2
        Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field
        sage: L is L2
        True
    """
    try:
        if key in _cache:
            return _cache[key]   # put () here to re-enable weakrefs
    except TypeError as msg:
        raise TypeError('key = %s\n%s'%(key,msg))
    return None

def _save_in_cache(key, R):
    """
    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import _save_in_cache, _get_from_cache
        sage: L = LaurentPolynomialRing(QQ,2,'x')
        sage: _save_in_cache('testkey', L)
        sage: _get_from_cache('testkey')
        Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field
        sage: _ is L
        True
    """
    try:
        # We disable weakrefs since they cause segfault at the end of doctesting.
        #weakref.ref(R)
        _cache[key] = R
    except TypeError as msg:
        raise TypeError('key = %s\n%s'%(key,msg))

def _single_variate(base_ring, names, sparse):
    """
    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import _single_variate
        sage: _single_variate(QQ, ('x',), False)
        Univariate Laurent Polynomial Ring in x over Rational Field
    """
    names = normalize_names(1, names)
    key = (base_ring, names, sparse)
    P = _get_from_cache(key)
    if P is not None:
        return P
    prepend_string = "qk"
    while True:
        if prepend_string in names:
            prepend_string += 'k'
        else:
            break
    R = _single_variate_poly(base_ring, names, sparse, None)
    P = LaurentPolynomialRing_univariate(R, names)
    _save_in_cache(key, P)
    return P

def _multi_variate(base_ring, names, n, sparse, order):
    """
    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import _multi_variate
        sage: _multi_variate(QQ, ('x','y'), 2, False, 'degrevlex')
        Multivariate Laurent Polynomial Ring in x, y over Rational Field
    """
    # We need to come up with a name for the inverse that is easy to search
    # for in a string *and* doesn't overlap with the name that we already have.
    # For now, I'm going to use a name mangling with checking method.
    names = normalize_names(n, names)

    from .term_order import TermOrder
    order = TermOrder(order, n)

    if isinstance(names, list):
        names = tuple(names)
    elif isinstance(names, str):
        if ',' in names:
            names = tuple(names.split(','))

    key = (base_ring, names, n, sparse, order)
    P = _get_from_cache(key)
    if P is not None:
        return P
    prepend_string = "qk"
    while True:
        for a in names:
            if prepend_string in a:
                prepend_string += 'k'
                break
        else:
            break
    R = _multi_variate_poly(base_ring, names, n, sparse, order, None)
    P = LaurentPolynomialRing_mpair(R, prepend_string, names)
    _save_in_cache(key, P)
    return P

class LaurentPolynomialRing_generic(CommutativeRing, ParentWithGens):
    """
    Laurent polynomial ring (base class).

    EXAMPLES:

    This base class inherits from :class:`~sage.rings.ring.CommutativeRing`.
    Since :trac:`11900`, it is also initialised as such::

        sage: R.<x1,x2> = LaurentPolynomialRing(QQ)
        sage: R.category()
        Category of commutative rings
        sage: TestSuite(R).run()

    """
    def __init__(self, R, prepend_string, names):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,2,'x')
            sage: R == loads(dumps(R))
            True
        """
        self._n = R.ngens()
        self._R = R
        self._prepend_string = prepend_string
        CommutativeRing.__init__(self, R.base_ring(), names=names)
        self._populate_coercion_lists_(element_constructor=self._element_constructor_,
                                       init_no_parent=True)


    def __repr__(self):
        """
        TESTS::

            sage: LaurentPolynomialRing(QQ,2,'x').__repr__()
            'Multivariate Laurent Polynomial Ring in x0, x1 over Rational Field'
            sage: LaurentPolynomialRing(QQ,1,'x').__repr__()
            'Univariate Laurent Polynomial Ring in x over Rational Field'
        """
        if self._n == 1:
            return "Univariate Laurent Polynomial Ring in %s over %s"%(self._R.variable_name(), self._R.base_ring())
        else:
            return "Multivariate Laurent Polynomial Ring in %s over %s"%(", ".join(self._R.variable_names()), self._R.base_ring())

    def ngens(self):
        """
        Returns the number of generators of self.

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
        Returns the construction of self.

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
        """
        if str(p) == self._names[0] and self._n == 1:
            from sage.rings.laurent_series_ring import LaurentSeriesRing
            return LaurentSeriesRing(self.base_ring(), name=self._names[0])
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
            Composite map:
              From: Rational Field
              To:   Multivariate Laurent Polynomial Ring in x, y over Rational Field
              Defn:   Polynomial base injection morphism:
                      From: Rational Field
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                    then
                      Call morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Laurent Polynomial Ring in x, y over Rational Field

        Let us check that coercion between Laurent Polynomials over
        different base rings works (:trac:`15345`)::

            sage: R = LaurentPolynomialRing(ZZ, 'x')
            sage: T = LaurentPolynomialRing(QQ, 'x')
            sage: R.gen() + 3*T.gen()
            4*x
        """
        if R is self._R or (isinstance(R, LaurentPolynomialRing_generic)
            and self._R.has_coerce_map_from(R._R)):
            from sage.structure.coerce_maps import CallableConvertMap
            return CallableConvertMap(R, self, self._element_constructor_,
                                      parent_as_first_arg=False)
        elif isinstance(R, LaurentPolynomialRing_generic) and \
             R.variable_names() == self.variable_names() and \
             self.base_ring().has_coerce_map_from(R.base_ring()):
            return True

        f = self._R.coerce_map_from(R)
        if f is not None:
            from sage.categories.homset import Hom
            from sage.categories.morphism import CallMorphism
            return CallMorphism(Hom(self._R, self)) * f

    def __cmp__(left, right):
        """
        EXAMPLES::

            sage: R = LaurentPolynomialRing(QQ,'x,y,z')
            sage: P = LaurentPolynomialRing(ZZ,'x,y,z')
            sage: Q = LaurentPolynomialRing(QQ,'x,y')

            sage: cmp(R,R)
            0
            sage: cmp(R,Q) == 0
            False
            sage: cmp(Q,P) == 0
            False
            sage: cmp(R,P) == 0
            False
        """
        c = cmp(type(left), type(right))
        if c == 0:
            c = cmp(left._R, right._R)
        return c

    def _latex_(self):
        """
        EXAMPLES::

            sage: latex(LaurentPolynomialRing(QQ,2,'x'))
            \Bold{Q}[x_{0}^{\pm 1}, x_{1}^{\pm 1}]
        """
        vars = ', '.join([a + '^{\pm 1}' for a in self.latex_variable_names()])
        return "%s[%s]"%(latex(self.base_ring()), vars)

    def _ideal_class_(self, n=0):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x')._ideal_class_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # One may eventually want ideals in these guys.
        raise NotImplementedError

    def ideal(self):
        """
        EXAMPLES::

            sage: LaurentPolynomialRing(QQ,2,'x').ideal()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: L.<x,y> = LaurentPolynomialRing(QQ)
            sage: L._is_valid_homomorphism_(QQ, (1/2, 3/2))
            True
        """
        if not codomain.has_coerce_map_from(self.base_ring()):
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
        """
        if base_ring is None:
            base_ring = self.base_ring()
        if names is None:
            names = self.variable_names()
        if self._n == 1:
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
    def __init__(self, R, names):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,'x')
            sage: type(L)
            <class 'sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_univariate_with_category'>
            sage: L == loads(dumps(L))
            True
        """
        if R.ngens() != 1:
            raise ValueError("must be 1 generator")
        if not R.base_ring().is_integral_domain():
            raise ValueError("base ring must be an integral domain")
        LaurentPolynomialRing_generic.__init__(self, R, '', names)

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: L(1/2)
            1/2
        """
        return LaurentPolynomial_univariate(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ, 'x')
            sage: loads(dumps(L)) == L
            True
        """
        return LaurentPolynomialRing_univariate, (self._R, self._names)

class LaurentPolynomialRing_mpair(LaurentPolynomialRing_generic):
    def __init__(self, R, prepend_string, names):
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
        LaurentPolynomialRing_generic.__init__(self, R, prepend_string, names)

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,2,'x')
            sage: L(1/2)
            1/2
        """
        return LaurentPolynomial_mpair(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = LaurentPolynomialRing(QQ,2,'x')
            sage: loads(dumps(L)) == L
            True
        """
        return LaurentPolynomialRing_mpair, (self._R, self._prepend_string, self._names)


