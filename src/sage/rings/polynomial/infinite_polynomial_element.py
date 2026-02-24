r"""
Elements of Infinite Polynomial Rings

AUTHORS:

- Simon King <simon.king@nuigalway.ie>
- Mike Hansen <mhansen@gmail.com>

An Infinite Polynomial Ring has generators `x_\ast, y_\ast,...`, so
that the variables are of the form `x_0, x_1, x_2, ..., y_0, y_1,
y_2,...,...` (see :mod:`~sage.rings.polynomial.infinite_polynomial_ring`).
Using the generators, we can create elements as follows::

    sage: X.<x,y> = InfinitePolynomialRing(QQ)
    sage: a = x[3]
    sage: b = y[4]
    sage: a
    x_3
    sage: b
    y_4
    sage: c = a*b + a^3 - 2*b^4
    sage: c
    x_3^3 + x_3*y_4 - 2*y_4^4

Any Infinite Polynomial Ring ``X`` is equipped with a monomial ordering.
We only consider monomial orderings in which:

    ``X.gen(i)[m] > X.gen(j)[n]`` `\iff` ``i<j``, or ``i==j`` and ``m>n``

Under this restriction, the monomial ordering can be lexicographic
(default), degree lexicographic, or degree reverse lexicographic.
Here, the ordering is lexicographic, and elements can be compared
as usual::

    sage: X._order
    'lex'
    sage: a > b
    True

Note that, when a method is called that is not directly implemented
for 'InfinitePolynomial', it is tried to call this method for the
underlying *classical* polynomial. This holds, e.g., when applying the
``latex`` function::

    sage: latex(c)
    x_{3}^{3} + x_{3} y_{4} - 2 y_{4}^{4}

There is a permutation action on Infinite Polynomial Rings by
permuting the indices of the variables::

    sage: P = Permutation(((4,5),(2,3)))
    sage: c^P
    x_2^3 + x_2*y_5 - 2*y_5^4

Note that ``P(0)==0``, and thus variables of index zero are invariant
under the permutation action.  More generally, if ``P`` is any
callable object that accepts nonnegative integers as input and
returns nonnegative integers, then ``c^P`` means to apply ``P`` to
the variable indices occurring in ``c``.

If you want to substitute variables you can use the standard polynomial
methods, such as
:meth:`~sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_sparse.subs`::

    sage: R.<x,y> = InfinitePolynomialRing(QQ)
    sage: f = x[1] + x[1]*x[2]*x[3]
    sage: f.subs({x[1]: x[0]})
    x_3*x_2*x_0 + x_0
    sage: g = x[0] + x[1] + y[0]
    sage: g.subs({x[0]: y[0]})
    x_1 + 2*y_0

TESTS:

We test whether coercion works, even in complicated cases in which
finite polynomial rings are merged with infinite polynomial rings::

    sage: A.<a> = InfinitePolynomialRing(ZZ, implementation='sparse', order='degrevlex')
    sage: B.<b_2,b_1> = A[]
    sage: C.<b,c> = InfinitePolynomialRing(B, order='degrevlex')
    sage: C
    Infinite polynomial ring in b, c over Infinite polynomial ring in a over Integer Ring
    sage: 1/2*b_1*a[4] + c[3]
    1/2*a_4*b_1 + c_3
"""

# ****************************************************************************
#       Copyright (C) 2009 Simon King <king@mathematik.nuigalway.ie>
#                          and Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import copy

from sage.misc.cachefunc import cached_method
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.commutative_polynomial import CommutativePolynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.structure.element import RingElement, parent
from sage.structure.factorization import Factorization
from sage.structure.richcmp import richcmp


class InfinitePolynomial(CommutativePolynomial,
                         metaclass=InheritComparisonClasscallMetaclass):
    """
    Create an element of a Polynomial Ring with a Countably Infinite Number of Variables.

    Usually, an InfinitePolynomial is obtained by using the generators
    of an Infinite Polynomial Ring (see :mod:`~sage.rings.polynomial.infinite_polynomial_ring`)
    or by conversion.

    INPUT:

    - ``A`` -- an Infinite Polynomial Ring
    - ``p`` -- a *classical* polynomial that can be interpreted in ``A``

    ASSUMPTIONS:

    In the dense implementation, it must be ensured that the argument
    ``p`` coerces into ``A._P`` by a name preserving conversion map.

    In the sparse implementation, in the direct construction of an
    infinite polynomial, it is *not* tested whether the argument ``p``
    makes sense in ``A``.

    EXAMPLES::

        sage: from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        sage: X.<alpha> = InfinitePolynomialRing(ZZ)
        sage: P.<alpha_1,alpha_2> = ZZ[]

    Currently, ``P`` and ``X._P`` (the underlying polynomial ring of
    ``X``) both have two variables::

        sage: X._P
        Multivariate Polynomial Ring in alpha_1, alpha_0 over Integer Ring

    By default, a coercion from ``P`` to  ``X._P`` would not be name preserving.
    However, this is taken care for; a name preserving conversion is impossible,
    and by consequence an error is raised::

        sage: InfinitePolynomial(X, (alpha_1+alpha_2)^2)
        Traceback (most recent call last):
        ...
        TypeError: Could not find a mapping of the passed element to this ring.

    When extending the underlying polynomial ring, the construction of
    an infinite polynomial works::

        sage: alpha[2]
        alpha_2
        sage: InfinitePolynomial(X, (alpha_1+alpha_2)^2)
        alpha_2^2 + 2*alpha_2*alpha_1 + alpha_1^2

    In the sparse implementation, it is not checked whether the
    polynomial really belongs to the parent, and when it does not,
    the results may be unexpected due to coercions::

        sage: Y.<alpha,beta> = InfinitePolynomialRing(GF(2), implementation='sparse')
        sage: a = (alpha_1+alpha_2)^2
        sage: InfinitePolynomial(Y, a)
        alpha_0^2 + beta_0^2

    However, it is checked when doing a conversion::

        sage: Y(a)
        alpha_2^2 + alpha_1^2
    """

    @staticmethod
    def __classcall_private__(cls, A, p):
        r"""
        TESTS::

            sage: from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
            sage: X.<x,y> = InfinitePolynomialRing(ZZ, implementation='sparse')
            sage: xy = (x[0] + y[0]).polynomial()
            sage: xy.parent()
            Multivariate Polynomial Ring in x_1, x_0, y_1, y_0 over Integer Ring
            sage: sparse_xy = InfinitePolynomial(X, xy); sparse_xy
            x_0 + y_0
            sage: isinstance(sparse_xy, InfinitePolynomial)
            True
            sage: type(sparse_xy)
            <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_sparse'>
            sage: X.<x,y> = InfinitePolynomialRing(ZZ, implementation='dense')
            sage: dense_xy = InfinitePolynomial(X, xy); dense_xy
            x_0 + y_0
            sage: isinstance(dense_xy, InfinitePolynomial)
            True
            sage: type(dense_xy)
            <class 'sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_dense'>
        """
        from sage.structure.element import parent
        if hasattr(A, '_P'):
            if parent(p) is A._P or (A._P.base_ring().has_coerce_map_from(parent(p))):
                return InfinitePolynomial_dense(A, p)
            # MPolynomialRing_polydict is crab. So, in that case, use sage_eval
            from sage.rings.polynomial.multi_polynomial_ring import (
                MPolynomialRing_polydict,
            )
            if isinstance(A._P, MPolynomialRing_polydict):
                from sage.misc.sage_eval import sage_eval
                from sage.rings.polynomial.infinite_polynomial_ring import (
                    GenDictWithBasering,
                )
                p = sage_eval(repr(p), GenDictWithBasering(A._P, A._P.gens_dict()))
                return InfinitePolynomial_dense(A, p)
            else:
                # Now there remains to fight the oddities and bugs of libsingular.
                PP = p.parent()
                if A._P.has_coerce_map_from(PP):
                    if A._P.ngens() == PP.ngens():  # coercion is sometimes by position!
                        f = PP.hom(PP.variable_names(), A._P)
                        try:
                            return InfinitePolynomial_dense(A, f(p))
                        except (ValueError, TypeError):
                            # last desperate attempt: String conversion
                            from sage.misc.sage_eval import sage_eval
                            from sage.rings.polynomial.infinite_polynomial_ring import (
                                GenDictWithBasering,
                            )
                            # the base ring may be a function field, therefore
                            # we need GenDictWithBasering
                            return InfinitePolynomial_dense(A, sage_eval(repr(p), GenDictWithBasering(A._P, A._P.gens_dict())))
                    return InfinitePolynomial_dense(A, A._P(p))
                # there is no coercion, so, we set up a name-preserving map.
                SV = set(repr(x) for x in p.variables())
                f = PP.hom([x if x in SV else 0 for x in PP.variable_names()], A._P)
                try:
                    return InfinitePolynomial_dense(A, f(p))
                except (ValueError, TypeError):
                    # last desperate attempt: String conversion
                    from sage.misc.sage_eval import sage_eval
                    from sage.rings.polynomial.infinite_polynomial_ring import (
                        GenDictWithBasering,
                    )
                    # the base ring may be a function field, therefore
                    # we need GenDictWithBasering
                    return InfinitePolynomial_dense(A, sage_eval(repr(p), GenDictWithBasering(A._P, A._P.gens_dict())))
        return InfinitePolynomial_sparse(A, p)

    # Construction and other basic methods
    # We assume that p is good input. Type checking etc. is now done
    # in the _element_constructor_ of the parent.
    def __init__(self, A, p):
        """
        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[1] + x[2]
            sage: a == loads(dumps(a))
            True
        """
        # Despite the above comment, it can still happen that p is in
        # the wrong ring and we get here without going through
        # _element_constructor_.  See trac 22514 for examples.
        # So a little extra checking is done here.
        if not isinstance(p, MPolynomial) or p.base_ring() is not A.base_ring():
            # coerce to a convenient multivariate polynomial ring
            p = A._minP(p)

        self._has_footprint = False
        self._footprint = {}
        self._p = p
        RingElement.__init__(self, A)

    def _repr_(self):
        """
        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: str(x[1] + x[2])  # indirect doctest
            'x_2 + x_1'
        """
        return repr(self._p)

    def __hash__(self):
        """
        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1]
            sage: b = 1 + 4*x[1]
            sage: hash(a) != hash(b)
            True
        """
        return hash(self._p)

    def polynomial(self):
        """
        Return the underlying polynomial.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(GF(7))
            sage: p = x[2]*y[1] + 3*y[0]
            sage: p
            x_2*y_1 + 3*y_0
            sage: p.polynomial()
            x_2*y_1 + 3*y_0
            sage: p.polynomial().parent()
            Multivariate Polynomial Ring in x_2, x_1, x_0, y_2, y_1, y_0
             over Finite Field of size 7
            sage: p.parent()
            Infinite polynomial ring in x, y over Finite Field of size 7
        """
        return self._p

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: X.<alpha> = InfinitePolynomialRing(QQ)
            sage: latex(alpha[3]*alpha[2]^2)  # indirect doctest
            \alpha_{3} \alpha_{2}^{2}
        """
        return self._p._latex_()

    def subs(self, fixed=None, **kwargs):
        """
        Substitute variables in ``self``.

        INPUT:

        - ``fixed`` -- (optional) ``dict`` with ``{variable: value}`` pairs
        - ``**kwargs`` -- named parameters

        OUTPUT: the resulting substitution

        EXAMPLES::

            sage: R.<x,y> = InfinitePolynomialRing(QQ)
            sage: f = x[1] + x[1]*x[2]*x[3]

        Passing ``fixed={x[1]: x[0]}``. Note that the keys may be given
        using the generators of the infinite polynomial ring
        or as a string::

            sage: f.subs({x[1]: x[0]})
            x_3*x_2*x_0 + x_0
            sage: f.subs({'x_1': x[0]})
            x_3*x_2*x_0 + x_0

        Passing the variables as names parameters::

            sage: f.subs(x_1=y[1])
            x_3*x_2*y_1 + y_1
            sage: f.subs(x_1=y[1], x_2=2)
            2*x_3*y_1 + y_1

        The substitution returns the original polynomial if you try
        to substitute a variable not present::

            sage: g = x[0] + x[1]
            sage: g.subs({y[0]: x[0]})
            x_1 + x_0

        The substitution can also handle matrices::

            sage: # needs sage.modules
            sage: M = matrix([[1,0], [0,2]])
            sage: N = matrix([[0,3], [4,0]])
            sage: g = x[0]^2 + 3*x[1]
            sage: g.subs({'x_0': M})
            [3*x_1 + 1         0]
            [        0 3*x_1 + 4]
            sage: g.subs({x[0]: M, x[1]: N})
            [ 1  9]
            [12  4]

        If you pass both ``fixed`` and ``kwargs``, any conflicts
        will defer to ``fixed``::

            sage: R.<x,y> = InfinitePolynomialRing(QQ)
            sage: f = x[0]
            sage: f.subs({x[0]: 1})
            1
            sage: f.subs(x_0=5)
            5
            sage: f.subs({x[0]: 1}, x_0=5)
            1

        TESTS::

            sage: # needs sage.modules
            sage: g.subs(fixed=x[0], x_1=N)
            Traceback (most recent call last):
            ...
            ValueError: fixed must be a dict
        """
        if fixed:
            if not isinstance(fixed, dict):
                raise ValueError('fixed must be a dict')
            kwargs.update(fixed)
        try:
            return self(**kwargs)
        except TypeError:
            str_kwargs = {str(k): v for k, v in kwargs.items()}
            return self(**str_kwargs)

    def ring(self):
        """
        The ring which ``self`` belongs to.

        This is the same as ``self.parent()``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(ZZ, implementation='sparse')
            sage: p = x[100]*y[1]^3*x[1]^2 + 2*x[10]*y[30]
            sage: p.ring()
            Infinite polynomial ring in x, y over Integer Ring
        """
        return self.parent()

    def is_unit(self):
        r"""
        Answer whether ``self`` is a unit.

        EXAMPLES::

            sage: R1.<x,y> = InfinitePolynomialRing(ZZ)
            sage: R2.<a,b> = InfinitePolynomialRing(QQ)
            sage: (1 + x[2]).is_unit()
            False
            sage: R1(1).is_unit()
            True
            sage: R1(2).is_unit()
            False
            sage: R2(2).is_unit()
            True
            sage: (1 + a[2]).is_unit()
            False

        Check that :issue:`22454` is fixed::

            sage: _.<x> = InfinitePolynomialRing(Zmod(4))
            sage: (1 + 2*x[0]).is_unit()
            True
            sage: (x[0]*x[1]).is_unit()
            False
            sage: _.<x> = InfinitePolynomialRing(Zmod(900))
            sage: (7+150*x[0] + 30*x[1] + 120*x[1]*x[100]).is_unit()
            True

        TESTS::

            sage: R.<x> = InfinitePolynomialRing(ZZ.quotient_ring(8))
            sage: [R(i).is_unit() for i in range(8)]
            [False, True, False, True, False, True, False, True]
        """
        return self._p.is_unit()

    def is_nilpotent(self):
        r"""
        Return ``True`` if ``self`` is nilpotent, i.e., some power of ``self``
        is 0.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQbar)                                 # needs sage.rings.number_field
            sage: (x[0] + x[1]).is_nilpotent()                                          # needs sage.rings.number_field
            False
            sage: R(0).is_nilpotent()                                                   # needs sage.rings.number_field
            True
            sage: _.<x> = InfinitePolynomialRing(Zmod(4))
            sage: (2*x[0]).is_nilpotent()
            True
            sage: (2+x[4]*x[7]).is_nilpotent()
            False
            sage: _.<y> = InfinitePolynomialRing(Zmod(100))
            sage: (5+2*y[0] + 10*(y[0]^2+y[1]^2)).is_nilpotent()
            False
            sage: (10*y[2] + 20*y[5] - 30*y[2]*y[5] + 70*(y[2]^2+y[5]^2)).is_nilpotent()
            True
        """
        return self._p.is_nilpotent()

    def is_constant(self):
        """
        Return ``True`` if ``self`` is a constant and ``False`` otherwise.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQbar)
            sage: f = 3*x[3]^2 - 2*x[1] + 5
            sage: f.is_constant()
            False
            sage: g = 10*f^0
            sage: g.is_constant()
            True
        """
        return self._p.is_constant()

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate polynomial.

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: f = 3*x[3]^2 - 2*x[2]*y[2] + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x[3]^2
            sage: f.constant_coefficient()
            0
        """
        return self._p.constant_coefficient()

    def is_monomial(self):
        """
        Return whether ``self`` is a monomial.

        A monomial is a product of generators with coefficient 1.
        Use :meth:`is_term` to allow an arbitrary coefficient.

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[2]*y[2]
            sage: p.is_monomial()
            False
            sage: (p/2).is_monomial()
            True
        """
        return self._p.is_monomial()

    def is_term(self):
        """
        Return whether ``self`` is a term.

        A term is a product of generators times a non-zero element of
        the base ring.  Use :meth:`is_monomial` to check whether the
        element is a term with coefficient 1.

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[2]*y[2] + 5
            sage: p.is_term()
            False
            sage: (p-5).is_term()
            True
        """
        return self._p.is_term()

    def degree(self, x=None, std_grading=False):
        """
        Return the degree of ``self`` in ``x``, where ``x`` must
        be one of the generators for the parent of ``self``.

        INPUT:

        - ``x`` -- a generator of the parent of ``self``. If ``x`` is
          not specified (or is ``None``), return the total degree,
          which is the maximum degree of any monomial. Note that a
          weighted term ordering alters the grading of the generators
          of the ring; see the tests below.  To avoid this behavior,
          set the optional argument ``std_grading=True``.

        OUTPUT: integer

        EXAMPLES::

            sage: R.<x, y> = InfinitePolynomialRing(QQ)
            sage: p = x[3]*x[2]^2*y[5]^3 + 1
            sage: p.degree()
            6
            sage: p.degree(y[5])
            3
            sage: p.degree(y[3])
            0

        """
        if x is not None:
            x = x._p
        return self._p.degree(x=x, std_grading=std_grading)

    def numerator(self):
        r"""
        Return a numerator of ``self``, computed as ``self * self.denominator()``.

        .. WARNING::

           This is not the numerator of the rational function
           defined by ``self``, which would always be ``self`` since it is a
           polynomial.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: p = 2/3*x[1] + 4/9*x[2] - 2*x[1]*x[3]
            sage: num = p.numerator(); num
            -18*x_3*x_1 + 4*x_2 + 6*x_1

        TESTS::

            sage: num.parent()
            Infinite polynomial ring in x over Rational Field

        Check that :issue:`37756` is fixed::

            sage: R.<a> = InfinitePolynomialRing(QQ)
            sage: P.<x,y> = QQ[]
            sage: FF = P.fraction_field()
            sage: FF(a[0])  # known bug
            Traceback (most recent call last):
            ...
            TypeError: Could not find a mapping of the passed element to this ring.
        """
        P = self.parent()
        return InfinitePolynomial(P, self._p.numerator())

    def denominator(self):
        r"""
        Return a denominator of ``self``.

        First, the lcm of the denominators of the entries of ``self``
        is computed and returned. If this computation fails, the unit
        of the parent of ``self`` is returned.

        Note that some subclasses may implement their own denominator
        function.

        .. WARNING::

           This is not the denominator of the rational function
           defined by ``self``, which would always be 1 since ``self`` is a
           polynomial.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: p = 2/3*x[1] + 4/9*x[2] - 2*x[1]*x[3]
            sage: d = p.denominator(); d
            9

        TESTS::

            sage: d.parent()
            Infinite polynomial ring in x over Rational Field
        """
        P = self.parent()
        return InfinitePolynomial(P, self._p.denominator())

    @cached_method
    def variables(self):
        """
        Return the variables occurring in ``self`` (tuple of elements of some polynomial ring).

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: p = x[1] + x[2] - 2*x[1]*x[3]
            sage: p.variables()
            (x_3, x_2, x_1)
            sage: x[1].variables()
            (x_1,)
            sage: X(1).variables()
            ()
        """
        if hasattr(self._p, 'variables'):
            P = self.parent()
            return tuple(InfinitePolynomial(P, v) for v in self._p.variables())
        return ()

    def monomials(self):
        """
        Return the list of monomials in ``self``.

        The returned list is decreasingly ordered by the term ordering of
        ``self.parent()``.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: p = x[1]^3 + x[2] - 2*x[1]*x[3]
            sage: p.monomials()
            [x_3*x_1, x_2, x_1^3]

            sage: X.<x> = InfinitePolynomialRing(QQ, order='deglex')
            sage: p = x[1]^3 + x[2] - 2*x[1]*x[3]
            sage: p.monomials()
            [x_1^3, x_3*x_1, x_2]
        """
        P = self.parent()
        return [InfinitePolynomial(P, m) for m in self._p.monomials()]

    def monomial_coefficients(self, copy=None):
        """
        Return underlying dictionary with keys the exponents and values
        the coefficients of this polynomial.

        EXAMPLES::

            sage: Z.<z> = InfinitePolynomialRing(QQ)
            sage: f = 2*z[2]*z[3] + 3*z[1]^2*z[5]
            sage: f.monomial_coefficients()
            {2*B[1] + B[5]: 3, B[2] + B[3]: 2}

            sage: X.<x, y> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: g = 2*x[0]*x[2] + 3*x[1]^2*y[1]
            sage: g.monomial_coefficients()
            {B[(0, 0)] + B[(0, 2)]: 2, 2*B[(0, 1)] + B[(1, 1)]: 3}

        TESTS::

            sage: g*x[10]
            2*x_10*x_2*x_0 + 3*x_10*x_1^2*y_1
            sage: g.monomial_coefficients()
            {B[(0, 0)] + B[(0, 2)]: 2, 2*B[(0, 1)] + B[(1, 1)]: 3}

            sage: X.<x, y> = InfinitePolynomialRing(QQ)
            sage: g = 2*x[0]*x[2] + 3*x[1]^2*y[1]
            sage: g.monomial_coefficients()
            {B[(0, 0)] + B[(0, 2)]: 2, 2*B[(0, 1)] + B[(1, 1)]: 3}

            sage: Z.<z> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: f = 2*z[2]*z[3] + 3*z[1]^2*z[5]
            sage: f.monomial_coefficients()
            {2*B[1] + B[5]: 3, B[2] + B[3]: 2}
        """
        P = self.parent()
        B = P._indices
        p = self._p
        R_gens = p.parent().gens()
        if P.ngens() > 1:
            def gen_to_index(g):
                j, i = P.varname_key(str(g))
                return -j, i
        else:
            def gen_to_index(g):
                return P.varname_key(str(g))[1]

        return {B._from_dict({gen_to_index(g): e
                              for e, g in zip(m, R_gens)},
                             coerce=False, remove_zeros=True): c
                for m, c in p.monomial_coefficients().items()}

    def exponents(self):
        """
        Return the exponents of the monomials appearing in this
        polynomial.

        EXAMPLES::

            sage: Z.<z> = InfinitePolynomialRing(QQ)
            sage: f = 2*z[2]*z[3] + 3*z[1]^2*z[5]
            sage: f.exponents()
            [2*B[1] + B[5], B[2] + B[3]]

            sage: X.<x, y> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: g = 2*x[0]*x[2] + 3*x[1]^2*y[1]
            sage: g.exponents()  # random
            [B[(0, 0)] + B[(0, 2)], 2*B[(0, 1)] + B[(1, 1)]]

        TESTS::

            sage: X.<x, y> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: g*x[10]
            2*x_10*x_2*x_0 + 3*x_10*x_1^2*y_1

            sage: g.exponents()  # random
            [B[(0, 0)] + B[(0, 2)], 2*B[(0, 1)] + B[(1, 1)]]

            sage: X.<x, y> = InfinitePolynomialRing(QQ)
            sage: g = 2*x[0]*x[2] + 3*x[1]^2*y[1]
            sage: g.exponents()
            [B[(0, 0)] + B[(0, 2)], 2*B[(0, 1)] + B[(1, 1)]]

            sage: Z.<z> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: f = 2*z[2]*z[3] + 3*z[1]^2*z[5]
            sage: f.exponents()
            [2*B[1] + B[5], B[2] + B[3]]
        """
        return list(self.monomial_coefficients())

    def degrees(self):
        """
        Return an index element corresponding to the maximal
        degree of each variable in this polynomial.

        EXAMPLES::

            sage: Z.<z> = InfinitePolynomialRing(QQ)
            sage: f = 2*z[2]*z[3] + 3*z[1]^2*z[5]^4
            sage: f.degrees()
            {B[1]: 2, B[2]: 1, B[3]: 1, B[5]: 4}

            sage: X.<x, y> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: g = 2*x[0]*x[2] + 3*x[1]^2*y[1]^4
            sage: g.degrees()
            {B[(0, 0)]: 1, B[(0, 1)]: 2, B[(0, 2)]: 1, B[(1, 1)]: 4}

        TESTS::

            sage: g*x[10]
            2*x_10*x_2*x_0 + 3*x_10*x_1^2*y_1^4
            sage: g.degrees()
            {B[(0, 0)]: 1, B[(0, 1)]: 2, B[(0, 2)]: 1, B[(1, 1)]: 4}

            sage: X.<x, y> = InfinitePolynomialRing(QQ)
            sage: g = 2*x[0]*x[2] + 3*x[1]^2*y[1]^4
            sage: g.degrees()
            {B[(0, 0)]: 1, B[(0, 1)]: 2, B[(0, 2)]: 1, B[(1, 1)]: 4}

            sage: Z.<z> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: f = 2*z[2]*z[3] + 3*z[1]^2*z[5]^4
            sage: f.degrees()
            {B[1]: 2, B[2]: 1, B[3]: 1, B[5]: 4}
        """
        B = self.parent()._indices
        degs = {}
        for e in self.exponents():
            for m, d in e.monomial_coefficients().items():
                m = B.monomial(m)
                degs[m] = max(d, degs.get(m, 0))
        return degs

    @cached_method
    def max_index(self):
        r"""
        Return the maximal index of a variable occurring in ``self``, or -1 if ``self`` is scalar.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = x[1]^2 + y[2]^2 + x[1]*x[2]*y[3] + x[1]*y[4]
            sage: p.max_index()
            4
            sage: x[0].max_index()
            0
            sage: X(10).max_index()
            -1
        """
        return max([Integer(str(X).split('_')[1]) for X in self.variables()]+[-1])

    def _rmul_(self, left):
        """
        TESTS::

            sage: R.<alpha,beta> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: R.from_base_ring(4)   # indirect doctest
            4
        """
        return type(self)(self.parent(), left * self._p)

    def _lmul_(self, right):
        """
        TESTS::

            sage: R.<alpha,beta> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: alpha[3]*4   # indirect doctest
            4*alpha_3
        """
        return type(self)(self.parent(), self._p * right)

    def _div_(self, x):
        r"""
        Division of Infinite Polynomials.

        EXAMPLES:

        Division by a rational over `\QQ`::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: x[0]/2
            1/2*x_0

        Division by an integer over `\ZZ`::

            sage: R.<x> = InfinitePolynomialRing(ZZ, implementation='sparse')
            sage: p = x[3] + x[2]
            sage: q = p/2
            sage: q
            1/2*x_3 + 1/2*x_2
            sage: q.parent()
            Infinite polynomial ring in x over Rational Field

        Division by a nonzero element::

            sage: R.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: 1/x[1]
            1/x_1
            sage: (x[0]/x[0])
            x_0/x_0
            sage: qt = 1/x[2] + 2/x[1]; qt
            (2*x_2 + x_1)/(x_2*x_1)
            sage: qt.parent()
            Fraction Field of Infinite polynomial ring in x over Rational Field

            sage: z = 1/(x[2]*(x[1]+x[2]))+1/(x[1]*(x[1]+x[2]))
            sage: z.parent()
            Fraction Field of Infinite polynomial ring in x over Rational Field
            sage: factor(z)                                                             # needs sage.libs.singular
            x_1^-1 * x_2^-1
        """
        if not x.variables():
            p = self.base_ring()(x._p)
            divisor = self.base_ring().one() / p  # use induction
            OUTP = self.parent().tensor_with_ring(divisor.base_ring())
            return OUTP(self) * OUTP(divisor)
        else:
            from sage.rings.fraction_field_element import FractionFieldElement
            field = self.parent().fraction_field()
            # there remains a problem in reduction
            return FractionFieldElement(field, self, x, reduce=False)

    def factor(self, proof=None):
        """
        Return the factorization of this polynomial.

        INPUT:

        - ``proof`` -- ignored

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQbar)
            sage: factor(x[3]^2 - x[1]^2)
            (x_3 - x_1) * (x_3 + x_1)

        TESTS::

            sage: P = factor(x[3]^2 - x[1]^2)[0][0].parent()
            sage: P
            Infinite polynomial ring in x over Algebraic Field
            sage: P == R
            True
        """
        P = self.parent()
        f = self._p.factor(proof=proof)
        return Factorization([(InfinitePolynomial(P, p), e) for p, e in f],
                             unit=f.unit(),
                             cr=f._cr(),
                             sort=False,
                             simplify=False)

    @cached_method
    def lm(self):
        """
        The leading monomial of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30] + x[10]*y[1]^3*x[1]^2
            sage: p.lm()
            x_10*x_1^2*y_1^3
        """
        if hasattr(self._p, 'lm'):
            return InfinitePolynomial(self.parent(), self._p.lm())
        if self._p == 0:
            return self
        if hasattr(self._p, 'variable_name'):  # if it is univariate
            return InfinitePolynomial(self.parent(),
                                      self._p.parent().gen() ** max(self._p.exponents()))
        return self  # if it is scalar

    @cached_method
    def lc(self):
        """
        The coefficient of the leading term of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30] + 3*x[10]*y[1]^3*x[1]^2
            sage: p.lc()
            3
        """
        if hasattr(self._p, 'lc'):
            return self._p.lc()
        if hasattr(self._p, 'variable_name'):  # univariate case
            return self._p.leading_coefficient()
        # scalar case
        return self._p

    @cached_method
    def lt(self):
        """
        The leading term (= product of coefficient and monomial) of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30] + 3*x[10]*y[1]^3*x[1]^2
            sage: p.lt()
            3*x_10*x_1^2*y_1^3
        """
        if hasattr(self._p, 'lt'):
            return InfinitePolynomial(self.parent(), self._p.lt())
        if self._p == 0:
            return self
        if hasattr(self._p, 'variable_name'):  # if it is univariate
            return InfinitePolynomial(self.parent(), self._p.leading_coefficient()*self._p.parent().gen()**max(self._p.exponents()))
        return self  # if it is scalar

    def tail(self):
        """
        The tail of ``self`` (this is ``self`` minus its leading term).

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30] + 3*x[10]*y[1]^3*x[1]^2
            sage: p.tail()
            2*x_10*y_30
        """
        return self-self.lt()

    def squeezed(self):
        """
        Reduce the variable indices occurring in ``self``.

        OUTPUT:

        Apply a permutation to ``self`` that does not change the order of
        the variable indices of ``self`` but squeezes them into the range
        1,2,...

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: p = x[1]*y[100] + x[50]*y[1000]
            sage: p.squeezed()
            x_2*y_4 + x_1*y_3
        """
        Indices = set([0] + [Integer(str(Y).split('_')[1])
                             for Y in self.variables()])
        Indices = sorted(Indices)

        def P(n):
            return Indices.index(n) if n in Indices else n

        return self**P

    def footprint(self):
        """
        Leading exponents sorted by index and generator.

        OUTPUT: ``D``; dictionary whose keys are the occurring variable indices

        ``D[s]`` is a list ``[i_1,...,i_n]``, where ``i_j`` gives the
        exponent of ``self.parent().gen(j)[s]`` in the leading
        term of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = x[30]*y[1]^3*x[1]^2 + 2*x[10]*y[30]
            sage: sorted(p.footprint().items())
            [(1, [2, 3]), (30, [1, 0])]

        TESTS:

        This is a test whether it also works when the underlying polynomial ring is
        not implemented in libsingular::

            sage: X.<x> = InfinitePolynomialRing(ZZ)
            sage: Y.<y,z> = X[]
            sage: Z.<a> = InfinitePolynomialRing(Y)
            sage: Z
            Infinite polynomial ring in a over Multivariate Polynomial Ring in y, z over Infinite polynomial ring in x over Integer Ring
            sage: type(Z._P)
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_domain_with_category'>
            sage: p = a[12]^3*a[2]^7*a[4] + a[4]*a[2]
            sage: sorted(p.footprint().items())
            [(2, [7]), (4, [1]), (12, [3])]
        """
        if not self._has_footprint:
            P = self.parent()
            l = len(P._names)
            # get the pairs (shift,exponent) of the leading monomial, indexed by the variable names
            Vars = self._p.parent().variable_names()
            from sage.rings.polynomial.multi_polynomial import MPolynomial_libsingular
            if isinstance(self._p, MPolynomial_libsingular):
                L = [(Vars[i].split('_'), e) for i, e in enumerate(self._p.lm().exponents(as_ETuples=False)[0]) if e]
            elif hasattr(self._p, 'lm'):
                # self._p  is multivariate, but not libsingular, hence,
                # exponents is slow and does not accept the optional argument as_ETuples.
                # Thus, fall back to regular expressions
                L = P._find_varpowers.findall(repr(self.lm()._p))
                L = [((x[0:2]), int(x[2]) if x[2] else 1) for x in L]
            else:  # it is a univariate polynomial -- this should never happen, but just in case...
                L = [(Vars[0].split('_'), self._p.degree())]
            for t in L:
                n = t[0][0]       # the variable *n*ame
                s = int(t[0][1])  # the variable *s*hift
                if s not in self._footprint:
                    self._footprint[s] = [0]*l
                self._footprint[s][P._name_dict[n]] = t[1]   # the exponent
            self._has_footprint = True
        return self._footprint

    def symmetric_cancellation_order(self, other):
        """
        Comparison of leading terms by Symmetric Cancellation Order, `<_{sc}`.

        INPUT:

        - ``self``, ``other`` -- two Infinite Polynomials

        ASSUMPTION:

        Both Infinite Polynomials are nonzero.

        OUTPUT:

        ``(c, sigma, w)``, where

        * c = -1,0,1, or None if the leading monomial of ``self`` is smaller, equal,
          greater, or incomparable with respect to ``other`` in the monomial
          ordering of the Infinite Polynomial Ring
        * sigma is a permutation witnessing
          ``self`` `<_{sc}` ``other`` (resp. ``self`` `>_{sc}` ``other``)
          or is 1 if ``self.lm()==other.lm()``
        * w is 1 or is a term so that
          ``w*self.lt()^sigma == other.lt()`` if `c\\le 0`, and
          ``w*other.lt()^sigma == self.lt()`` if `c=1`

        THEORY:

        If the Symmetric Cancellation Order is a well-quasi-ordering
        then computation of Groebner bases always terminates. This is
        the case, e.g., if the monomial order is lexicographic. For
        that reason, lexicographic order is our default order.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: (x[2]*x[1]).symmetric_cancellation_order(x[2]^2)
            (None, 1, 1)
            sage: (x[2]*x[1]).symmetric_cancellation_order(x[2]*x[3]*y[1])
            (-1, [2, 3, 1], y_1)
            sage: (x[2]*x[1]*y[1]).symmetric_cancellation_order(x[2]*x[3]*y[1])
            (None, 1, 1)
            sage: (x[2]*x[1]*y[1]).symmetric_cancellation_order(x[2]*x[3]*y[2])
            (-1, [2, 3, 1], 1)
        """
        PARENT = self.parent()
        other = PARENT(other)
        slt = self.lt()
        olt = other.lt()
        if self.lm() == other.lm():
            if olt == 0:
                return (0, 1, 1)
            return (0, 1, self.lc() / other.lc())
        if self.lm() < other.lm():
            rawcmp = -1
            Fsmall = {k: list(v) for k, v in self.footprint().items()}
            Fbig = {k: list(v) for k, v in other.footprint().items()}
            ltsmall = slt
            ltbig = olt
        else:
            rawcmp = 1
            Fbig = {k: list(v) for k, v in self.footprint().items()}
            Fsmall = {k: list(v) for k, v in other.footprint().items()}
            ltbig = slt
            ltsmall = olt
        # Case 1: one of the Infinite Polynomials is scalar.
        if not Fsmall:
            return (rawcmp, 1, ltbig/ltsmall)
        # "not Fbig" is now impossible, because we only consider *global* monomial orderings.
        # These are the occurring shifts:
        Lsmall = sorted(Fsmall.keys())
        Lbig = sorted(Fbig.keys())
        P = list(range(Lbig[-1] + 1))
        gens = list(range(PARENT.ngens()))
        if Lsmall[0] == 0:
            if 0 not in Fbig:
                return (None, 1, 1)
            Lsmall.pop(0)
            Lbig.pop(0)
            ExpoSmall = Fsmall[0]
            ExpoBig = Fbig[0]
            for k in gens:
                if ExpoBig[k] < ExpoSmall[k]:
                    return (None, 1, 1)
                ExpoBig[k] -= ExpoSmall[k]
        lenBig = len(Lbig)
        j = -1  # will have Lbig[j] -- a shift of the bigger polynomial
        for i in Lsmall:  # i is a shift of the smaller polynomial
            j += 1
            ExpoSmall = Fsmall[i]
            while j < lenBig:
                found = False
                if Lbig[j] >= i:
                    ExpoBigSave = list(Fbig[Lbig[j]])
                    ExpoBig = Fbig[Lbig[j]]
                    found = True
                    for k in gens:
                        if ExpoBig[k] < ExpoSmall[k]:
                            found = False
                            Fbig[Lbig[j]] = ExpoBigSave
                            break
                        ExpoBig[k] -= ExpoSmall[k]
                if found:
                    break
                j += 1
            if j == lenBig:
                # no "increasing" permutation transforms
                # the smaller monomial into a factor of
                # the bigger monomial
                return (None, 1, 1)
            tmp = P[i]
            P[i] = Lbig[j]
            P[Lbig[j]] = tmp
        # now, P defines an 'up-shift' permutation, slt^P divides olt, and
        # Fbig contains the exponents for olt/slt^P.
        OUT = PARENT(PARENT._base(ltbig.lc() / ltsmall.lc()))
        for shift, Expo in Fbig.items():
            for g in gens:
                if Expo[g]:
                    OUT *= PARENT.gen(g)[shift] ** Expo[g]
        from sage.combinat.permutation import Permutation
        return (rawcmp, Permutation(P[1:]), OUT)

    # Essentials for Buchberger
    def reduce(self, I, tailreduce=False, report=None):
        r"""
        Symmetrical reduction of ``self`` with respect to a symmetric ideal (or list of Infinite Polynomials).

        INPUT:

        - ``I`` -- a :class:`~sage.rings.polynomial.symmetric_ideal.SymmetricIdeal` or a list
          of Infinite Polynomials
        - ``tailreduce`` -- boolean (default: ``False``); *tail reduction* is performed if this
          parameter is ``True``.
        - ``report`` -- object (default: ``None``); if not ``None``, some information on the
          progress of computation is printed, since reduction of huge polynomials may take
          a long time

        OUTPUT: symmetrical reduction of ``self`` with respect to ``I``, possibly with tail reduction

        THEORY:

        Reducing an element `p` of an Infinite Polynomial Ring `X` by
        some other element `q` means the following:

        1. Let `M` and `N` be the leading terms of `p` and `q`.
        2. Test whether there is a permutation `P` that does not does not diminish the variable
           indices occurring in `N` and preserves their order, so that there is some term `T\in X`
           with `TN^P = M`. If there is no such permutation, return `p`
        3. Replace `p` by `p-T q^P` and continue with step 1.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = y[1]^2*y[3] + y[2]*x[3]^3
            sage: p.reduce([y[2]*x[1]^2])
            x_3^3*y_2 + y_3*y_1^2

        The preceding is correct: If a permutation turns
        ``y[2]*x[1]^2`` into a factor of the leading monomial
        ``y[2]*x[3]^3`` of ``p``, then it interchanges the variable
        indices 1 and 2; this is not allowed in a symmetric
        reduction. However, reduction by ``y[1]*x[2]^2`` works, since
        one can change variable index 1 into 2 and 2 into 3::

            sage: p.reduce([y[1]*x[2]^2])                                               # needs sage.libs.singular
            y_3*y_1^2

        The next example shows that tail reduction is not done, unless
        it is explicitly advised.  The input can also be a Symmetric
        Ideal::

            sage: I = (y[3])*X
            sage: p.reduce(I)
            x_3^3*y_2 + y_3*y_1^2
            sage: p.reduce(I, tailreduce=True)                                          # needs sage.libs.singular
            x_3^3*y_2

        Last, we demonstrate the ``report`` option::

            sage: p = x[1]^2 + y[2]^2 + x[1]*x[2]*y[3] + x[1]*y[4]
            sage: p.reduce(I, tailreduce=True, report=True)                             # needs sage.libs.singular
            :T[2]:>
            >
            x_1^2 + y_2^2

        The output ':' means that there was one reduction of the
        leading monomial. 'T[2]' means that a tail reduction was
        performed on a polynomial with two terms. At '>', one round of
        the reduction process is finished (there could only be several
        non-trivial rounds if `I` was generated by more than one
        polynomial).
        """
        from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
        if hasattr(I, 'gens'):
            I = I.gens()
        if (not I):
            return self
        I = list(I)
        S = SymmetricReductionStrategy(self.parent(), I, tailreduce)
        return S.reduce(self, report=report)

    # Further methods
    def stretch(self, k):
        r"""
        Stretch ``self`` by a given factor.

        INPUT:

        - ``k`` -- integer

        OUTPUT: replace `v_n` with `v_{n\cdot k}` for all generators `v_\ast`
        occurring in ``self``

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1] + x[2]
            sage: a.stretch(2)
            x_4 + x_2 + x_0

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1] + y[0]*y[1]; a
            x_1 + x_0 + y_1*y_0
            sage: a.stretch(2)
            x_2 + x_0 + y_2*y_0

        TESTS:

        The following would hardly work in a dense implementation,
        because an underlying polynomial ring with 6001 variables
        would be created. This is avoided in the sparse
        implementation::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: a = x[2] + x[3]
            sage: a.stretch(2000)
            x_6000 + x_4000
        """
        def P(n):
            return k*n
        return self ** P

    def __iter__(self):
        """
        Return an iterator over all pairs ``(coefficient, monomial)``
        of this polynomial.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + 2*x[1] + y[0]*y[1]
            sage: list(a)
            [(2, x_1), (1, x_0), (1, y_1*y_0)]
        """
        return iter((coefficient,
                     self.__class__(self.parent(), monomial))
                    for coefficient, monomial in self._p)


class InfinitePolynomial_sparse(InfinitePolynomial):
    """
    Element of a sparse Polynomial Ring with a Countably Infinite Number of Variables.

    INPUT:

    - ``A`` -- an Infinite Polynomial Ring in sparse implementation
    - ``p`` -- a *classical* polynomial that can be interpreted in ``A``

    Of course, one should not directly invoke this class, but rather
    construct elements of ``A`` in the usual way.

    EXAMPLES::

        sage: A.<a> = QQ[]
        sage: B.<b,c> = InfinitePolynomialRing(A, implementation='sparse')
        sage: p = a*b[100] + 1/2*c[4]
        sage: p
        a*b_100 + 1/2*c_4
        sage: p.parent()
        Infinite polynomial ring in b, c
         over Univariate Polynomial Ring in a over Rational Field
        sage: p.polynomial().parent()
        Multivariate Polynomial Ring in b_100, b_0, c_4, c_0
         over Univariate Polynomial Ring in a over Rational Field
    """
    def __call__(self, *args, **kwargs):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: a = x[0] + x[1]
            sage: a(x_0=2,x_1=x[1])
            x_1 + 2
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field
            sage: a(x_1=3)
            x_0 + 3
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field
            sage: a(x_1=x[100])
            x_100 + x_0

            sage: M = matrix([[1,1], [2,0]])                                            # needs sage.modules
            sage: a(x_1=M)                                                              # needs sage.modules
            [x_0 + 1       1]
            [      2     x_0]
        """
        # Replace any InfinitePolynomials by their underlying polynomials
        if hasattr(self._p, 'variables'):
            V = [str(x) for x in self._p.variables()]
        else:
            V = []
        for kw, value in kwargs.items():
            if isinstance(value, InfinitePolynomial):
                kwargs[kw] = value._p
                V.append(kw)
                if hasattr(value._p, 'variables'):
                    V.extend([str(x) for x in value._p.variables()])
        args = list(args)
        for i, arg in enumerate(args):
            if isinstance(arg, InfinitePolynomial):
                args[i] = arg._p
                if hasattr(arg._p, 'variables'):
                    V.extend([str(x) for x in arg._p.variables()])
        V = list(set(V))
        V.sort(key=self.parent().varname_key, reverse=True)
        if V:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(), V, order=self.parent()._order)
        else:
            return self
        res = R(self._p)(*args, **kwargs)
        try:
            from sage.misc.sage_eval import sage_eval
            return sage_eval(repr(res), self.parent().gens_dict())
        except Exception:
            return res

    def _common_polynomial_ring(self, x):
        r"""
        Return the polynomial ring that contains all variables of
        ``self`` and of ``x``.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: x[10]._common_polynomial_ring(x[5])
            Multivariate Polynomial Ring in x_10, x_5, x_0 over Rational Field
        """
        VarList = set(self._p.parent().variable_names())
        VarList.update(x._p.parent().variable_names())
        VarList = sorted(VarList, key=self.parent().varname_key, reverse=True)
        if VarList:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(), VarList, order=self.parent()._order)
        else:
            # TODO: this is currently dead code
            R = self._p.base_ring()
        return R

    # Basic arithmetics
    def _add_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: x[1] + x[2]  # indirect doctest
            x_2 + x_1

        Check adding from a different parent::

            sage: Y.<x_0> = PolynomialRing(QQ)
            sage: x[0] - x_0
            0

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar, implementation='sparse')
            sage: x[1] + x[3]
            x_3 + x_1
        """
        try:
            result = self._p + x._p
        except TypeError:
            R = self._common_polynomial_ring(x)
            result = R(self._p) + R(x._p)

        return InfinitePolynomial_sparse(self.parent(), result)

    def _mul_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(ZZ)
            sage: x[2] * x[1]  # indirect doctest
            x_2*x_1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar, implementation='sparse')
            sage: x[1] * x[3]
            x_3*x_1
        """
        try:
            result = self._p * x._p
        except TypeError:
            R = self._common_polynomial_ring(x)
            result = R(self._p) * R(x._p)

        return InfinitePolynomial_sparse(self.parent(), result)

    def _sub_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: x[2] - x[1] # indirect doctest
            x_2 - x_1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar, implementation='sparse')
            sage: x[1] - x[3]
            -x_3 + x_1
        """
        try:
            result = self._p - x._p
        except TypeError:
            R = self._common_polynomial_ring(x)
            result = R(self._p) - R(x._p)

        return InfinitePolynomial_sparse(self.parent(), result)

    def _floordiv_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(ZZ, implementation="sparse")
            sage: x[2] // x[2]  # indirect doctest
            1
            sage: (x[2]^2 - 1) // (x[2] + 1)
            x_2 - 1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar, implementation='sparse')
            sage: (x[2]^2 - 1) // (x[2] + 1)
            x_2 - 1
        """
        try:
            result = self._p // x._p
        except TypeError:
            R = self._common_polynomial_ring(x)
            result = R(self._p) // R(x._p)

        return InfinitePolynomial_sparse(self.parent(), result)

    def __pow__(self, n):
        """
        Exponentiation by an integer, or action by a callable object.

        NOTE:

        The callable object must accept nonnegative integers as input
        and return nonnegative integers. Typical use case is a
        permutation, that will result in the corresponding permutation
        of variables.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: p = x[10]*y[2] + 2*x[1]*y[3]
            sage: P = Permutation(((1,2),(3,4,5)))
            sage: p^P # indirect doctest
            x_10*y_1 + 2*x_2*y_4

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar, implementation='sparse')
            sage: (x[3] + x[1])^2
            x_3^2 + 2*x_3*x_1 + x_1^2
        """
        P = self.parent()
        if callable(n):
            if (self._p.parent() == self._p.base_ring()):
                return self
            if not (hasattr(self._p, 'variables') and self._p.variables()):
                return self
            if hasattr(n, 'to_cycles') and hasattr(n, '__len__'):  # duck typing Permutation
                # auxiliary function, necessary since n(m) raises an error if m>len(n)
                l = len(n)

                def p(m):
                    return n(m) if 0 < m <= l else m
            else:  # Permutation group element
                p = n

            def q(s):
                return s[0] + '_' + str(p(ZZ(s[1])))

            newVars = [q(X.split('_')) for X in self._p.parent().variable_names()]
            if not newVars:
                return self
            copyVars = copy.copy(newVars)
            newVars = list(set(list(self._p.parent().variable_names())+newVars))
            newVars.sort(key=self.parent().varname_key, reverse=True)
            if newVars == list(self._p.parent().variable_names()):
                newR = self._p.parent()
            else:
                from sage.rings.polynomial.polynomial_ring_constructor import (
                    PolynomialRing,
                )
                newR = PolynomialRing(self._p.base_ring(), newVars, order=P._order)
            mapR = self._p.parent().hom(copyVars, newR)
            return InfinitePolynomial_sparse(self.parent(), mapR(self._p))
        return InfinitePolynomial_sparse(self.parent(), self._p**n)

    # Basic tools for Buchberger algorithm:
    # order, leading term/monomial, symmetric cancellation order

    def _richcmp_(self, x, op):
        r"""
        Comparison of Infinite Polynomials.

        NOTE:

        Let x and y be generators of the parent of ``self``. We only consider
        monomial orderings in which x[m] > y[n] iff x appears earlier in the
        list of generators than y, or x==y and m>n

        Under this restriction, the monomial ordering can be 'lex' (default),
        'degrevlex' or 'deglex'.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: a = x[10]^3
            sage: b = x[1] + x[2]
            sage: c = x[1] + x[2]
            sage: d = y[1] + x[2]
            sage: a == a # indirect doctest
            True
            sage: b == c # indirect doctest
            True
            sage: a == b # indirect doctest
            False
            sage: c > d # indirect doctest
            True

        TESTS:

        A classical and an infinite sparse polynomial ring. Note that
        the Sage coercion system allows comparison only if a common
        parent for the two rings can be constructed. This is why we
        have to have the order 'degrevlex'::

            sage: X.<x,y> = InfinitePolynomialRing(ZZ, order='degrevlex', implementation='sparse')
            sage: Y.<z,x_3,x_1> = QQ[]
            sage: x[3] == x_3 # indirect doctest
            True

        Two infinite polynomial rings in different implementation and
        order::

            sage: Y = InfinitePolynomialRing(QQ, ['x','y'], order='deglex', implementation='dense')
            sage: x[2] == Y(x[2]) # indirect doctest
            True

        An example in which a previous version had failed::

            sage: X.<x,y> = InfinitePolynomialRing(GF(3), order='degrevlex', implementation='sparse')
            sage: p = Y('x_3*x_0^2 + x_0*y_3*y_0')
            sage: q = Y('x_1*x_0^2 + x_0*y_1*y_0')
            sage: p < q
            False
        """
        # We can assume that self.parent() is x.parent(),
        # but of course the underlying polynomial rings
        # may be widely different, and the sage coercion
        # system can't guess what order we want.
        from sage.structure.element import parent
        R1 = parent(self._p)
        R2 = parent(x._p)
        if ((hasattr(R1, 'has_coerce_map_from') and R1.has_coerce_map_from(R2))
            or (hasattr(R2, 'has_coerce_map_from') and R2.has_coerce_map_from(R1))):
            return richcmp(self._p, x._p, op)
        R = self._common_polynomial_ring(x)
        if (self._p.parent() is self._p.base_ring()) or not self._p.parent().gens():
            fself = self._p.base_ring()
        else:
            fself = self._p.parent().hom(self._p.parent().variable_names(), R)
        if (x._p.parent() is x._p.base_ring()) or not x._p.parent().gens():
            fx = x._p.base_ring()
        else:
            fx = x._p.parent().hom(x._p.parent().variable_names(), R)
        return richcmp(fself(self._p), fx(x._p), op)

    def gcd(self, x):
        """
        Compute the greatest common divisor.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: p1 = x[0] + x[1]^2
            sage: gcd(p1, p1 + 3)
            1
            sage: gcd(p1, p1) == p1
            True
        """
        try:
            result = self._p.gcd(x._p)
        except (ValueError, TypeError):
            R = self._common_polynomial_ring(x)
            result = R(self._p).gcd(R(x._p))

        return InfinitePolynomial_sparse(self.parent(), result)

    def quo_rem(self, x):
        """
        Return quotient and remainder.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: p = 1 + 3*x[0]*x[1] + 2*x[2]
            sage: q = x[0] - 1
            sage: p.quo_rem(q)
            (3*x_1, 2*x_2 + 3*x_1 + 1)

        TESTS::

            sage: R.<x> = InfinitePolynomialRing(QQbar, implementation="sparse")
            sage: p = 1 + 3*x[0]*x[1] + 2*x[2]
            sage: q = x[0] - 1
            sage: p.quo_rem(q)
            (3*x_1, 2*x_2 + 3*x_1 + 1)
        """
        try:
            result = self._p.quo_rem(x._p)
        except (ValueError, TypeError):
            R = self._common_polynomial_ring(x)
            result = R(self._p).quo_rem(R(x._p))

        return (InfinitePolynomial_sparse(self.parent(), result[0]),
                InfinitePolynomial_sparse(self.parent(), result[1]))

    def monomial_coefficient(self, mon):
        """
        Return the base ring element that is the coefficient of ``mon``
        in ``self``.

        This function contrasts with the function :meth:`coefficient`,
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` -- a monomial in the parent of ``self``

        OUTPUT: coefficient in base ring

        .. SEEALSO::

            For coefficients in a base ring of fewer variables,
            look at :meth:`coefficient`.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: f = 2*x[0]*x[2] + 3*x[1]^2
            sage: c = f.monomial_coefficient(x[1]^2); c
            3
            sage: c.parent()
            Rational Field

            sage: c = f.coefficient(x[2]); c
            2*x_0
            sage: c.parent()
            Infinite polynomial ring in x over Rational Field

        TESTS::

            sage: R.<a> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: m = a[1]*a[2]
            sage: p = a[1]*a[2] + a[1]*a[2]*a[4]
            sage: p.monomial_coefficient(m)
            1

            sage: p = R(25)
            sage: p.monomial_coefficient(R(1))
            25
        """
        P = self.parent()
        if parent(mon) is not P:
            raise TypeError("mon must be a monomial in the parent of self")

        VarList = set(self._p.parent().variable_names())
        if any(str(v) not in VarList for v in mon._p.variables()):
            return P.zero()

        if VarList:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

            VarList.update(mon._p.parent().variable_names())
            VarList = sorted(VarList, key=P.varname_key, reverse=True)
            R = PolynomialRing(self._p.base_ring(), VarList, order=P._order)
            return R(self._p).monomial_coefficient(R(mon._p))

        return self._p.constant_coefficient()

    def coefficient(self, x):
        """
        Return the coefficient of a monomial in ``self``.

        INPUT:

        - a monomial (element of the parent of ``self``) or
        - a dictionary that describes a monomial (the keys
          are variables of the parent of ``self``, the values
          are the corresponding exponents)

        EXAMPLES:

        We can get the coefficient in front of monomials::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: a = 2*x[0]*x[1] + x[1] + x[2]
            sage: a.coefficient(x[0])
            2*x_1
            sage: a.coefficient(x[1])
            2*x_0 + 1
            sage: a.coefficient(x[2])
            1
            sage: a.coefficient(x[0]*x[1])
            2

        We can also pass in a dictionary::

            sage: a.coefficient({x[0]:1, x[1]:1})
            2

        TESTS:

        Check that :issue:`40504` is fixed::

            sage: R.<a> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: p = a[1]*a[2] + a[3]*a[4]
            sage: p.coefficient(a[5])
            0

            sage: R.<a> = InfinitePolynomialRing(QQ, implementation="sparse")
            sage: P.<x> = LazyPowerSeriesRing(R)
            sage: f1 = P(lambda n: -a[n] if n else 1)
            sage: diff(log(f1))[3]
            -4*a_4 - 4*a_3*a_1 - 2*a_2^2 - 4*a_2*a_1^2 - a_1^4
            sage: diff(log(f1))[3].coefficient(a[3]*a[1])
            -4
        """
        if isinstance(x, dict):
            if x:
                I = iter(x)
                K = next(I)
                del x[K]
                return self.coefficient(K).coefficient(x)

            return self

        try:
            result = self._p.coefficient(x._p)
        except TypeError:
            R = self._common_polynomial_ring(x)
            result = R(self._p).coefficient(R(x._p))

        return InfinitePolynomial_sparse(self.parent(), result)


class InfinitePolynomial_dense(InfinitePolynomial):
    """
    Element of a dense Polynomial Ring with a Countably Infinite Number of Variables.

    INPUT:

    - ``A`` -- an Infinite Polynomial Ring in dense implementation
    - ``p`` -- a *classical* polynomial that can be interpreted in ``A``

    Of course, one should not directly invoke this class, but rather
    construct elements of ``A`` in the usual way.
    """

    def __call__(self, *args, **kwargs):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1]
            sage: a(x_0=2, x_1=x[1])
            x_1 + 2
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field
            sage: a(x_1=3)
            x_0 + 3
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field

            sage: a(x_1=x[100])
            x_100 + x_0
        """
        # Replace any InfinitePolynomials by their underlying polynomials
        for kw, value in kwargs.items():
            if isinstance(value, InfinitePolynomial):
                kwargs[kw] = value._p
        args = list(args)
        for i, arg in enumerate(args):
            if isinstance(arg, InfinitePolynomial):
                args[i] = arg._p
        self._p = self.parent().polynomial_ring()(self._p)
        res = self._p(*args, **kwargs)
        try:
            return self.parent()(res)
        except ValueError:
            return res

    def _richcmp_(self, x, op):
        r"""
        TESTS:

        A classical and an infinite polynomial ring::

            sage: X.<x,y> = InfinitePolynomialRing(ZZ, order='degrevlex')
            sage: Y.<z,x_3,x_1> = QQ[]
            sage: x[3] == x_3
            True

        Two infinite polynomial rings with different order and
        implementation::

            sage: Y = InfinitePolynomialRing(QQ,['x','y'], order='deglex', implementation='sparse')
            sage: x[2] == Y(x[2])
            True

        An example in which a previous version had failed::

            sage: X.<x,y> = InfinitePolynomialRing(GF(3), order='degrevlex', implementation='dense')
            sage: p = Y('x_3*x_0^2 + x_0*y_3*y_0')
            sage: q = Y('x_1*x_0^2 + x_0*y_1*y_0')
            sage: p < q
            False
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return richcmp(self._p, x._p, op)

    # Basic arithmetics
    def _add_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[1] + x[2]  # indirect doctest
            x_2 + x_1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar)
            sage: x[1] + x[3]
            x_3 + x_1
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(P, self._p + x._p)

    def _mul_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[2]*x[1]  # indirect doctest
            x_2*x_1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar)
            sage: x[1]*x[3]
            x_3*x_1
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(P, self._p * x._p)

    def _sub_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[2] - x[1] # indirect doctest
            x_2 - x_1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar)
            sage: x[1] - x[3]
            -x_3 + x_1
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(P, self._p - x._p)

    def _floordiv_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(ZZ)
            sage: x[2] // x[2]  # indirect doctest
            1
            sage: (x[2]^2 - 1) // (x[2] + 1)
            x_2 - 1

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQbar)
            sage: (x[2]^2 - 1) // (x[2] + 1)
            x_2 - 1
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(P, self._p // x._p)

    def __pow__(self, n):
        """
        Exponentiation by an integer, or action by a callable object.

        NOTE:

        The callable object must accept nonnegative integers as input
        and return nonnegative integers. Typical use case is a
        permutation, that will result in the corresponding permutation
        of variables.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: x[10]^3
            x_10^3
            sage: p = x[10]*y[2] + 2*x[1]*y[3]
            sage: P = Permutation(((1,2),(3,4,5)))
            sage: p^P
            x_10*y_1 + 2*x_2*y_4
        """
        P = self.parent()
        if callable(n):
            if (self._p.parent() == self._p.base_ring()):
                return self
            if not (hasattr(self._p, 'variables') and self._p.variables()):
                return self
            if hasattr(n, 'to_cycles') and hasattr(n, '__len__'):  # duck typing Permutation
                # auxiliary function, necessary since n(m) raises an error if m>len(n)
                l = len(n)

                def p(m):
                    return n(m) if 0 < m <= l else m
            else:  # Permutation group element
                p = n

            # determine whether the maximal index must be raised
            oldMax = P._max
            newMax = max([p(X) for X in range(oldMax+1)]+[oldMax])
            if newMax > P._max:
                P.gen()[newMax]
            self._p = P._P(self._p)
            # next, determine the images of variable names
            PP = P._P
            PPgens = PP.gens()

            newVars = []
            sh = PP.ngens() // P.ngens() - 1
            blocklength = sh
            nM = sh + 1
            for i in range(P.ngens()):
                newVars.extend([PPgens[sh-p(j)] for j in range(blocklength, -1, -1)])
                sh += nM
            mapR = PP.hom(newVars, PP)
            return InfinitePolynomial_dense(P, mapR(self._p))

        # else, n is supposed to be an integer
        return InfinitePolynomial_dense(P, self._p**n)

    def gcd(self, x):
        """
        Compute the greatest common divisor.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQ)
            sage: p1 = x[0] + x[1]^2
            sage: gcd(p1, p1 + 3)
            1
            sage: gcd(p1, p1) == p1
            True
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(P, self._p.gcd(x._p))

    def quo_rem(self, x):
        """
        Return quotient and remainder.

        EXAMPLES::

            sage: R.<x> = InfinitePolynomialRing(QQ)
            sage: p = 1 + 3*x[0]*x[1] + 2*x[2]
            sage: q = x[0] - 1
            sage: p.quo_rem(q)
            (3*x_1, 2*x_2 + 3*x_1 + 1)

        TESTS::

            sage: R.<x> = InfinitePolynomialRing(QQbar)
            sage: p = 1 + 3*x[0]*x[1] + 2*x[2]
            sage: q = x[0] - 1
            sage: p.quo_rem(q)
            (3*x_1, 2*x_2 + 3*x_1 + 1)
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        result = (self._p).quo_rem(x._p)

        return (InfinitePolynomial_dense(P, result[0]),
                InfinitePolynomial_dense(P, result[1]))

    def monomial_coefficient(self, mon):
        """
        Return the base ring element that is the coefficient of ``mon``
        in ``self``.

        This function contrasts with the function :meth:`coefficient`,
        which returns the coefficient of a monomial viewing this
        polynomial in a polynomial ring over a base ring having fewer
        variables.

        INPUT:

        - ``mon`` -- a monomial in the parent of ``self``

        OUTPUT: coefficient in base ring

        .. SEEALSO::

            For coefficients in a base ring of fewer variables,
            look at :meth:`coefficient`.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: f = 2*x[0]*x[2] + 3*x[1]^2
            sage: c = f.monomial_coefficient(x[1]^2); c
            3
            sage: c.parent()
            Rational Field

            sage: c = f.coefficient(x[2]); c
            2*x_0
            sage: c.parent()
            Infinite polynomial ring in x over Rational Field

        TESTS::

            sage: R.<a> = InfinitePolynomialRing(QQ)
            sage: m = a[1]*a[2]
            sage: p = a[1]*a[2] + a[1]*a[2]*a[4]
            sage: p.monomial_coefficient(m)
            1

            sage: p = R(25)
            sage: p.monomial_coefficient(R(1))
            25
        """
        P = self.parent()
        if parent(mon) is not P:
            raise TypeError("mon must be a monomial in the parent of self")
        self._p = P._P(self._p)
        mon._p = P._P(mon._p)
        return self._p.monomial_coefficient(mon._p)

    def coefficient(self, monomial):
        """
        Return the coefficient of a monomial in ``self``.

        INPUT:

        - a monomial (element of the parent of ``self``) or
        - a dictionary that describes a monomial (the keys
          are variables of the parent of ``self``, the values
          are the corresponding exponents)

        EXAMPLES:

        We can get the coefficient in front of monomials::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = 2*x[0]*x[1] + x[1] + x[2]
            sage: a.coefficient(x[0])
            2*x_1
            sage: a.coefficient(x[1])
            2*x_0 + 1
            sage: a.coefficient(x[2])
            1
            sage: a.coefficient(x[0]*x[1])
            2

        We can also pass in a dictionary::

            sage: a.coefficient({x[0]:1, x[1]:1})
            2

        TESTS:

        Check that :issue:`40504` is fixed::

            sage: R.<a> = InfinitePolynomialRing(QQ)
            sage: p = a[1]*a[2] + a[3]*a[4]
            sage: p.coefficient(a[5])
            0

            sage: R.<a> = InfinitePolynomialRing(QQ)
            sage: P.<x> = LazyPowerSeriesRing(R)
            sage: f1 = P(lambda n: -a[n] if n else 1)
            sage: diff(log(f1))[3]
            -4*a_4 - 4*a_3*a_1 - 2*a_2^2 - 4*a_2*a_1^2 - a_1^4
            sage: diff(log(f1))[3].coefficient(a[3]*a[1])
            -4
        """
        P = self.parent()
        if not self._p:
            return P.zero()
        self._p = P._P(self._p)
        if isinstance(monomial, dict):
            x = {P._P(v): d for v, d in monomial.items()}
        elif parent(monomial) is P:
            monomial._p = P._P(monomial._p)
            x = monomial._p
        else:
            raise TypeError(f"{monomial} must be a monomial in the parent of self")

        return InfinitePolynomial_dense(P, self._p.coefficient(x))
