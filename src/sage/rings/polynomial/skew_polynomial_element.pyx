r"""
Univariate Skew Polynomials

This module provides the
:class:`~sage.rings.polynomial.skew_polynomial_element.SkewPolynomial`,
which constructs a single univariate skew polynomial over commutative
base rings and an automorphism over the base ring. Skew polynomials are
non-commutative and so principal methods such as gcd, lcm, monic,
multiplication, and division are given in left and right forms.

The generic implementation of dense skew polynomials is
:class:`~sage.rings.polynomial.skew_polynomial_element.SkewPolynomial_generic_dense`.
The classes 
:class:`~sage.rings.polynomial.skew_polynomial_element.ConstantSkewPolynomialSection`
and :class:`~sage.rings.polynomial.skew_polynomial_element.SkewPolynomialBaseringInjection`
handle conversion from a skew polynomial ring to its base ring and vice versa respectively.

.. WARNING::

    The current semantics of
    :meth:`~sage.rings.polynomial.skew_polynomial_element.SkewPolynomial.__call__`
    are experimental, so a warning is thrown when a skew polynomial is evaluated
    for the first time in a session. See the method documentation for details.

    TESTS::

        sage: R.<t> = QQ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: a = 2*(t + x) + 1
        sage: a(t^2)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/13215 for details.
        2*t^3 + 3*t^2 + 4*t + 2
        sage: a(t)
        2*t^2 + 3*t + 2

AUTHORS:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and
  refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and
  doctest errors

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from __future__ import print_function, absolute_import, division

import re
from copy import copy
from sage.rings.infinity import infinity
from sage.structure.factorization import Factorization
from sage.structure.element cimport Element, RingElement, AlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from sage.structure.parent_gens cimport ParentWithGens
from sage.misc.abstract_method import abstract_method
from sage.categories.homset import Hom
from sage.categories.fields import Fields
from sage.rings.integer cimport Integer
from cpython.object cimport PyObject_RichCompare
from sage.categories.map cimport Map
from sage.rings.morphism cimport Morphism, RingHomomorphism
from sage.structure.element import coerce_binop
from sage.misc.superseded import experimental

cdef class SkewPolynomial(AlgebraElement):
    r"""
    Abstract base class for skew polynomials.

    This class must be inherited from and have key methods overridden.

    .. RUBRIC:: Definition

    Let `R` be a commutative ring equipped with an automorphism `\sigma`.

    Then, a skew polynomial is given by the equation:

    .. MATH::

        F(X) = a_{n} X^{n} + \cdots + a_0,

    where the coefficients `a_i \in R` and `X` is a formal variable.

    Addition between two skew polynomials is defined by the usual addition
    operation and the modified multiplication is defined by the rule
    `X a = \sigma(a) X` for all `a` in `R`. Skew polynomials are thus
    non-commutative and the degree of a product is equal to the sum of the
    degrees of the factors.

    Let `a` and `b` be two skew polynomials in the same ring `S`.
    The *left (resp. right) euclidean division* of `a` by `b` is a couple
    `(q,r)` of elements in `S` such that

    -  `a = q b + r` (resp. `a = b q + r`)

    -  the degree of `r` is less than the degree of `b`

    `q` (resp. `r`) is called the *quotient* (resp. the remainder)
    of this euclidean division.

    .. RUBRIC:: Properties

    Keeping the previous notation, if the leading coefficient of `b`
    is a unit (e.g. if `b` is monic) then the quotient and the remainder
    in the *right* euclidean division exist and are unique.

    The same result holds for the *left* euclidean division if in addition
    the twist map defining the skew polynomial ring is invertible.

    .. RUBRIC:: Evaluation

    The value of a given a skew polynomial `p(x) = \sum_{i=0}^d a_i x^i`
    at `r` is calculated using the formula:

    .. MATH::

        p(r) = \sum_{i=0}^d a_i \sigma^i(r)

    where `\sigma` is the base ring automorphism. This is called
    the *operator evaluation* method.

    EXAMPLES:

    We illustrate some functionalities implemented in this class.

    We create the skew polynomial ring::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    and some elements in it::

        sage: a = t + x + 1; a
        x + t + 1
        sage: b = S([t^2,t+1,1]); b
        x^2 + (t + 1)*x + t^2
        sage: c = S.random_element(degree=3,monic=True); c
        x^3 + (-95*t^2 + t + 2)*x^2 + (-t^2 + t)*x + 2*t - 8

    Ring operations are supported::

        sage: a + b
        x^2 + (t + 2)*x + t^2 + t + 1
        sage: a - b
        -x^2 - t*x - t^2 + t + 1

        sage: a * b
        x^3 + (2*t + 3)*x^2 + (2*t^2 + 4*t + 2)*x + t^3 + t^2
        sage: b * a
        x^3 + (2*t + 4)*x^2 + (2*t^2 + 3*t + 2)*x + t^3 + t^2
        sage: a * b == b * a
        False

        sage: b^2
        x^4 + (2*t + 4)*x^3 + (3*t^2 + 7*t + 6)*x^2
         + (2*t^3 + 4*t^2 + 3*t + 1)*x + t^4
        sage: b^2 == b*b
        True

    Sage also implements arithmetic over skew polynomial rings. You will find
    below a short panorama::

        sage: q,r = c.right_quo_rem(b)
        sage: q
        x - 95*t^2
        sage: r
        (95*t^3 + 93*t^2 - t - 1)*x + 95*t^4 + 2*t - 8
        sage: c == q*b + r
        True

    The operators ``//`` and ``%`` give respectively the quotient
    and the remainder of the *right* euclidean division::

        sage: q == c // b
        True
        sage: r == c % b
        True

    Left euclidean division won't work over our current `S` because Sage can't
    invert the twist map::

        sage: q,r = c.left_quo_rem(b)
        Traceback (most recent call last):
        ...
        NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
            Defn: t |--> t + 1

    Here we can see the effect of the operator evaluation compared to the usual
    polynomial evaluation::

        sage: a = x^2
        sage: a(t)
        t + 2

    Here is a working example over a finite field::

        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: a = x^4 + (4*t + 1)*x^3 + (t^2 + 3*t + 3)*x^2 + (3*t^2 + 2*t + 2)*x + (3*t^2 + 3*t + 1)
        sage: b = (2*t^2 + 3)*x^2 + (3*t^2 + 1)*x + 4*t + 2
        sage: q,r = a.left_quo_rem(b)
        sage: q
        (4*t^2 + t + 1)*x^2 + (2*t^2 + 2*t + 2)*x + 2*t^2 + 4*t + 3
        sage: r
        (t + 2)*x + 3*t^2 + 2*t + 4
        sage: a == b*q + r
        True

    Once we have euclidean divisions, we have for free gcd and lcm
    (at least if the base ring is a field)::

        sage: a = (x + t) * (x + t^2)^2
        sage: b = (x + t) * (t*x + t + 1) * (x + t^2)
        sage: a.right_gcd(b)
        x + t^2
        sage: a.left_gcd(b)
        x + t

    The left lcm has the following meaning: given skew polynomials `a` and `b`,
    their left lcm is the least degree polynomial `c = ua = vb` for some skew
    polynomials `u, v`. Such a `c` always exist if the base ring is a field::

        sage: c = a.left_lcm(b); c
        x^5 + (4*t^2 + t + 3)*x^4 + (3*t^2 + 4*t)*x^3 + 2*t^2*x^2 + (2*t^2 + t)*x + 4*t^2 + 4
        sage: c.is_right_divisible_by(a)
        True
        sage: c.is_right_divisible_by(b)
        True

    The right lcm is defined similarly as the least degree polynomial `c = au =
    bv` for some `u,v`::

        sage: d = a.right_lcm(b); d
        x^5 + (t^2 + 1)*x^4 + (3*t^2 + 3*t + 3)*x^3 + (3*t^2 + t + 2)*x^2 + (4*t^2 + 3*t)*x + 4*t + 4
        sage: d.is_left_divisible_by(a)
        True
        sage: d.is_left_divisible_by(b)
        True

    .. SEEALSO::

        - :mod:`sage.rings.polynomial.skew_polynomial_ring`
        - :mod:`sage.rings.polynomial.skew_polynomial_ring_constructor`
    """
    def __init__(self, parent, construct=False):
        r"""
        Initialize ``self``.

        INPUT:

        - ``parent`` -- parent of ``self``

        - ``construct`` -- boolean (default: ``False``)

        TESTS::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: P = x + t
            sage: TestSuite(P).run()
            sage: Q = S([1, t, t+2])
            sage: TestSuite(Q).run()
        """
        AlgebraElement.__init__(self, parent)

    cdef long _hash_c(self):
        raise NotImplementedError

    def __hash__(self):
        r"""
        Return hash of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: hash(a) == hash(a)
            True
        """
        return self._hash_c()

    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right):
        raise NotImplementedError
    cdef void _inplace_pow(self, Py_ssize_t n):
        raise NotImplementedError

    cpdef int degree(self):
        r"""
        Return the degree of ``self``.

        By convention, the zero skew polynomial has degree `-1`.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t*x^3 + t^2*x + 1
            sage: a.degree()
            3
            sage: S.zero().degree()
            -1
            sage: S(5).degree()
            0
        """

    cdef SkewPolynomial _new_c(self, list coeffs, Parent P, char check=0):
        r"""
        Fast creation of a new skew polynomial

        .. NOTE::

            Override this function in classes which inherit
            from SkewPolynomial.
        """
        l = P(coeffs)
        return l

    cpdef SkewPolynomial _new_constant_poly(self, RingElement a, Parent P, char check=0):
        r"""
        Fast creation of a new constant skew polynomial

        EXAMPLES::

            sage: from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: SkewPolynomialBaseringInjection(k, k['x', Frob]) #indirect doctest
            Skew Polynomial base injection morphism:
              From: Finite Field in t of size 5^3
              To:   Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        if a:
            n = self._new_c([a],P,check)
        else:
            n = self._new_c([],P)
        return n

    def __call__(self, eval_pt):
        r"""
        Evaluate ``self`` at ``eval_pt`` using operator evaluation.

        Given a skew polynomial `p(x) = \sum_{i=0}^d a_i * x^i`, we define
        the evaluation `p(r)` to be `\sum_{i=0}^d a_i * \sigma^i(r)`, where
        `\sigma` is the twist map of the skew polynomial ring.

        INPUT:

        - ``eval_pt`` -- element of the base ring of ``self``

        OUTPUT:

        The operator evaluation of ``self`` at ``eval_pt``.

        .. TODO::

            Currently, only "operator evaluation" of skew polynomials is
            implemented (see :meth:`.operator_eval`).
            There are two other notions of evaluation of a skew polynomial
            `p(x)` at some element `a` of the base ring. First, the value
            of the polynomial can be defined as the remainder of the right
            division of `p(x)` by `x-a`. Second, the value can be given by
            the formula, `p(a) = \sum_{i=0}^{m-1} B_{i} * p(\beta_{i})`
            where `m` is the degree of the base ring (`F_{q^m}`) of the skew
            polynomial ring, `B_{i}` is the `i`-th element in the vector
            representation of `a` in `F_{q}` and`\beta_{i}` is the `i`-th
            element of the corresponding basis of `F_{q^m}` over `F_{q}`.
            
            The current calling convention might change in the future to
            accommodate these. Therefore, the current method has been
            marked as experimental.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x + 1
            sage: a(t^2)
            t^3 + 3*t^2 + t
            sage: b = x^2 + t*x^3 + t^2*x + 1
            sage: b(2*t + 3)
            2*t^3 + 7*t^2 + 13*t + 10
        """
        return self._call(eval_pt)

    @experimental(trac_number=13215)
    def _call(self, eval_pt):
        r"""
        Helper function for the :meth:`__call__` method to accommodate
        the ``@experimental`` decorator.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x',Frob]
            sage: a = 3*t^2*x^2 + (t + 1)*x + 2
            sage: a(t) #indirect test
            2*t^2 + 2*t + 3
        """
        return self.operator_eval(eval_pt)

    cpdef operator_eval(self, eval_pt):
        r"""
        Evaluate ``self`` at ``eval_pt`` by the operator evaluation
        method.

        INPUT:

        - ``eval_pt`` -- element of the base ring of ``self``

        OUTPUT:

        The value of the polynomial at the point specified by the argument.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x',Frob]
            sage: a = 3*t^2*x^2 + (t + 1)*x + 2
            sage: a(t) #indirect test
            2*t^2 + 2*t + 3
            sage: a.operator_eval(t)
            2*t^2 + 2*t + 3

        Evaluation points outside the base ring is usually not possible due to the twist map::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x + 1
            sage: a.operator_eval(1/t)
            Traceback (most recent call last):
            ...
            TypeError: 1/t fails to convert into the map's domain Univariate Polynomial Ring in t over Rational Field, but a `pushforward` method is not properly implemented
        """
        cdef RingHomomorphism sigma = self._parent.twist_map()
        cdef list coefficients = self.list()
        cdef RingElement ret = self.base_ring().zero()
        cdef RingElement a = eval_pt
        for c in coefficients:
            ret += c * a
            a = sigma(a)
        return ret

    def __setitem__(self, n, value):
        r"""
        Set the ``n``-th coefficient of ``self``.

        This always raises an ``IndexError``, since polynomials are immutable in
        Sage.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: a[1] = t + 1
            Traceback (most recent call last):
            ...
            IndexError: skew polynomials are immutable
        """
        raise IndexError("skew polynomials are immutable")

    def square(self):
        r"""
        Return the square of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t; a
            x + t
            sage: a.square()
            x^2 + (2*t + 1)*x + t^2
            sage: a.square() == a*a
            True
        """
        return self * self

    def conjugate(self, n):
        r"""
        Return ``self`` conjugated by `x^n`, where `x` is the
        variable of ``self``.

        The conjugate is obtained from ``self`` by applying the `n`-th iterate
        of the twist map to each of its coefficients.

        INPUT:

        - `n` -- an integer, the power of conjugation

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^3 + (t^2 + 1)*x^2 + 2*t
            sage: b = a.conjugate(2); b
            (t + 2)*x^3 + (t^2 + 4*t + 5)*x^2 + 2*t + 4
            sage: x^2*a == b*x^2
            True

        In principle, negative values for `n` are allowed, but Sage needs to be
        able to invert the twist map::

            sage: b = a.conjugate(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t + 1

        Here is a working example::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: u = T.random_element(); u
            (2*t^2 + 3)*y^2 + (4*t^2 + t + 4)*y + 2*t^2 + 2
            sage: v = u.conjugate(-1); v
            (3*t^2 + t)*y^2 + (4*t^2 + 2*t + 4)*y + 3*t^2 + t + 4
            sage: u*y == y*v
            True
        """
        r = self._new_c([self._parent.twist_map(n)(x) for x in self.list()],
                        self._parent, 0)
        return r

    def constant_coefficient(self):
        r"""
        Return the constant coefficient (i.e. the coefficient of term
        of degree `0`) of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t^2 + 2
            sage: a.constant_coefficient()
            t^2 + 2
        """
        if not self:
            return self.base_ring().zero()
        else:
            return self[0]

    def leading_coefficient(self):
        r"""
        Return the coefficient of the highest-degree monomial of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (t+1)*x^5 + t^2*x^3 + x
            sage: a.leading_coefficient()
            t + 1
        """
        cdef int d = self.degree()
        if d == -1:
            raise ValueError("the skew polynomial must not be 0")
        return self[d]

    def is_unit(self):
        r"""
        Return ``True`` if this skew polynomial is a unit.

        When the base ring `R` is an integral domain, then a skew polynomial `f`
        is a unit if and only if degree of `f` is `0` and `f` is then a unit in
        `R`.

        .. NOTE::

            The case when `R` is not an integral domain is not yet implemented.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + (t+1)*x^5 + t^2*x^3 - x^5
            sage: a.is_unit()
            False
        """
        # TODO: Sage does not yet have support for finding order of
        #       automorphisms. Once that is available, general case can
        #       be implemented. Reference: http://bit.ly/29Vidu7
        if self._parent.base_ring().is_integral_domain():
            if self.degree() == 0 and self[0].is_unit():
                return True
            else:
                return False
        else:
            raise NotImplementedError("is_unit is not implemented for skew polynomial rings "
                                      "over base rings which are not integral domains.")

    def is_nilpotent(self):
        r"""
        Check if ``self`` is nilpotent.

        Given a commutative ring `R` and a base ring automorphism `\sigma`
        of order `n`, an element `f` of `R[X, \sigma]` is nilpotent if
        and only if all coefficients of `f^n` are nilpotent in `R`.

        .. NOTE::

            The paper "Nilpotents and units in skew polynomial rings
            over commutative rings" by M. Rimmer and K.R. Pearson describes
            the method to check whether a given skew polynomial is nilpotent.
            That method however, requires one to know the order of the
            automorphism which is not available in Sage. This method is thus
            not yet implemented. 

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.is_nilpotent()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_monic(self):
        r"""
        Return ``True`` if this skew polynomial is monic.

        The zero polynomial is by definition not monic.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: a.is_monic()
            True
            sage: a = 0*x
            sage: a.is_monic()
            False
            sage: a = t*x^3 + x^4 + (t+1)*x^2
            sage: a.is_monic()
            True
            sage: a = (t^2 + 2*t)*x^2 + x^3 + t^10*x^5
            sage: a.is_monic()
            False
        """
        return not self.is_zero() and self[self.degree()] == 1

    def left_monic(self):
        r"""
        Return the unique monic skew polynomial `m` which divides ``self`` on
        the left and has the same degree.

        Given a skew polynomial `p` of degree `n`, its left monic is given by
        `m = p \sigma^{-n}(1/k)`, where `k` is the leading coefficient of
        `p`, i.e. by the appropriate scalar multiplication on the right.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.left_monic(); b
            x^3 + (4*t^2 + 3*t)*x^2 + (4*t + 2)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: b.degree() == a.degree()
            True
            sage: b.is_left_divisible_by(a)
            True
            sage: twist = S.twist_map(-a.degree())
            sage: a == b * twist(a.leading_coefficient())
            True

        Note that `b` does not divide `a` on the right::

            sage: a.is_right_divisible_by(b)
            False

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x
            sage: a.left_monic()
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient is not a unit
        """
        try:
            a = self.base_ring()(~self.leading_coefficient())
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient is not a unit")
        r = self * self._parent.twist_map(-self.degree())(a)
        return r

    def right_monic(self):
        r"""
        Return the unique monic skew polynomial `m` which divides ``self`` on
        the right and has the same degree.

        Given a skew polynomial `p` of degree `n`, its left monic is given by
        `m = (1/k) * p`, where `k` is the leading coefficient of `p`, i.e. by
        the appropriate scalar multiplication on the left.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.right_monic(); b
            x^3 + (2*t^2 + 3*t + 4)*x^2 + (3*t^2 + 4*t + 1)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: b.degree() == a.degree()
            True
            sage: b.is_right_divisible_by(a)
            True
            sage: a == a.leading_coefficient() * b
            True

        Note that `b` does not divide `a` on the right::

            sage: a.is_left_divisible_by(b)
            False

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x
            sage: a.right_monic()
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient is not a unit
        """
        try:
            a = self.base_ring()(~self.leading_coefficient())
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient is not a unit")
        r = a * self
        return r

    cpdef _mod_(self, other):
        r"""
        Return the remainder in the *right* euclidean division of
        ``self`` by ``other```.

        TESTS::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: b = x^2 + 2*t*x + 2
            sage: a = (x+t)*b + t*x + 1
            sage: a % b
            t*x + 1

            sage: (a*t).right_quo_rem(b*t)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        _,r = self.right_quo_rem(other)
        return r

    cpdef _floordiv_(self, right):
        r"""
        Return the quotient of the *right* euclidean division of
        ``self`` by ``right``.

        The algorithm fails if the leading coefficient of the divisor
        (``right``) is not invertible.

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: b = x^2 + t
            sage: a = (x^2 + t*x + 1)*b + t^3*x
            sage: a // b
            x^2 + t*x + 1

            sage: (t*a) // (t*b)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible

        """
        q,_ = self.right_quo_rem(right)
        return q

    cpdef _div_(self, right):
        r"""
        Not Implemented.

        To implement this, localization of Ore rings is needed, see
        :trac:`13215`.

        Use the operator `//` even for exact division.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^5 + (t + 2)*x^2 + t^2
            sage: b = x^3 + 4*t
            sage: c = a*b

            sage: c / b
            Traceback (most recent call last):
            ...
            NotImplementedError: localization of Ore rings not yet implemented

            sage: c // b == a
            True
        """
        # Should this actually return something in the fraction field like
        #   we do elsewhere in Sage? - TCS
        raise NotImplementedError("localization of Ore rings not yet implemented")

    def is_left_divisible_by(self, other):
        r"""
        Check if ``self`` is divisible by ``other`` on the left.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a*b
            sage: c.is_left_divisible_by(a)
            True
            sage: c.is_left_divisible_by(b)
            False

        Divisibility by `0` does not make sense::

            sage: c.is_left_divisible_by(S(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero is not valid
        """
        _, r = self.left_quo_rem(other)
        return r.is_zero()

    def is_right_divisible_by(self, other):
        r"""
        Check if ``self`` is divisible by ``other`` on the right.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a*b
            sage: c.is_right_divisible_by(a)
            False
            sage: c.is_right_divisible_by(b)
            True

        Divisibility by `0` does not make sense::

            sage: c.is_right_divisible_by(S(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero is not valid

        This function does not work if the leading coefficient of the divisor
        is not a unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + 2*x + t
            sage: b = (t+1)*x + t^2
            sage: c = a*b
            sage: c.is_right_divisible_by(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        _, r = self.right_quo_rem(other)
        return r.is_zero()

    def left_divides(self, other):
        r"""
        Check if ``self`` divides ``other`` on the left.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a*b
            sage: a.left_divides(c)
            True
            sage: b.left_divides(c)
            False

        Divisibility by `0` does not make sense::

            sage: S(0).left_divides(c)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero is not valid
        """
        _, r = other.left_quo_rem(self)
        return r.is_zero()

    def right_divides(self, other):
        r"""
        Check if ``self`` divides ``other`` on the right.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a*b
            sage: a.right_divides(c)
            False
            sage: b.right_divides(c)
            True

        Divisibility by `0` does not make sense::

            sage: S(0).right_divides(c)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero is not valid

        This function does not work if the leading coefficient of the divisor
        is not a unit::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + 2*x + t
            sage: b = (t+1)*x + t^2
            sage: c = a*b
            sage: b.right_divides(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        _, r = other.right_quo_rem(self)
        return r.is_zero()

    @coerce_binop
    def left_xgcd(self, other, monic=True):
        r"""
        Return the left gcd of ``self`` and ``other`` along with the
        coefficients for the linear combination.

        If `a` is ``self`` and `b` is ``other``, then there are skew polynomials
        `u` and `v` such that `g = a u + b v`, where `g` is the left gcd of `a`
        and `b`. This method returns `(g, u, v)`.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``). Return whether the left gcd
          should be normalized to be monic.

        OUTPUT:

        - The left gcd of ``self`` and ``other``, that is a skew polynomial
          `g` with the following property: any skew polynomial is
          divisible on the left by `g` iff it is divisible on the left
          by both ``self`` and ``other``.
          If monic is ``True``, `g` is in addition monic. (With this
          extra condition, it is uniquely determined.)

        - Two skew polynomials `u` and `v` such that:

          .. MATH::

              g = a * u + b * v,

          where `s` is ``self`` and `b` is ``other``.

        .. NOTE::

            Works only if following two conditions are fulfilled
            (otherwise left gcd do not exist in general):
            1) the base ring is a field and 2) the twist map on
            this field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: g,u,v = a.left_xgcd(b); g
            x + t
            sage: a*u + b*v == g
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: g,u,v = a.left_xgcd(b, monic=False); g
            2*t*x + 4*t + 2
            sage: a*u + b*v == g
            True

        The base ring must be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_xgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map must be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_xgcd(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Fraction Field of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t^2
        """
        if self.base_ring() not in Fields:
            raise TypeError("the base ring must be a field")
        G = self
        U = self._parent.one()
        if other.is_zero():
            V = self._parent.zero()
        else:
            V1 = self._parent.zero()
            V3 = other
            while not V3.is_zero():
                Q,R = G.left_quo_rem(V3)
                T = U - V1*Q
                U = V1
                G = V3
                V1 = T
                V3 = R
            V, _ = (G - self*U).left_quo_rem(other)
        if monic:
            lc = ~G.leading_coefficient()
            lc = self._parent.twist_map(-G.degree())(lc)
            G = G * lc
            U = U * lc
            V = V * lc
        return G,U,V

    @coerce_binop
    def right_xgcd(self, other, monic=True):
        r"""
        Return the right gcd of ``self`` and ``other`` along with the
        coefficients for the linear combination.

        If `a` is ``self`` and `b` is ``other``, then there are skew polynomials
        `u` and `v` such that `g = u a + v b`, where `g` is the right gcd of `a`
        and `b`. This method returns `(g, u, v)`.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``). Return whether the right gcd
          should be normalized to be monic.

        OUTPUT:

        - The right gcd of ``self`` and ``other``, that is a skew polynomial
          `g` with the following property: any skew polynomial is
          divisible on the right by `g` iff it is divisible on the right
          by both ``self`` and ``other``.
          If monic is ``True``, `g` is in addition monic. (With this
          extra condition, it is uniquely determined.)

        - Two skew polynomials `u` and `v` such that:

          .. MATH::

              g = u * a + v * b

          where `a` is ``self`` and `b` is ``other``.

        .. NOTE::

            Works only if the base ring is a field (otherwise right
            gcd do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: g,u,v = a.right_xgcd(b); g
            x + t
            sage: u*a + v*b == g
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: g,u,v = a.right_xgcd(b,monic=False); g
            (4*t^2 + 4*t + 1)*x + 4*t^2 + 4*t + 3
            sage: u*a + v*b == g
            True

        The base ring must be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.right_xgcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
        """
        if self.base_ring() not in Fields:
            raise TypeError("the base ring must be a field")
        G = self
        U = self._parent.one()
        if other.is_zero():
            V = self._parent.zero()
        else:
            V1 = self._parent.zero()
            V3 = other
            while not V3.is_zero():
                Q, R = G.right_quo_rem(V3)
                T = U - Q*V1
                U = V1
                G = V3
                V1 = T
                V3 = R
            V,_ = (G - U*self).right_quo_rem(other)
        if monic:
            lc = ~G.leading_coefficient()
            G = lc * G
            U = lc * U
            V = lc * V
        return G,U,V

    @coerce_binop
    def right_gcd(self, other, monic=True):
        r"""
        Return the right gcd of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``). Return whether the right gcd
          should be normalized to be monic.

        OUTPUT:

        The right gcd of ``self`` and ``other``, that is a skew polynomial
        `g` with the following property: any skew polynomial is
        divisible on the right by `g` iff it is divisible on the right
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field (otherwise right
            gcd do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.right_gcd(b)
            x + t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.right_gcd(b,monic=False)
            (4*t^2 + 4*t + 1)*x + 4*t^2 + 4*t + 3

        The base ring need to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.right_gcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
        """
        if self.base_ring() not in Fields:
            raise TypeError("the base ring must be a field")
        if other.is_zero():
            return self
        A = self
        B = other
        while not B.is_zero():
            A, B = B, A % B
        if monic:
            A = A.right_monic()
        return A

    @coerce_binop
    def left_gcd(self, other, monic=True):
        r"""
        Return the left gcd of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``). Return whether the left gcd
          should be normalized to be monic.

        OUTPUT:

        The left gcd of ``self`` and ``other``, that is a skew polynomial
        `g` with the following property: any skew polynomial is
        divisible on the left by `g` iff it is divisible on the left
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if two following conditions are fulfilled
            (otherwise left gcd do not exist in general):
            1) the base ring is a field and 2) the twist map on
            this field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_gcd(b)
            x + t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.left_gcd(b,monic=False)
            2*t*x + 4*t + 2

        The base ring needs to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_gcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map needs to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_gcd(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Fraction Field of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t^2
        """
        if self.base_ring() not in Fields:
            raise TypeError("the base ring must be a field")
        if other.is_zero():
            return self
        A = self
        B = other
        while not B.is_zero():
            A_ = A
            A = B
            _, B = A_.left_quo_rem(B)
        if monic:
            A = A.left_monic()
        return A

    @coerce_binop
    def left_lcm(self, other, monic=True):
        r"""
        Return the left lcm of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``). Return whether the left lcm
          should be normalized to be monic.

        OUTPUT:

        The left lcm of ``self`` and ``other``, that is a skew polynomial
        `g` with the following property: any skew polynomial divides
        `g` on the *right* iff it divides both ``self`` and ``other``
        on the *right*.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field (otherwise left
            lcm do not exist in general).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t^2) * (x + t)
            sage: b = 2 * (x^2 + t + 1) * (x * t)
            sage: c = a.left_lcm(b); c
            x^5 + (2*t^2 + t + 4)*x^4 + (3*t^2 + 4)*x^3 + (3*t^2 + 3*t + 2)*x^2 + (t^2 + t + 2)*x
            sage: c.is_right_divisible_by(a)
            True
            sage: c.is_right_divisible_by(b)
            True
            sage: a.degree() + b.degree() == c.degree() + a.right_gcd(b).degree()
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.left_lcm(b,monic=False)
            (t^2 + t)*x^5 + (4*t^2 + 4*t + 1)*x^4 + (t + 1)*x^3 + (t^2 + 2)*x^2 + (3*t + 4)*x

        The base ring needs to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t^2) * (x + t)
            sage: b = 2 * (x^2 + t + 1) * (x * t)
            sage: a.left_lcm(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field
        """
        if self.base_ring() not in Fields:
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError("division by zero is not valid")
        U = self._parent.one()
        G = self
        V1 = self._parent.zero()
        V3 = other
        while not V3.is_zero():
            Q, R = G.right_quo_rem(V3)
            T = U - Q*V1
            U = V1
            G = V3
            V1 = T
            V3 = R
        V1 = V1 * self
        if monic:
            V1 = V1.right_monic()
        return V1

    @coerce_binop
    def right_lcm(self, other, monic=True):
        r"""
        Return the right lcm of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``
        
        - ``monic`` -- boolean (default: ``True``). Return whether the right lcm
          should be normalized to be monic.

        OUTPUT:

        The right lcm of ``self`` and ``other``, that is a skew polynomial
        `g` with the following property: any skew polynomial divides
        `g` on the *left* iff it divides both ``self`` and ``other``
        on the *left*.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if two following conditions are fulfilled
            (otherwise right lcm do not exist in general):
            1) the base ring is a field and 2) the twist map on
            this field is bijective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: c = a.right_lcm(b); c
            x^4 + (2*t^2 + t + 2)*x^3 + (3*t^2 + 4*t + 1)*x^2 + (3*t^2 + 4*t + 1)*x + t^2 + 4
            sage: c.is_left_divisible_by(a)
            True
            sage: c.is_left_divisible_by(b)
            True
            sage: a.degree() + b.degree() == c.degree() + a.left_gcd(b).degree()
            True

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.right_lcm(b,monic=False)
            2*t*x^4 + (3*t + 1)*x^3 + (4*t^2 + 4*t + 3)*x^2
             + (3*t^2 + 4*t + 2)*x + 3*t^2 + 2*t + 3

        The base ring needs to be a field::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: a.right_lcm(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twist map needs to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: a.right_lcm(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Fraction Field of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t^2
        """
        if self.base_ring() not in Fields:
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError("division by zero is not valid")
        R = self.parent()
        U = R.one()
        G = self
        V1 = R.zero()
        V3 = other
        while not V3.is_zero():
            Q, R = G.left_quo_rem(V3)
            T = U - V1*Q
            U = V1
            G = V3
            V1 = T
            V3 = R
        V1 = self * V1
        if monic:
            V1 = V1.left_monic()
        return V1

    def _repr_(self, name=None):
        r"""
        Return string representation of this skew polynomial.

        INPUT:

        - ``name`` -- the name of the variable (default: the
          name given when the skew polynomial ring was created)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t;
            sage: a._repr_()
            '(1/2*t + 1/2)*x + t^2'
            sage: a._repr_(name='y')
            '(1/2*t + 1/2)*y + t^2'
        """
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = self.list()
        for n in reversed(range(m)):
            x = coeffs[n]
            if x:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*',' ', s)
        s = re.sub(r' -1(\.0+)?\*',' -', s)
        if s == " ":
            return "0"
        return s[1:]

    def _latex_(self, name=None):
        r"""
        Return a latex representation of this skew polynomial.

        INPUT:

        - ``name`` -- the name of the variable (default: the
          name given when the skew polynomial ring was created)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t;
            sage: a._latex_()
            '\\left(\\frac{1}{2} t + \\frac{1}{2}\\right) x + t^{2}'
            sage: a._latex_(name='y')
            '\\left(\\frac{1}{2} t + \\frac{1}{2}\\right) y + t^{2}'
        """
        s = " "
        coeffs = self.list()
        m = len(coeffs)
        if name is None:
            name = self.parent().latex_variable_names()[0]
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        for n in reversed(range(m)):
            x = self[n]
            x = y = x._latex_()
            if x != '0':
                if n != m-1:
                    s += " + "
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "|%s^{%s}"%(name,n)
                elif n==1:
                    var = "|%s"%name
                else:
                    var = ""
                s += "%s %s"%(x,var)
        s = s.replace(" + -", " - ")
        s = re.sub(" 1(\.0+)? \|"," ", s)
        s = re.sub(" -1(\.0+)? \|", " -", s)
        s = s.replace("|","")
        if s == " ":
            return "0"
        return s[1:].lstrip().rstrip()

    def _is_atomic(self):
        r"""
        Check ``self`` is a single monomial whose leading coefficient
        is atomic in the base ring.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S([t+1])._is_atomic()
            False
            sage: S([1])._is_atomic()
            True
        """
        return (self.degree() == self.valuation() and
                self.leading_coefficient()._is_atomic())

    def __nonzero__(self):
        r"""
        Test whether ``self`` is nonzero.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + 1
            sage: bool(a)
            True
            sage: b = S.zero()
            sage: bool(b)
            False
        """
        return not self.is_zero()

    def base_ring(self):
        r"""
        Return the base ring of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element()
            sage: a.base_ring()
            Univariate Polynomial Ring in t over Integer Ring
            sage: a.base_ring() is R
            True
        """
        return self.parent().base_ring()

    def shift(self, n):
        r"""
        Return ``self`` multiplied on the right by the power `x^n`.

        If `n` is negative, terms below `x^n` will be discarded.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^5 + t^4*x^4 + t^2*x^2 + t^10
            sage: a.shift(0)
            x^5 + t^4*x^4 + t^2*x^2 + t^10
            sage: a.shift(-1)
            x^4 + t^4*x^3 + t^2*x
            sage: a.shift(-5)
            1
            sage: a.shift(2)
            x^7 + t^4*x^6 + t^2*x^4 + t^10*x^2

        One can also use the infix shift operator::

            sage: a >> 2
            x^3 + t^4*x^2 + t^2
            sage: a << 2
            x^7 + t^4*x^6 + t^2*x^4 + t^10*x^2
        """
        if n == 0 or self.degree() < 0:
            return self
        if n > 0:
            return self._parent(n*[self.base_ring().zero()] + self.list(), check=False)
        if n < 0:
            if n > self.degree():
                return self._parent([])
            else:
                return self._parent(self.list()[-n:], check=False)

    def __lshift__(self, k):
        r"""
        Return ``self`` multiplied on the right by the power `x^k`.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^5 + t^4*x^4 + t^2*x^2 + t^10
            sage: a << 2
            x^7 + t^4*x^6 + t^2*x^4 + t^10*x^2
        """
        return self.shift(k)

    def __rshift__(self, k):
        r"""
        Return ``self`` multiplied on the right by the power `x^(-k)`.

        If `n` is negative, terms below `x^n` will be discarded.
        
        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^5 + t^4*x^4 + t^2*x^2 + t^10
            sage: a >> 2
            x^3 + t^4*x^2 + t^2
        """
        return self.shift(-k)

    def change_variable_name(self, var):
        r"""
        Change the name of the variable of ``self``.

        This will create the skew polynomial ring with the new name but same
        base ring and twist map. The returned skew polynomial will be an element
        of that skew polynomial ring.

        INPUT:

        - ``var`` -- the name of the new variable

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x', sigma]
            sage: a = x^3 + (2*t + 1)*x  + t^2 + 3*t + 5
            sage: b = a.change_variable_name('y'); b
            y^3 + (2*t + 1)*y  + t^2 + 3*t + 5

        Note that a new parent is created at the same time::

            sage: b.parent()
            Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring
             twisted by t |--> t + 1
        """
        parent = self._parent
        R = parent.base_ring()[var,parent.twist_map()]
        return R(self.list())

    def is_term(self):
        r"""
        Return ``True`` if ``self`` is an element of the base ring times a
        power of the generator.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.is_term()
            True
            sage: R(1).is_term()
            True
            sage: (3*x^5).is_term()
            True
            sage: (1+3*x^5).is_term()
            False

        If you want to test that ``self`` also has leading coefficient 1, use
        :meth:`is_monomial()` instead::

            sage: (3*x^5).is_monomial()
            False
        """
        return len(self.exponents()) == 1

    def is_monomial(self):
        r"""
        Return ``True`` if ``self`` is a monomial, i.e., a power of
        the generator.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.is_monomial()
            True
            sage: (x+1).is_monomial()
            False
            sage: (x^2).is_monomial()
            True
            sage: S(1).is_monomial()
            True

        The coefficient must be 1::

            sage: (2*x^5).is_monomial()
            False
            sage: S(t).is_monomial()
            False

        To allow a non-1 leading coefficient, use is_term()::

            sage: (2*x^5).is_term()
            True
            sage: S(t).is_term()
            True
        """
        return self.is_term() and self.leading_coefficient() == 1

    cpdef list coefficients(self, sparse=True):
        r"""
        Return the coefficients of the monomials appearing in ``self``.

        If ``sparse=True`` (the default), return only the non-zero coefficients.
        Otherwise, return the same value as ``self.list()``.

        .. NOTE::

            This should be overridden in subclasses.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.coefficients()
            [t^2 + 1, t + 1, 1]
            sage: a.coefficients(sparse=False)
            [t^2 + 1, 0, t + 1, 0, 1]
        """
        raise NotImplementedError

    def number_of_terms(self):
        r"""
        Return the number of non-zero coefficients of ``self``.

        This is also known as the weight, hamming weight or sparsity.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.number_of_terms()
            3

        This is also an alias for ``hamming_weight``::

            sage: a.hamming_weight()
            3
        """
        return len(self.coefficients())

    # alias hamming_weight for number_of_terms:
    hamming_weight = number_of_terms

    def __copy__(self):
        r"""
        Return a "copy" of ``self``.

        In Sage, since skew polynomials are immutable, this just returns
        ``self`` again.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: b = copy(a)
            sage: b is a
            True
        """
        return self

    cpdef bint is_zero(self):
        r"""
        Return ``True`` if ``self`` is the zero polynomial.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + 1
            sage: a.is_zero()
            False
            sage: b = S.zero()
            sage: b.is_zero()
            True
        """
        return self.degree() == -1

    cpdef bint is_one(self):
        r"""
        Test whether this polynomial is `1`.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: R(1).is_one()
            True
            sage: (x + 3).is_one()
            False
        """
        return self.degree() == 0 and self[0].is_one()

    @coerce_binop
    def right_mod(self, other):
        r"""
        Return the remainder of right division of ``self`` by ``other``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + t*x^2
            sage: b = x + 1
            sage: a % b
            t + 1
            sage: (x^3 + x - 1).right_mod(x^2 - 1)
            2*x - 1
        """
        return self % other

    @coerce_binop
    def left_mod(self, other):
        r"""
        Return the remainder of left division of ``self`` by ``other``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = 1 + t*x^2
            sage: b = x + 1
            sage: a.left_mod(b)
            2*t^2 + 4*t
        """
        (_,r) = self.left_quo_rem(other)
        return r

    def is_constant(self):
        r"""
        Return whether ``self`` is a constant polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: R(2).is_constant()
            True
            sage: (x + 1).is_constant()
            False
        """
        return self.degree() <= 0

    def exponents(self):
        r"""
        Return the exponents of the monomials appearing in ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.exponents()
            [0, 2, 4]
        """
        return [i for i in range(self.degree()+1) if bool(self[i])]

    def prec(self):
        r"""
        Return the precision of ``self``.

        This is always infinity, since polynomials are of infinite precision by
        definition (there is no big-oh).

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: x.prec()
            +Infinity
        """
        return infinity

    def padded_list(self, n=None):
        r"""
        Return list of coefficients of ``self`` up to (but not including)
        degree `n`.

        Includes `0`s in the list on the right so that the list always has length
        exactly `n`.

        INPUT:

        - ``n`` -- (default: ``None``); if given, an integer that
          is at least `0`

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + t*x^3 + t^2*x^5
            sage: a.padded_list()
            [1, 0, 0, t, 0, t^2]
            sage: a.padded_list(10)
            [1, 0, 0, t, 0, t^2, 0, 0, 0, 0]
            sage: len(a.padded_list(10))
            10
            sage: a.padded_list(3)
            [1, 0, 0]
            sage: a.padded_list(0)
            []
            sage: a.padded_list(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be at least 0
        """
        v = self.list()
        if n is None:
            return v
        if n < 0:
            raise ValueError("n must be at least 0")
        if len(v) < n:
            z = self._parent.base_ring().zero()
            return v + [z]*(n - len(v))
        else:
            return v[:int(n)]

    def variable_name(self):
        r"""
        Return the string name of the variable used in ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: a.variable_name()
            'x'
        """
        return self.parent().variable_name()


cdef class SkewPolynomial_generic_dense(SkewPolynomial):
    r"""
    Generic implementation of dense skew polynomial supporting any valid base
    ring and twist map.
    """

    def __init__(self, parent, x=None, int check=1, int construct=0, **kwds):
        r"""
        Construct a skew polynomial over the given parent with the given
        coefficients.

        INPUT:

        - ``parent`` -- parent of ``self``

        - ``x`` -- list of coefficients from which ``self`` can be constructed

        - ``check`` -- flag variable to normalize the polynomial

        - ``construct`` -- boolean (default: ``False``)

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]

        We create a skew polynomial from a list::

            sage: S([t,1])
            x + t

        from another skew polynomial::

            sage: S(x^2 + t)
            x^2 + t

        from a constant::

            sage: x = S(t^2 + 1); x
            t^2 + 1
            sage: x.parent() is S
            True
        """
        SkewPolynomial.__init__(self, parent)
        if x is None:
            self._coeffs = []
            return

        R = parent.base_ring()
        if isinstance(x, list):
            if check:
                self._coeffs = [R(t) for t in x]
                self.__normalize()
            else:
                self._coeffs = x
            return

        if isinstance(x, SkewPolynomial):
            if (<Element>x)._parent is self._parent:
                x = list(x.list())
            elif R.has_coerce_map_from((<Element>x)._parent):
                try:
                    if x.is_zero():
                        self._coeffs = []
                        return
                except (AttributeError, TypeError):
                    pass
                x = [x]
            else:
                self._coeffs = [R(a, **kwds) for a in x.list()]
                if check:
                    self.__normalize()
                return

        elif isinstance(x, int) and x == 0:
            self._coeffs = []
            return

        elif isinstance(x, dict):
            x = self._dict_to_list(x, R.zero())

        elif not isinstance(x, list):
            x = [x]
        if check:
            self._coeffs = [R(z, **kwds) for z in x]
            self.__normalize()
        else:
            self._coeffs = x

    def __reduce__(self):
        r"""
        Return the generic dense skew polynomial corresponding to the
        current parameters provided ``self``.

        EXAMPLES:

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: loads(dumps(x)) == x
            True
            sage: loads(dumps(x))
            x
        """
        return (self._parent, (self._coeffs,))

    cdef long _hash_c(self):
        r"""
        This hash incorporates the name of the variable.

        .. NOTE::

            This is an internal method. Use :meth:`__hash__` instead.
        """
        #todo - come up with a way to create hashes of zero that
        #       that do not incorrectly indicate that the element is 0.
        cdef long result = 0
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash = 0
        cdef int i
        for i from 0 <= i < len(self._coeffs):
            if i == 1:
                var_name_hash = hash((<ParentWithGens>self._parent)._names[0])
            c_hash = hash(self._coeffs[i])
            if c_hash != 0:
                if i == 0:
                    result += c_hash
                else:
                    result_mon = c_hash
                    result_mon = (1000003 * result_mon) ^ var_name_hash
                    result_mon = (1000003 * result_mon) ^ i
                    result += result_mon
        if result == -1:
            return -2
        return result

    cpdef _richcmp_(left, right, int op):
        r"""
        Compare the two skew polynomials ``self`` and ``other``.

        We order polynomials first by degree, then in dictionary order
        starting with the coefficient of largest degree.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: b = (2*t^2)*x + t + 1
            sage: a > b
            True
            sage: a < b
            False
        """
        cdef x = (<SkewPolynomial_generic_dense>left)._coeffs
        cdef y = (<SkewPolynomial_generic_dense>right)._coeffs
        return PyObject_RichCompare(x, y, op)

    def __iter__(self):
        r"""
        Iterate over the list of coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: P = S([1, 2, 3])
            sage: [y for y in iter(P)]
            [1, 2, 3]
        """
        return iter((<SkewPolynomial_generic_dense>self)._coeffs)

    def __getitem__(self, n):
        r"""
        Return the `n`-th coefficient of ``self``.

        INPUT:

        - ``n`` -- an integer

        OUTPUT:

        - the ``n``-th coefficient of ``self``

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^2 + (t + 3/7)*x + t^2
            sage: a[1]
            t + 3/7
            sage: a[3]
            0
        """
        try:
            l = (<SkewPolynomial_generic_dense>self)._coeffs[n]
            return l
        except IndexError:
            return self.base_ring().zero()

    cpdef list list(self):
        r"""
        Return a list of the coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: l = a.list(); l
            [t^2 + 1, 0, t + 1, 0, 1]

        Note that `l` is a list, it is mutable, and each call to the list
        method returns a new list::

            sage: type(l)
            <... 'list'>
            sage: l[0] = 5
            sage: a.list()
            [t^2 + 1, 0, t + 1, 0, 1]
        """
        return list((<SkewPolynomial_generic_dense>self)._coeffs) # This creates a shallow copy

    cpdef dict dict(self):
        r"""
        Return a dictionary representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2012 + t*x^1006 + t^3 + 2*t
            sage: a.dict()
            {0: t^3 + 2*t, 1006: t, 2012: 1}
        """
        cdef dict X = {}
        cdef list Y = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef int i
        for i in range(len(Y)):
            c = Y[i]
            if c:
                X[i] = c
        return X

    cpdef int degree(self):
        r"""
        Return the degree of ``self``.

        By convention, the zero skew polynomial has degree `-1`.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t*x^3 + t^2*x + 1
            sage: a.degree()
            3

        By convention, the degree of `0` is `-1`::

            sage: S(0).degree()
            -1
        """
        return len(self._coeffs) - 1

    cpdef _add_(self, right):
        r"""
        Add two polynomials.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(monic=True); a
            x^2 + (-12*t^2 + 1/2*t - 1/95)*x - 1/2*t^2 - 4
            sage: b = -S.random_element(monic=True); b
            -x^2 + (5/2*t - 2/3)*x + 1/4*t^2 - 1/2*t + 1
            sage: c = a+b; c
            (-12*t^2 + 3*t - 193/285)*x - 1/4*t^2 - 1/2*t - 3
            sage: c.degree()
            1
        """
        cdef Py_ssize_t i, min
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef list y = (<SkewPolynomial_generic_dense>right)._coeffs
        cdef Py_ssize_t dx = len(x), dy = len(y)

        if dx > dy:
            r = self._new_c([x[i] + y[i] for i from 0 <= i < dy] + x[dy:], self._parent, 0)
        elif dx < dy:
            r = self._new_c([x[i] + y[i] for i from 0 <= i < dx] + y[dx:], self._parent, 0)
        else:
            r = self._new_c([x[i] + y[i] for i in range(dx)], self._parent, 1)
        return r

    cpdef _sub_(self, right):
        r"""
        Subtract polynomial ``right`` from ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(monic=True); a
            x^2 + (-12*t^2 + 1/2*t - 1/95)*x - 1/2*t^2 - 4
            sage: b = S.random_element(monic=True); b
            x^2 + (-5/2*t + 2/3)*x - 1/4*t^2 + 1/2*t - 1
            sage: c = a-b; c
            (-12*t^2 + 3*t - 193/285)*x - 1/4*t^2 - 1/2*t - 3
            sage: c.degree()
            1
        """
        cdef Py_ssize_t i, min
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef list y = (<SkewPolynomial_generic_dense>right)._coeffs
        cdef Py_ssize_t dx = len(x), dy = len(y)
        cdef RingElement c

        if dx > dy:
            r = self._new_c([x[i] - y[i] for i from 0 <= i < dy] + x[dy:], self._parent, 0)
        elif dx < dy:
            r = self._new_c([x[i] - y[i] for i from 0 <= i < dx] + [ -c for c in y[dx:] ], self._parent, 0)
        else:
            r = self._new_c([x[i] - y[i] for i from 0 <= i < dx], self._parent, 1)
        return r

    cpdef _neg_(self):
        r"""
        Return the negative of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^2 + x - 3
            sage: -a
            -t*x^2 - x + 3
        """
        c = self._new_c([-x for x in (<SkewPolynomial_generic_dense>self)._coeffs],
                        self._parent, 0)
        return c

    cpdef ModuleElement _lmul_(self, RingElement right):
        r"""
        Multiply ``self`` on the right by scalar.

        INPUT:

        - ``right`` -- an element of the base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: b = t
            sage: a * b
            (t + 1)*x + t^2
            sage: a * b == b * a
            False
        """
        if right == 0:
            return self._parent.zero()
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef Py_ssize_t i
        twist_map = self._parent._map
        r = self._new_c([ (twist_map**i)(right)*x[i] for i from 0 <= i < len(x) ],
                        self._parent, 0)
        return r

    cpdef ModuleElement _rmul_(self, RingElement left):
        r"""
        Multiply ``self`` on the left by scalar.

        INPUT:

        - ``left`` -- an element of the base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t
            sage: b = x + t
            sage: a * b
            t*x + t^2
            sage: a * b == b * a
            False
        """
        if left == 0:
            return self.parent().zero()
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef Py_ssize_t i
        r = self._new_c([ left*x[i] for i from 0 <= i < len(x) ], self._parent, 0)
        return r

    cpdef _mul_(self, right):
        r"""
        Multiply ``self`` on the right by a skew polynomial.

        INPUT:

        - ``right`` -- a skew polynomial in the same ring as ``self``

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t; a
            x^2 + t
            sage: b = x^2 + (t + 1)*x; b
            x^2 + (t + 1)*x
            sage: a * b
            x^4 + (t + 3)*x^3 + t*x^2 + (t^2 + t)*x
            sage: a * b == b * a
            False

        TESTS::

            sage: S(0)*a, (S(0)*a).list()
            (0, [])
        """
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef list y = (<SkewPolynomial_generic_dense>right)._coeffs
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t dx = len(x)-1, dy = len(y)-1
        parent = self._parent
        if dx == -1:
            return self # = zero
        elif dy == -1:
            return right # = zero
        elif dx == 0:
            c = x[0]
            r = self._new_c([c*a for a in y], parent, 0)
            return r
        cdef list coeffs = []
        for k from 0 <= k <= dx+dy:
            start = 0 if k <= dy else k-dy
            end = k if k <= dx else dx
            sum = x[start] * parent.twist_map(start)(y[k-start])
            for i from start < i <= end:
                sum += x[i] * parent.twist_map(i)(y[k-i])
            coeffs.append(sum)
        r = self._new_c(coeffs, parent, 0)
        return r

    cdef SkewPolynomial _new_c(self, list coeffs, Parent P, char check=0):
        r"""
        Fast creation of a new skew polynomial given a list of coefficients.

        .. WARNING::

            The list ``coeffs`` is stored internally in the newly created skew
            polynomial, so this must not be modified after calling this method.

        TESTS::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^3 + x^4 + (t+1)*x^2
            sage: a.truncate(4) #indirect doctest
            t*x^3 + (t + 1)*x^2
        """
        cdef type t = type(self)
        cdef SkewPolynomial_generic_dense f = t.__new__(t)
        f._parent = P
        f._coeffs = coeffs
        if check:
            f.__normalize()
        return f

    cdef void __normalize(self):
        r"""
        Remove higher order `0`-coefficients from the representation of ``self``.

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S #indirect doctest
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
        """
        cdef list x = self._coeffs
        cdef Py_ssize_t n = len(x) - 1
        while n >= 0 and not x[n]:
            del x[n]
            n -= 1

    def valuation(self):
        r"""
        Return the minimal degree of a non-zero monomial of ``self``.

        By convention, the zero skew polynomial has valuation `+\infty`.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t*x^3 + t^2*x
            sage: a.valuation()
            1

        By convention, the valuation of `0` is `+\infty`::

            sage: S(0).valuation()
            +Infinity
        """
        cdef list x = self._coeffs
        if not x:
            return infinity
        cdef Py_ssize_t v = 0
        while x[v].is_zero() and v < len(x):
            v += 1
        return v

    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right):
        r"""
        Replace ``self`` by ``self*right`` (only for internal use).

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: a.left_power_mod(100,modulus)  # indirect doctest
            (4*t^2 + t + 1)*x^2 + (t^2 + 4*t + 1)*x + 3*t^2 + 3*t
        """
        cdef list x = self._coeffs
        cdef list y = right._coeffs
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t d1 = len(x)-1, d2 = len(y)-1
        parent = self._parent
        if d2 == -1:
            self._coeffs = [ ]
        elif d1 >= 0:
            for k from d1 < k <= d1+d2:
                start = 0 if k <= d2 else k-d2
                sum = x[start] * parent.twist_map(start)(y[k-start])
                for i from start < i <= d1:
                    sum += x[i] * parent.twist_map(i)(y[k-i])
                x.append(sum)
            for k from d1 >= k >= 0:
                start = 0 if k <= d2 else k-d2
                end = k if k <= d1 else d1
                sum = x[start] * parent.twist_map(start)(y[k-start])
                for i from start < i <= end:
                    sum += x[i] * parent.twist_map(i)(y[k-i])
                x[k] = sum

    cdef void _inplace_pow(self, Py_ssize_t n):
        r"""
        Replace ``self`` by ``self**n`` (only for internal use).

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: a.left_power_mod(100,modulus)  # indirect doctest
            (4*t^2 + t + 1)*x^2 + (t^2 + 4*t + 1)*x + 3*t^2 + 3*t
        """
        while n & 1 == 0:
            self._inplace_rmul(self)
            n = n >> 1
        cdef SkewPolynomial_generic_dense selfpow = <SkewPolynomial_generic_dense>self._new_c(list(self._coeffs), self._parent)
        n = n >> 1
        while n != 0:
            selfpow._inplace_rmul(selfpow)
            if n&1 == 1:
                self._inplace_rmul(selfpow)
            n = n >> 1

    def truncate(self, n):
        r"""
        Return the polynomial resulting from discarding all monomials of degree
        at least `n`.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x^3 + x^4 + (t+1)*x^2
            sage: a.truncate(4)
            t*x^3 + (t + 1)*x^2
            sage: a.truncate(3)
            (t + 1)*x^2
        """
        return self._new_c(self._coeffs[:n], self._parent, 1)

    @coerce_binop
    def left_quo_rem(self, other):
        r"""
        Return the quotient and remainder of the left euclidean
        division of ``self`` by ``other``.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        - the quotient and the remainder of the left euclidean
          division of this skew polynomial by ``other``

        .. NOTE::

            This will fail if the leading coefficient of ``other`` is not a unit
            or if Sage can't invert the twist map.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = (3*t^2 + 4*t + 2)*x^2 + (2*t^2 + 4*t + 3)*x + 2*t^2 + t + 1
            sage: q,r = a.left_quo_rem(b)
            sage: a == b*q + r
            True

        In the following example, Sage does not know the inverse
        of the twist map::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = (-2*t^2 - t + 1)*x^3 + (-t^2 + t)*x^2 + (-12*t - 2)*x - t^2 - 95*t + 1
            sage: b = x^2 + (5*t - 6)*x - 4*t^2 + 4*t - 1
            sage: a.left_quo_rem(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring
                Defn: t |--> t + 1
        """
        cdef list a = list(self._coeffs)
        cdef list b = (<SkewPolynomial_generic_dense?>other)._coeffs
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = self.degree(), db = other.degree()
        if db < 0:
            raise ZeroDivisionError("division by zero is not valid")
        if da < db:
            return (self._new_c([], self._parent), self)
        try:
            inv = self.base_ring()(~b[db])
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient of the divisor is not invertible")
        cdef list q = [ ]
        parent = self._parent
        for i from da-db >= i >= 0:
            try:
                c = parent.twist_map(-db)(inv*a[i+db])
                for j from 0 <= j < db:
                    a[i+j] -= b[j] * parent.twist_map(j)(c)
            except Exception:
                raise NotImplementedError("inversion of the twist map %s" % parent.twist_map())
            q.append(c)
        q.reverse()
        return (self._new_c(q, parent), self._new_c(a[:db], parent, 1))

    @coerce_binop
    def right_quo_rem(self, other):
        r"""
        Return the quotient and remainder of the right euclidean
        division of ``self`` by ``other``.

        INPUT:

        - ``other`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        - the quotient and the remainder of the left euclidean
          division of this skew polynomial by ``other``

        .. NOTE::

            This will fail if the leading coefficient of the divisor
            is not a unit.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(degree=4); a
            t^2*x^4 + (-12*t^2 - 2*t - 1)*x^3 + (-95*t^2 + t + 2)*x^2 + (-t^2 + t)*x + 2*t - 8
            sage: b = S.random_element(monic=True); b
            x^2 + (4*t^2 - t - 2)*x - t^2 + t - 1
            sage: q,r = a.right_quo_rem(b)
            sage: a == q*b + r
            True

        The leading coefficient of the divisor need to be invertible::

            sage: c = S.random_element(); c
            (-4*t^2 + t)*x^2 - 2*t^2*x + 5*t^2 - 6*t - 4
            sage: a.right_quo_rem(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        cdef list a = list(self._coeffs)
        cdef list b = (<SkewPolynomial_generic_dense?>other)._coeffs
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = self.degree(), db = other.degree()
        parent = self._parent
        if db < 0:
            raise ZeroDivisionError("division by zero is not valid")
        if da < db:
            return (self._new_c([],parent), self)
        try:
            inv = self.base_ring()(~b[db])
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient of the divisor"
                                      " is not invertible")
        cdef list q = [ ]
        parent = self._parent
        for i from da-db >= i >= 0:
            c = parent.twist_map(i)(inv) * a[i+db]
            for j from 0 <= j < db:
                a[i+j] -= c * parent.twist_map(i)(b[j])
            q.append(c)
        q.reverse()
        return (self._new_c(q, parent), self._new_c(a[:db], parent, 1))

    cpdef left_power_mod(self, exp, modulus):
        r"""
        Return the remainder of ``self**exp`` in the left euclidean division
        by ``modulus``.

        INPUT:

        - ``exp`` -- an Integer

        - ``modulus`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Remainder of ``self**exp`` in the left euclidean division
        by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: a.left_power_mod(100,modulus)
            (4*t^2 + t + 1)*x^2 + (t^2 + 4*t + 1)*x + 3*t^2 + 3*t
        """
        cdef SkewPolynomial_generic_dense r
        if not isinstance(exp, Integer):
            try:
                exp = Integer(exp)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        if len(self._coeffs) <= 1:
            return self.parent()(self._coeffs[0]**exp)
        if exp == 0:
            return self.parent().one()
        if exp < 0:
            return (~self).left_power_mod(-exp, modulus)

        if self == self.parent().gen():
            P = self.parent()
            R = P.base_ring()
            v = [R.zero()]*exp + [R.one()]
            r = <SkewPolynomial_generic_dense>self._parent(v)
        else:
            r = <SkewPolynomial_generic_dense>self._new_c(list(self._coeffs), self._parent)
            r._inplace_pow(exp)

        if modulus:
            _, r = r.left_quo_rem(modulus)
        return r

    cpdef right_power_mod(self, exp, modulus):
        r"""
        Return the remainder of ``self**exp`` in the right euclidean division
        by ``modulus``.

        INPUT:

        - ``exp`` -- an Integer

        - ``modulus`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Remainder of ``self**exp`` in the right euclidean division
        by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: b = a^10  # short form for ``a._pow_(10)``
            sage: b == a*a*a*a*a*a*a*a*a*a
            True
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: br = a.right_power_mod(10,modulus); br
            (t^2 + t)*x^2 + (3*t^2 + 1)*x + t^2 + t
            sage: rq, rr = b.right_quo_rem(modulus)
            sage: br == rr
            True
            sage: a.right_power_mod(100,modulus)
            (2*t^2 + 3)*x^2 + (t^2 + 4*t + 2)*x + t^2 + 2*t + 1
        """
        cdef SkewPolynomial_generic_dense r
        if not isinstance(exp, Integer):
            try:
                exp = Integer(exp)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        if len(self._coeffs) <= 1:
            return self.parent()(self._coeffs[0]**exp)
        if exp == 0:
            return self.parent().one()
        if exp < 0:
            return (~self).right_power_mod(-exp, modulus)

        if self == self.parent().gen():
            P = self.parent()
            R = P.base_ring()
            v = [R.zero()]*exp + [R.one()]
            r = <SkewPolynomial_generic_dense>self._parent(v)
        else:
            r = <SkewPolynomial_generic_dense>self._new_c(list(self._coeffs), self._parent)
            r._inplace_pow(exp)

        if modulus:
            _, r = r.right_quo_rem(modulus)
        return r

    def __pow__(self, exp, modulus):
        r"""
        Return the remainder of ``self**exp`` in the left euclidean
        division by ``modulus``.

        INPUT:

        - ``exp`` -- an Integer

        - ``modulus`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Remainder of ``self**exp`` in the right euclidean division
        by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        .. SEEALSO::

            :meth:`~sage.rings.polynomial.skew_polynomial_element._pow_`

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: b = a^10
            sage: b == a*a*a*a*a*a*a*a*a*a
            True
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: bmod = a.right_power_mod(10,modulus); bmod
            (t^2 + t)*x^2 + (3*t^2 + 1)*x + t^2 + t
            sage: rq, rr = b.right_quo_rem(modulus)
            sage: bmod == rr
            True
        """
        return self.right_power_mod(exp, modulus)

    cpdef list coefficients(self, sparse=True):
        r"""
        Return the coefficients of the monomials appearing in ``self``.

        If ``sparse=True`` (the default), return only the non-zero coefficients.
        Otherwise, return the same value as ``self.list()``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = 1 + x^4 + (t+1)*x^2 + t^2
            sage: a.coefficients()
            [t^2 + 1, t + 1, 1]
            sage: a.coefficients(sparse=False)
            [t^2 + 1, 0, t + 1, 0, 1]
        """
        zero = self.parent().base_ring().zero()
        if sparse:
            return [c for c in self._coeffs if not c.is_zero()]
        else:
            return self._coeffs

cdef class ConstantSkewPolynomialSection(Map):
    r"""
    Representation of the canonical homomorphism from the constants of a skew
    polynomial ring to the base ring.

    This class is necessary for automatic coercion from zero-degree skew
    polynomial ring into the base ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.skew_polynomial_element import ConstantSkewPolynomialSection
        sage: R.<t> = QQ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: m = ConstantSkewPolynomialSection(S, R); m
        Generic map:
            From: Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
            To:   Univariate Polynomial Ring in t over Rational Field
    """
    cpdef Element _call_(self, x):
        r"""
        Return the corresponding element of the base ring if ``self`` is a
        constant skew polynomial. Otherwise, it fails.
        
        TESTS::
            
            sage: from sage.rings.polynomial.skew_polynomial_element import ConstantSkewPolynomialSection
            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: m = ConstantSkewPolynomialSection(S, R); m
            Generic map:
                From: Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
                To:   Univariate Polynomial Ring in t over Rational Field
            sage: m(S([0,1])-S([0,1]))
            0
            sage: m(S([3,1])-S([0,1]))
            3
            sage: m(S([0,1])-S([0,t]))
            Traceback (most recent call last):
            ...
            TypeError: not a constant polynomial
        """
        if x.degree() <= 0:
            try:
                return <Element>(x.constant_coefficient())
            except AttributeError:
                return <Element>((<SkewPolynomial>x).constant_coefficient())
        else:
            raise TypeError("not a constant polynomial")

cdef class SkewPolynomialBaseringInjection(Morphism):
    r"""
    Representation of the canonical homomorphism from a ring `R` into a skew
    polynomial ring over `R`.

    This class is necessary for automatic coercion from the base ring to the skew
    polynomial ring.

    .. SEEALSO::

        :class:~sage.rings.polynomial.polynomial_element.PolynomialBaseringInjection`

    EXAMPLES::

        sage: R.<t> = QQ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: S.coerce_map_from(S.base_ring()) #indirect doctest
        Skew Polynomial base injection morphism:
          From: Univariate Polynomial Ring in t over Rational Field
          To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
    """

    def __init__(self, domain, codomain):
        r"""
        Construct a Skew Polynomial Basering Injection.

        INPUT:

        - ``domain`` -- a ring `R`. This will be the domain of the injection.

        - ``codomain`` -- a skew polynomial ring over ``domain``. This will be
          the codomain.

        TESTS::

            sage: from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: SkewPolynomialBaseringInjection(k, k['x', Frob])
            Skew Polynomial base injection morphism:
              From: Finite Field in t of size 5^3
              To:   Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: R.<t> = QQ[]
            sage: SkewPolynomialBaseringInjection(QQ, k['x', Frob])
            Traceback (most recent call last):
            ...
            AssertionError: the domain of the injection must be the base ring of the skew polynomial ring
        """
        assert codomain.base_ring() is domain, \
            "the domain of the injection must be the base ring of the skew polynomial ring"
        Morphism.__init__(self, Hom(domain,codomain))
        self._an_element = codomain.gen()
        self._repr_type_str = "Skew Polynomial base injection"
        self._new_constant_poly_ = self._an_element._new_constant_poly

    def an_element(self):
        r"""
        Return an element of the codomain of the ring homomorphism.

        EXAMPLES::

            sage: from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: m = SkewPolynomialBaseringInjection(k, k['x', Frob])
            sage: m.an_element()
            x
        """
        return self._an_element

    cpdef Element _call_(self, e):
        r"""
        Return the corresponding skew polynomial to the element from the
        base ring according to ``self``.

        INPUT:

        - ``e`` -- element belonging to the base ring according to ``self``

        OUTPUT:

        The skew polynomial corresponding to `e` according to ``self``.

        TESTS::

            sage: from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: m = SkewPolynomialBaseringInjection(k, k['x', Frob])
            sage: m(4)
            4
            sage: parent(m(4))
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        try:
            return self._codomain._element_constructor_(e)
        except AttributeError:
            return self._codomain(e)

    def section(self):
        r"""
        Return the canonical homomorphism from the constants of a skew
        polynomial ring to the base ring according to ``self``.

        TESTS::

            sage: from sage.rings.polynomial.skew_polynomial_element import SkewPolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: m = SkewPolynomialBaseringInjection(k, k['x', Frob])
            sage: m.section()
            Generic map:
            From: Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            To:   Finite Field in t of size 5^3
        """
        return ConstantSkewPolynomialSection(self._codomain, self.domain())
