# sage.doctest: needs sage.combinat
r"""
Univariate Ore polynomials

This module provides the
:class:`~sage.rings.polynomial.skew_polynomial_element.OrePolynomial`,
which constructs a single univariate Ore polynomial over a commutative
base equipped with an endomorphism and/or a derivation.
It provides generic implementation of standard arithmetical operations
on Ore polynomials as addition, multiplication, gcd, lcm, etc.

The generic implementation of dense Ore polynomials is
:class:`~sage.rings.polynomial.ore_polynomial_element.OrePolynomial_generic_dense`.
The classes
:class:`~sage.rings.polynomial.ore_polynomial_element.ConstantOrePolynomialSection`
and :class:`~sage.rings.polynomial.ore_polynomial_element.OrePolynomialBaseringInjection`
handle conversion from an Ore polynomial ring to its base ring and vice versa.

AUTHORS:

- Xavier Caruso (2020-05)
"""

# ***************************************************************************
#    Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import re
from cysignals.signals cimport sig_check
from sage.structure.element import coerce_binop

from sage.rings.infinity import infinity
from sage.matrix.special import companion_matrix
from sage.structure.element cimport Element, RingElement, AlgebraElement
from sage.structure.parent cimport Parent
from sage.categories.homset import Hom
from sage.rings.ring import _Fields
from sage.rings.integer cimport Integer
from cpython.object cimport PyObject_RichCompare
from sage.categories.map cimport Map
from sage.rings.morphism cimport Morphism
from sage.rings.polynomial.polynomial_element cimport _dict_to_list


cdef class OrePolynomial(AlgebraElement):
    r"""
    Abstract base class for Ore polynomials.

    This class must be inherited from and have key methods overridden.

    .. RUBRIC:: Definition

    Let `R` be a commutative ring equipped with an automorphism `\sigma`
    and a `\sigma`-derivation `\partial`.

    An Ore polynomial is given by the equation:

    .. MATH::

        F(X) = a_{n} X^{n} + \cdots + a_0,

    where the coefficients `a_i \in R` and `X` is a formal variable.

    Addition between two Ore polynomials is defined by the usual addition
    operation and the modified multiplication is defined by the rule
    `X a = \sigma(a) X + \partial(a)` for all `a` in `R`.
    Ore polynomials are thus non-commutative and the degree of a product
    is equal to the sum of the degrees of the factors.

    Let `a` and `b` be two Ore polynomials in the same ring `S`.
    The *right (resp. left) Euclidean division* of `a` by `b` is a couple
    `(q,r)` of elements in `S` such that

    - `a = q b + r` (resp. `a = b q + r`)

    - the degree of `r` is less than the degree of `b`

    `q` (resp. `r`) is called the *quotient* (resp. the remainder)
    of this Euclidean division.

    .. RUBRIC:: Properties

    Keeping the previous notation, if the leading coefficient of `b`
    is a unit (e.g. if `b` is monic) then the quotient and the remainder
    in the *right* Euclidean division exist and are unique.

    The same result holds for the *left* Euclidean division if in addition
    the twisting morphism defining the Ore polynomial ring is invertible.

    EXAMPLES:

    We illustrate some functionalities implemented in this class.

    We create the Ore polynomial ring (here the derivation is zero)::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]; S
        Ore Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring
         twisted by t |--> t + 1

    and some elements in it::

        sage: a = t + x + 1; a
        x + t + 1
        sage: b = S([t^2,t+1,1]); b
        x^2 + (t + 1)*x + t^2
        sage: c = S.random_element(degree=3,monic=True)
        sage: c.parent() is S
        True

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

    Sage also implements arithmetic over Ore polynomial rings. You will find
    below a short panorama::

        sage: q,r = c.right_quo_rem(b)
        sage: c == q*b + r
        True

    The operators ``//`` and ``%`` give respectively the quotient
    and the remainder of the *right* Euclidean division::

        sage: q == c // b
        True
        sage: r == c % b
        True

    Here we can see the effect of the operator evaluation compared to the usual
    polynomial evaluation::

        sage: a = x^2
        sage: a(t)
        doctest:...: FutureWarning: This class/method/function is marked as experimental.
        It, its functionality or its interface might change without a formal deprecation.
        See https://github.com/sagemath/sage/issues/13215 for details.
        t + 2

    Here is another example over a finite field::

        sage: # needs sage.rings.finite_rings
        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: a = x^4 + (4*t + 1)*x^3 + (t^2 + 3*t + 3)*x^2 + (3*t^2 + 2*t + 2)*x + (3*t^2 + 3*t + 1)
        sage: b = (2*t^2 + 3)*x^2 + (3*t^2 + 1)*x + 4*t + 2
        sage: q, r = a.left_quo_rem(b)
        sage: q
        (4*t^2 + t + 1)*x^2 + (2*t^2 + 2*t + 2)*x + 2*t^2 + 4*t + 3
        sage: r
        (t + 2)*x + 3*t^2 + 2*t + 4
        sage: a == b*q + r
        True

    Once we have Euclidean divisions, we have for free gcd and lcm
    (at least if the base ring is a field)::

        sage: # needs sage.rings.finite_rings
        sage: a = (x + t) * (x + t^2)^2
        sage: b = (x + t) * (t*x + t + 1) * (x + t^2)
        sage: a.right_gcd(b)
        x + t^2
        sage: a.left_gcd(b)
        x + t

    The left lcm has the following meaning: given Ore polynomials `a` and `b`,
    their left lcm is the least degree polynomial `c = ua = vb` for some Ore
    polynomials `u, v`. Such a `c` always exist if the base ring is a field::

        sage: c = a.left_lcm(b); c                                                      # needs sage.rings.finite_rings
        x^5 + (4*t^2 + t + 3)*x^4 + (3*t^2 + 4*t)*x^3 + 2*t^2*x^2 + (2*t^2 + t)*x + 4*t^2 + 4
        sage: c.is_right_divisible_by(a)                                                # needs sage.rings.finite_rings
        True
        sage: c.is_right_divisible_by(b)                                                # needs sage.rings.finite_rings
        True

    The right lcm is defined similarly as the least degree polynomial `c = au =
    bv` for some `u,v`::

        sage: d = a.right_lcm(b); d                                                     # needs sage.rings.finite_rings
        x^5 + (t^2 + 1)*x^4 + (3*t^2 + 3*t + 3)*x^3 + (3*t^2 + t + 2)*x^2 + (4*t^2 + 3*t)*x + 4*t + 4
        sage: d.is_left_divisible_by(a)                                                 # needs sage.rings.finite_rings
        True
        sage: d.is_left_divisible_by(b)                                                 # needs sage.rings.finite_rings
        True

    .. SEEALSO::

        - :mod:`sage.rings.polynomial.ore_polynomial_ring`
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

    cdef long _hash_c(self) noexcept:
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

    cpdef Integer degree(self):
        r"""
        Return the degree of ``self``.

        By convention, the zero Ore polynomial has degree `-1`.

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
        raise NotImplementedError

    cdef OrePolynomial _new_c(self, list coeffs, Parent P, char check=0):
        r"""
        Fast creation of a new Ore polynomial.

        .. NOTE::

            Override this function in classes which inherit
            from OrePolynomial.
        """
        return P(coeffs)

    cpdef OrePolynomial _new_constant_poly(self, RingElement a, Parent P, char check=0):
        r"""
        Fast creation of a new constant Ore polynomial.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.polynomial.ore_polynomial_element import OrePolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: OrePolynomialBaseringInjection(k, k['x', Frob]) #indirect doctest
            Ore Polynomial base injection morphism:
              From: Finite Field in t of size 5^3
              To:   Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        if a:
            return self._new_c([a], P, check)
        return self._new_c([], P)

    def __setitem__(self, n, value):
        r"""
        Set the ``n``-th coefficient of ``self``.

        This always raises an :exc:`IndexError`, since polynomials are immutable in
        Sage.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: a[1] = t + 1
            Traceback (most recent call last):
            ...
            IndexError: Ore polynomials are immutable
        """
        raise IndexError("Ore polynomials are immutable")

    def square(self):
        r"""
        Return the square of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x', sigma]
            sage: a = x + t; a
            x + t
            sage: a.square()
            x^2 + (2*t + 1)*x + t^2
            sage: a.square() == a*a
            True

            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: (d + t).square()
            d^2 + 2*t*d + t^2 + 1
        """
        return self * self

    def constant_coefficient(self):
        r"""
        Return the constant coefficient (i.e., the coefficient of term
        of degree `0`) of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t^2 + 2
            sage: a.constant_coefficient()
            t^2 + 2
        """
        if self:
            return self[0]
        else:
            return self.base_ring().zero()

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

        By convention, the leading coefficient to the zero polynomial is
        zero::

            sage: S(0).leading_coefficient()
            0
        """
        cdef int d = self.degree()
        if d == -1:
            return self.base_ring()(0)
        return self[d]

    def is_unit(self):
        r"""
        Return ``True`` if this Ore polynomial is a unit.

        When the base ring `R` is an integral domain, then an Ore polynomial `f`
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
            raise NotImplementedError("is_unit is not implemented for Ore polynomial rings "
                                      "over base rings which are not integral domains")

    def is_nilpotent(self):
        r"""
        Check if ``self`` is nilpotent.

        .. NOTE::

            The paper "Nilpotents and units in skew polynomial rings
            over commutative rings" by M. Rimmer and K.R. Pearson describes
            a method to check whether a given skew polynomial is nilpotent.
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
        Return ``True`` if this Ore polynomial is monic.

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
        return self.leading_coefficient() == 1

    def left_monic(self):
        r"""
        Return the unique monic Ore polynomial `m` which divides this
        polynomial on the left and has the same degree.

        Given an Ore polynomial `P` of degree `n`, its left monic is given by
        `P \cdot \sigma^{-n}(1/k)`, where `k` is the leading coefficient of
        `P` and `\sigma` is the twisting morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.left_monic(); b
            x^3 + (4*t^2 + 3*t)*x^2 + (4*t + 2)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: # needs sage.rings.finite_rings
            sage: b.degree() == a.degree()
            True
            sage: a.is_left_divisible_by(b)
            True
            sage: twist = S.twisting_morphism(-a.degree())
            sage: a == b * twist(a.leading_coefficient())
            True

        Note that `b` does not divide `a` on the right::

            sage: a.is_right_divisible_by(b)                                            # needs sage.rings.finite_rings
            False

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: S.<x> = R['x', der]
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
        cdef Morphism m = self._parent.twisting_morphism(-self.degree())
        if m is None:
            return self * a
        return self * m(a)

    def right_monic(self):
        r"""
        Return the unique monic Ore polynomial which divides this polynomial
        on the right and has the same degree.

        Given an Ore polynomial `P` of degree `n`, its right monic is given by
        `(1/k) \cdot P`, where `k` is the leading coefficient of `P`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = a.right_monic(); b
            x^3 + (2*t^2 + 3*t + 4)*x^2 + (3*t^2 + 4*t + 1)*x + 2*t^2 + 4*t + 3

        Check list::

            sage: b.degree() == a.degree()                                              # needs sage.rings.finite_rings
            True
            sage: a.is_right_divisible_by(b)                                            # needs sage.rings.finite_rings
            True
            sage: a == a.leading_coefficient() * b                                      # needs sage.rings.finite_rings
            True

        Note that `b` does not divide `a` on the right::

            sage: a.is_left_divisible_by(b)                                             # needs sage.rings.finite_rings
            False

        This function does not work if the leading coefficient is not a
        unit::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: S.<x> = R['x', der]
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
        return a * self

    cpdef _mod_(self, other):
        r"""
        Return the remainder in the *right* Euclidean division of
        ``self`` by ``other``.

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
        return self.right_quo_rem(other)[1]

    cpdef _floordiv_(self, right):
        r"""
        Return the quotient of the *right* Euclidean division of
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
        q, _ = self.right_quo_rem(right)
        return q

    cpdef _div_(self, right):
        r"""
        Return the quotient of this Ore polynomial by ``right``
        in the fraction field.

        INPUT:

        - ``right`` -- an Ore polynomial

        EXAMPLES::

            sage: R.<t> = GF(11)[]
            sage: der = R.derivation()
            sage: S.<x> = R['x', der]
            sage: f = t/x
            sage: f
            (x + 10/t)^(-1) * t
            sage: f.parent()
            Ore Function Field in x over
             Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 11 twisted by d/dt
        """
        parent = self.parent().fraction_field()
        return parent(self) / parent(right)

    def is_left_divisible_by(self, other):
        r"""
        Check if ``self`` is divisible by ``other`` on the left.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: c.is_left_divisible_by(S(0))                                          # needs sage.rings.finite_rings
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

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: c.is_right_divisible_by(S(0))                                         # needs sage.rings.finite_rings
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

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a * b
            sage: a.left_divides(c)
            True
            sage: b.left_divides(c)
            False

        Divisibility by `0` does not make sense::

            sage: S(0).left_divides(c)                                                  # needs sage.rings.finite_rings
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

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^2 + t*x + t^2 + 3
            sage: b = x^3 + (t + 1)*x^2 + 1
            sage: c = a * b
            sage: a.right_divides(c)
            False
            sage: b.right_divides(c)
            True

        Divisibility by `0` does not make sense::

            sage: S(0).right_divides(c)                                                 # needs sage.rings.finite_rings
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

        If `a` is ``self`` and `b` is ``other``, then there are Ore polynomials
        `u` and `v` such that `g = a u + b v`, where `g` is the left gcd of `a`
        and `b`. This method returns `(g, u, v)`.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); return whether the left gcd
          should be normalized to be monic

        OUTPUT:

        - The left gcd of ``self`` and ``other``, that is an Ore polynomial
          `g` with the following property: any Ore polynomial is
          divisible on the left by `g` iff it is divisible on the left
          by both ``self`` and ``other``.
          If monic is ``True``, `g` is in addition monic. (With this
          extra condition, it is uniquely determined.)

        - Two Ore polynomials `u` and `v` such that:

          .. MATH::

              g = a \cdot u + b \cdot v,

          where `s` is ``self`` and `b` is ``other``.

        .. NOTE::

            Works only if following two conditions are fulfilled
            (otherwise left gcd do not exist in general):
            1) the base ring is a field and
            2) the twisting morphism is bijective.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: g,u,v = a.left_xgcd(b, monic=False); g                                # needs sage.rings.finite_rings
            2*t*x + 4*t + 2
            sage: a*u + b*v == g                                                        # needs sage.rings.finite_rings
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

        And the twisting morphism must be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_xgcd(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twisting morphism Ring endomorphism of Fraction Field of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t^2
        """
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        cdef OrePolynomial G = self
        cdef OrePolynomial U = self._parent.one()
        cdef OrePolynomial V, V1, V2, Q, R, T
        cdef Morphism m
        if not other:
            V = self._parent.zero()
        else:
            V1 = self._parent.zero()
            V3 = other
            while V3:
                Q, R = G._left_quo_rem(V3)
                T = U - V1 * Q
                U = V1
                G = V3
                V1 = T
                V3 = R
            V = G - self * U
            V, _ = V._left_quo_rem(other)
        if monic:
            lc = ~G.leading_coefficient()
            m = self._parent.twisting_morphism(-G.degree())
            if m is not None:
                lc = m(lc)
            G = G * lc
            U = U * lc
            V = V * lc
        return G, U, V

    cdef _left_quo_rem(self, OrePolynomial other):
        r"""
        Return the quotient and remainder of the left Euclidean
        division of ``self`` by ``other`` (C implementation).

        Must be implemented in subclasses.
        """
        raise NotImplementedError("left Euclidean division is not implemented")

    @coerce_binop
    def left_quo_rem(self, other):
        r"""
        Return the quotient and remainder of the left Euclidean
        division of ``self`` by ``other``.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        OUTPUT:

        - the quotient and the remainder of the left Euclidean
          division of this Ore polynomial by ``other``

        .. NOTE::

            This will fail if the leading coefficient of ``other`` is not a unit
            or if Sage can't invert the twisting morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (3*t^2 + 3*t + 2)*x^3 + (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: b = (3*t^2 + 4*t + 2)*x^2 + (2*t^2 + 4*t + 3)*x + 2*t^2 + t + 1
            sage: q,r = a.left_quo_rem(b)
            sage: a == b*q + r
            True

        In the following example, Sage does not know the inverse
        of the twisting morphism::

            sage: R.<t> = QQ[]
            sage: K = R.fraction_field()
            sage: sigma = K.hom([(t+1)/(t-1)])
            sage: S.<x> = K['x',sigma]
            sage: a = (-2*t^2 - t + 1)*x^3 + (-t^2 + t)*x^2 + (-12*t - 2)*x - t^2 - 95*t + 1
            sage: b = x^2 + (5*t - 6)*x - 4*t^2 + 4*t - 1
            sage: a.left_quo_rem(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twisting morphism Ring endomorphism of Fraction Field of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> (t + 1)/(t - 1)
        """
        if not other:
            raise ZeroDivisionError("division by zero is not valid")
        return self._left_quo_rem(other)

    cdef _right_quo_rem(self, OrePolynomial other):
        r"""
        Return the quotient and remainder of the right Euclidean
        division of ``self`` by ``other`` (C implementation).

        Must be implemented in subclasses.
        """
        raise NotImplementedError("right Euclidean division is not implemented")

    @coerce_binop
    def right_quo_rem(self, other):
        r"""
        Return the quotient and remainder of the right Euclidean
        division of ``self`` by ``other``.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        OUTPUT:

        - the quotient and the remainder of the right Euclidean
          division of this Ore polynomial by ``other``

        .. NOTE::

            This will fail if the leading coefficient of the divisor
            is not a unit.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(degree=4)
            sage: b = S.random_element(monic=True)
            sage: q,r = a.right_quo_rem(b)
            sage: a == q*b + r
            True

        The leading coefficient of the divisor need to be invertible::

            sage: a.right_quo_rem(S(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero is not valid
            sage: c = S.random_element()
            sage: while not c or c.leading_coefficient().is_unit():
            ....:     c = S.random_element()
            sage: while a.degree() < c.degree():
            ....:     a = S.random_element(degree=4)
            sage: a.right_quo_rem(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: the leading coefficient of the divisor is not invertible
        """
        if not other:
            raise ZeroDivisionError("division by zero is not valid")
        return self._right_quo_rem(other)

    @coerce_binop
    def right_xgcd(self, other, monic=True):
        r"""
        Return the right gcd of ``self`` and ``other`` along with the
        coefficients for the linear combination.

        If `a` is ``self`` and `b` is ``other``, then there are Ore polynomials
        `u` and `v` such that `g = u a + v b`, where `g` is the right gcd of `a`
        and `b`. This method returns `(g, u, v)`.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); return whether the right gcd
          should be normalized to be monic

        OUTPUT:

        - The right gcd of ``self`` and ``other``, that is an Ore polynomial
          `g` with the following property: any Ore polynomial is
          divisible on the right by `g` iff it is divisible on the right
          by both ``self`` and ``other``.
          If monic is ``True``, `g` is in addition monic. (With this
          extra condition, it is uniquely determined.)

        - Two Ore polynomials `u` and `v` such that:

          .. MATH::

              g = u \cdot a + v \cdot b

          where `a` is ``self`` and `b` is ``other``.

        .. NOTE::

            Works only if the base ring is a field (otherwise right
            gcd do not exist in general).

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: g,u,v = a.right_xgcd(b, monic=False); g                               # needs sage.rings.finite_rings
            (4*t^2 + 4*t + 1)*x + 4*t^2 + 4*t + 3
            sage: u*a + v*b == g                                                        # needs sage.rings.finite_rings
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
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        cdef OrePolynomial G = self
        cdef OrePolynomial U = self._parent.one()
        cdef OrePolynomial V, V1, V3, Q, R, T
        if other.is_zero():
            V = self._parent.zero()
        else:
            V1 = self._parent.zero()
            V3 = other
            while not V3.is_zero():
                Q, R = G._right_quo_rem(V3)
                T = U - Q * V1
                U = V1
                G = V3
                V1 = T
                V3 = R
            V = G - U * self
            V, _ = V._right_quo_rem(other)
        if monic:
            lc = ~G.leading_coefficient()
            G = lc * G
            U = lc * U
            V = lc * V
        return G, U, V

    @coerce_binop
    def right_gcd(self, other, monic=True):
        r"""
        Return the right gcd of ``self`` and ``other``.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); return whether the right gcd
          should be normalized to be monic

        OUTPUT:

        The right gcd of ``self`` and ``other``, that is an Ore polynomial
        `g` with the following property: any Ore polynomial is
        divisible on the right by `g` iff it is divisible on the right
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field (otherwise right
            gcd do not exist in general).

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x^2 + t*x + 1) * (x + t)
            sage: b = 2 * (x^3 + (t+1)*x^2 + t^2) * (x + t)
            sage: a.right_gcd(b)
            x + t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.right_gcd(b,monic=False)                                            # needs sage.rings.finite_rings
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
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        if other.is_zero():
            return self
        cdef OrePolynomial A = self
        cdef OrePolynomial B = other
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

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); return whether the left gcd
          should be normalized to be monic

        OUTPUT:

        The left gcd of ``self`` and ``other``, that is an Ore polynomial
        `g` with the following property: any Ore polynomial is
        divisible on the left by `g` iff it is divisible on the left
        by both ``self`` and ``other``.
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if following two conditions are fulfilled
            (otherwise left gcd do not exist in general):
            1) the base ring is a field and
            2) the twisting morphism is bijective.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_gcd(b)
            x + t

        Specifying ``monic=False``, we *can* get a nonmonic gcd::

            sage: a.left_gcd(b,monic=False)                                             # needs sage.rings.finite_rings
            2*t*x + 4*t + 2

        The base ring needs to be a field::

            sage: # needs sage.rings.finite_rings
            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t + 1])
            sage: S.<x> = R['x',sigma]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_gcd(b)
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be a field

        And the twisting morphism needs to be bijective::

            sage: # needs sage.rings.finite_rings
            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x^2 + t*x + 1)
            sage: b = 2 * (x + t) * (x^3 + (t+1)*x^2 + t^2)
            sage: a.left_gcd(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twisting morphism Ring endomorphism
            of Fraction Field of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t^2
        """
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        if other.is_zero():
            return self
        cdef OrePolynomial A_, A = self
        cdef OrePolynomial B = other
        while not B.is_zero():
            A_ = A
            A = B
            _, B = A_._left_quo_rem(B)
        if monic:
            A = A.left_monic()
        return A

    cdef OrePolynomial _left_lcm_cofactor(self, OrePolynomial other):
        r"""
        Return an Ore polynomial `U` such that `U P = c L`
        where `P` is this Ore polynomial (``self``), `L`
        is the left lcm of `P` and ``other`` and `c` is a
        constant

        TESTS::

            sage: # needs sage.misc.cython sage.rings.finite_rings
            sage: cython(
            ....: '''
            ....: from sage.rings.polynomial.ore_polynomial_element cimport OrePolynomial
            ....: def left_lcm_cofactor(OrePolynomial P, OrePolynomial Q):
            ....:     return P._left_lcm_cofactor(Q)
            ....: ''')
            sage: k.<a> = GF(7^5)
            sage: Frob = k.frobenius_endomorphism(3)
            sage: S.<x> = k['x', Frob]
            sage: D = S.random_element(degree=2)
            sage: P = S.random_element(degree=2) * D
            sage: Q = S.random_element(degree=2) * D
            sage: L = P.left_lcm(Q)
            sage: U = left_lcm_cofactor(P, Q)
            sage: (U*P).right_monic() == L
            True
        """
        cdef OrePolynomial Q, R, T
        cdef OrePolynomial U = <OrePolynomial>self._parent.one()
        cdef OrePolynomial V = <OrePolynomial>self._parent.zero()
        while other:
            Q, R = self._right_quo_rem(other)
            T = U - Q * V
            U = V
            V = T
            self = other
            other = R
        return V

    @coerce_binop
    def left_xlcm(self, other, monic=True):
        r"""
        Return the left lcm `L` of ``self`` and ``other`` together
        with two Ore polynomials `U` and `V` such that

        .. MATH::

            U \cdot \text{self} = V \cdot \text{other} = L.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: P = (x + t^2) * (x + t)
            sage: Q = 2 * (x^2 + t + 1) * (x * t)
            sage: L, U, V = P.left_xlcm(Q)
            sage: L
            x^5 + (2*t^2 + t + 4)*x^4 + (3*t^2 + 4)*x^3 + (3*t^2 + 3*t + 2)*x^2 + (t^2 + t + 2)*x

            sage: U * P == L                                                            # needs sage.rings.finite_rings
            True
            sage: V * Q == L                                                            # needs sage.rings.finite_rings
            True
        """
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError("division by zero is not valid")
        cdef OrePolynomial V1 = self._left_lcm_cofactor(other)
        cdef OrePolynomial L = V1 * self
        if monic:
            s = ~(L.leading_coefficient())
            L = s * L
            V1 = s * V1
        return L, V1, L // other

    cdef OrePolynomial _right_lcm_cofactor(self, OrePolynomial other):
        r"""
        Return an Ore polynomial `U` such that `P U = L c`
        where `P` is this Ore polynomial (``self``), `L`
        is the right lcm of `P` and ``other`` and `c` is a
        constant

        TESTS::

            sage: # needs sage.misc.cython sage.rings.finite_rings
            sage: cython(
            ....: '''
            ....: from sage.rings.polynomial.ore_polynomial_element cimport OrePolynomial
            ....: def right_lcm_cofactor(OrePolynomial P, OrePolynomial Q):
            ....:     return P._right_lcm_cofactor(Q)
            ....: ''')
            sage: k.<a> = GF(7^5)
            sage: Frob = k.frobenius_endomorphism(3)
            sage: S.<x> = k['x', Frob]
            sage: D = S.random_element(degree=2)
            sage: P = D * S.random_element(degree=2)
            sage: Q = D * S.random_element(degree=2)
            sage: L = P.right_lcm(Q)
            sage: U = right_lcm_cofactor(P, Q)
            sage: (P*U).left_monic() == L
            True
        """
        cdef OrePolynomial Q, R, T
        cdef OrePolynomial U = <OrePolynomial>self._parent.one()
        cdef OrePolynomial V = <OrePolynomial>self._parent.zero()
        while other:
            Q, R = self._left_quo_rem(other)
            T = U - V * Q
            U = V
            V = T
            self = other
            other = R
        return V

    @coerce_binop
    def right_xlcm(self, other, monic=True):
        r"""
        Return the right lcm `L` of ``self`` and ``other`` together
        with two Ore polynomials `U` and `V` such that

        .. MATH::

            \text{self} \cdot U = \text{other} \cdot V = L.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); whether the right lcm
          should be normalized to be monic

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: P = (x + t) * (x + t^2)
            sage: Q = 2 * (x + t) * (x^2 + t + 1)
            sage: L, U, V = P.right_xlcm(Q)
            sage: L
            x^4 + (2*t^2 + t + 2)*x^3 + (3*t^2 + 4*t + 1)*x^2 + (3*t^2 + 4*t + 1)*x + t^2 + 4
            sage: P * U == L
            True
            sage: Q * V == L
            True
        """
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError("division by zero is not valid")
        cdef OrePolynomial V1 = self._right_lcm_cofactor(other)
        cdef OrePolynomial L = self * V1
        if monic:
            s = self.base_ring()(~L.leading_coefficient())
            s = self._parent.twisting_morphism(-L.degree())(s)
            L = L * s
            V1 = V1 * s
        W1, _ = L._left_quo_rem(other)
        return L, V1, W1

    @coerce_binop
    def left_lcm(self, other, monic=True):
        r"""
        Return the left lcm of ``self`` and ``other``.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); return whether the left lcm
          should be normalized to be monic

        OUTPUT:

        The left lcm of ``self`` and ``other``, that is an Ore polynomial
        `l` with the following property: any Ore polynomial is a left multiple of `l` (right divisible by `l`)
        iff it is a left multiple of both ``self`` and ``other`` (right divisible by ``self`` and ``other``).
        If monic is ``True``, `l` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if the base ring is a field (otherwise left
            lcm do not exist in general).

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

        Specifying ``monic=False``, we *can* get a nonmonic lcm::

            sage: a.left_lcm(b,monic=False)                                             # needs sage.rings.finite_rings
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
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError("division by zero is not valid")
        L = self._left_lcm_cofactor(other) * self
        if monic:
            L = L.right_monic()
        return L

    @coerce_binop
    def right_lcm(self, other, monic=True):
        r"""
        Return the right lcm of ``self`` and ``other``.

        INPUT:

        - ``other`` -- an Ore polynomial in the same ring as ``self``

        - ``monic`` -- boolean (default: ``True``); return whether the right lcm
          should be normalized to be monic

        OUTPUT:

        The right lcm of ``self`` and ``other``, that is an Ore polynomial
        `l` with the following property: any Ore polynomial is a right multiple of `g` (left divisible by `l`)
        iff it is a right multiple of both ``self`` and ``other`` (left divisible by ``self`` and ``other``).
        If monic is ``True``, `g` is in addition monic. (With this
        extra condition, it is uniquely determined.)

        .. NOTE::

            Works only if two following conditions are fulfilled
            (otherwise right lcm do not exist in general):
            1) the base ring is a field and 2) the twisting morphism on
            this field is bijective.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
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

            sage: a.right_lcm(b,monic=False)                                            # needs sage.rings.finite_rings
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

        And the twisting morphism needs to be bijective::

            sage: FR = R.fraction_field()
            sage: f = FR.hom([FR(t)^2])
            sage: S.<x> = FR['x',f]
            sage: a = (x + t) * (x + t^2)
            sage: b = 2 * (x + t) * (x^2 + t + 1)
            sage: a.right_lcm(b)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twisting morphism Ring endomorphism of
            Fraction Field of Univariate Polynomial Ring in t over Rational Field
                Defn: t |--> t^2
        """
        if self.base_ring() not in _Fields:
            raise TypeError("the base ring must be a field")
        if self.is_zero() or other.is_zero():
            raise ZeroDivisionError("division by zero is not valid")
        L = self * self._right_lcm_cofactor(other)
        if monic:
            L = L.left_monic()
        return L

    def _repr_(self, name=None):
        r"""
        Return string representation of this Ore polynomial.

        INPUT:

        - ``name`` -- the name of the variable (default: the
          name given when the Ore polynomial ring was created)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t
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
                if n < m - 1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)" % x
                if n > 1:
                    var = "*%s^%s" % (name, n)
                elif n == 1:
                    var = "*%s" % name
                else:
                    var = ""
                s += "%s%s" % (x, var)
        s = s.replace(" + -", " - ")
        s = re.sub(r' 1(\.0+)?\*', ' ', s)
        s = re.sub(r' -1(\.0+)?\*', ' -', s)
        if s == " ":
            return "0"
        return s[1:]

    def _latex_(self, name=None):
        r"""
        Return a latex representation of this Ore polynomial.

        INPUT:

        - ``name`` -- the name of the variable (default: the
          name given when the Ore polynomial ring was created)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t^2 + 1/2*x*t
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
                if n != m - 1:
                    s += " + "
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "\\left(%s\\right)" % x
                if n > 1:
                    var = "|%s^{%s}" % (name, n)
                elif n == 1:
                    var = "|%s" % name
                else:
                    var = ""
                s += "%s %s" % (x, var)
        s = s.replace(" + -", " - ")
        s = re.sub(r" 1(\.0+)? \|", " ", s)
        s = re.sub(r" -1(\.0+)? \|", " -", s)
        s = s.replace("|", "")
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

    def __bool__(self):
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
        return self.degree() > -1

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
            return self._parent(n * [self.base_ring().zero()] + self.list(), check=False)
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

        This will create the Ore polynomial ring with the new name but same
        base ring, twisting morphism and twisting derivation. The returned
        Ore polynomial will be an element of that Ore polynomial ring.

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
            Ore Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring
             twisted by t |--> t + 1
        """
        R = self._parent.change_var(var)
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

        If ``sparse=True`` (the default), return only the nonzero coefficients.
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
        Return the number of nonzero coefficients of ``self``.

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

        In Sage, since Ore polynomials are immutable, this just returns
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

    cpdef bint is_zero(self) noexcept:
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

    cpdef bint is_one(self) noexcept:
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

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = 1 + t*x^2
            sage: b = x + 1
            sage: a.left_mod(b)
            2*t^2 + 4*t
        """
        _, r = self.left_quo_rem(other)
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
        return [i for i in range(self.degree() + 1) if bool(self[i])]

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

        - ``n`` -- (default: ``None``) if given, an integer that
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
            return v + [z] * (n - len(v))
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


cdef void lmul_gen(list A, Morphism m, d) noexcept:
    r"""
    If ``A`` is the list of coefficients of an Ore polynomial ``P``,
    replace it by the list of coefficients of ``X*P`` (where ``X``
    is the variable in the Ore polynomial ring).

    This is a helper function.

    INPUT:

    - ``A`` -- list of coefficients

    - ``m`` -- the twisting morphism of the Ore polynomial ring

    - ``d`` -- the twisting derivation of the Ore polynomial ring
    """
    if m is None:
        A.append(A[-1])
        for j in range(len(A) - 2, 0, -1):
            A[j] = d(A[j]) + A[j - 1]
    else:
        A.append(m(A[-1]))
        for j in range(len(A) - 2, 0, -1):
            A[j] = d(A[j]) + m(A[j - 1])
    A[0] = d(A[0])


cdef class OrePolynomial_generic_dense(OrePolynomial):
    r"""
    Generic implementation of dense Ore polynomial supporting any valid base
    ring, twisting morphism and twisting derivation.
    """
    def __init__(self, parent, x=None, int check=1, int construct=0, **kwds):
        r"""
        Construct an Ore polynomial over the given parent with the given
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

        We create an Ore polynomial from a list::

            sage: S([t,1])
            x + t

        from another Ore polynomial::

            sage: S(x^2 + t)
            x^2 + t

        from a constant::

            sage: x = S(t^2 + 1); x
            t^2 + 1
            sage: x.parent() is S
            True
        """
        OrePolynomial.__init__(self, parent)
        if x is None:
            self._coeffs = []
            return

        R = parent.base_ring()
        if isinstance(x, list):
            if check:
                self._coeffs = [R(t) for t in x]
                self._normalize()
            else:
                self._coeffs = x
            return

        if isinstance(x, OrePolynomial):
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
                    self._normalize()
                return

        elif isinstance(x, int) and x == 0:
            self._coeffs = []
            return

        elif isinstance(x, dict):
            x = _dict_to_list(x, R.zero())

        elif not isinstance(x, list):
            x = [x]
        if check:
            self._coeffs = [R(z, **kwds) for z in x]
            self._normalize()
        else:
            self._coeffs = x

    def __reduce__(self):
        r"""
        Return the generic dense Ore polynomial corresponding to the
        current parameters provided ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: loads(dumps(x)) == x
            True
            sage: loads(dumps(x))
            x
        """
        return (self._parent, (self._coeffs,))

    cdef long _hash_c(self) noexcept:
        r"""
        This hash incorporates the name of the variable.

        .. NOTE::

            This is an internal method. Use :meth:`__hash__` instead.
        """
        # todo - come up with a way to create hashes of zero that
        # that do not incorrectly indicate that the element is 0.
        cdef long result = 0
        cdef long result_mon
        cdef long c_hash
        cdef long var_name_hash = 0
        cdef int i
        for i in range(len(self._coeffs)):
            if i == 1:
                var_name_hash = hash(self._parent._names[0])
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

    cdef OrePolynomial _new_c(self, list coeffs, Parent P, char check=0):
        r"""
        Fast creation of a new Ore polynomial given a list of coefficients.

        .. WARNING::

            The list ``coeffs`` is stored internally in the newly created Ore
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
        cdef OrePolynomial_generic_dense f = t.__new__(t)
        f._parent = P
        f._coeffs = coeffs
        if check:
            f._normalize()
        return f

    cdef void _normalize(self) noexcept:
        r"""
        Remove higher order `0`-coefficients from the representation of ``self``.

        TESTS::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S #indirect doctest
            Ore Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
        """
        cdef list x = self._coeffs
        cdef Py_ssize_t n = len(x) - 1
        while n >= 0 and not x[n]:
            del x[n]
            n -= 1

    cpdef _richcmp_(left, right, int op):
        r"""
        Compare the two Ore polynomials ``self`` and ``other``.

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
        cdef x = (<OrePolynomial_generic_dense>left)._coeffs
        cdef y = (<OrePolynomial_generic_dense>right)._coeffs
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
        return iter((<OrePolynomial_generic_dense>self)._coeffs)

    def __getitem__(self, n):
        r"""
        Return the `n`-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer

        OUTPUT: the ``n``-th coefficient of ``self``

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
            return self._coeffs[n]
        except IndexError:
            return self.base_ring().zero()

    cpdef list list(self, bint copy=True):
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
        if copy:
            # This creates a shallow copy
            return list((<OrePolynomial_generic_dense>self)._coeffs)
        else:
            return (<OrePolynomial_generic_dense>self)._coeffs

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
        cdef list Y = (<OrePolynomial_generic_dense>self)._coeffs
        cdef int i
        for i in range(len(Y)):
            c = Y[i]
            if c:
                X[i] = c
        return X

    cpdef Integer degree(self):
        r"""
        Return the degree of ``self``.

        By convention, the zero Ore polynomial has degree `-1`.

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

        TESTS:

        We check that the degree is an ``Integer`` object (see #35519)::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t*x^3 + t^2*x + 1
            sage: isinstance(a.degree(), Integer)
            True

        ::

            sage: R.<t> = OrePolynomialRing(GF(5)['T'], GF(5)['T'].frobenius_endomorphism())
            sage: isinstance((t + 1).degree(), Integer)
            True
        """
        return Integer(len(self._coeffs) - 1)

    cpdef _add_(self, right):
        r"""
        Add two polynomials.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = S.random_element(monic=True)
            sage: b = -S.random_element(monic=True)
            sage: c = a + b
            sage: c - b == a
            True
            sage: c - b == a
            True
            sage: a + a == 2*a
            True
        """
        cdef Py_ssize_t i, min
        cdef list x = (<OrePolynomial_generic_dense>self)._coeffs
        cdef list y = (<OrePolynomial_generic_dense>right)._coeffs
        cdef Py_ssize_t dx = len(x), dy = len(y)
        if dx > dy:
            r = self._new_c([x[i] + y[i] for i in range(dy)] + x[dy:], self._parent, 0)
        elif dx < dy:
            r = self._new_c([x[i] + y[i] for i in range(dx)] + y[dx:], self._parent, 0)
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
            sage: a = S.random_element(monic=True)
            sage: b = S.random_element(monic=True)
            sage: c = a - b
            sage: c - a == -b
            True
            sage: c + b == a
            True
            sage: c - c == 0
            True
        """
        cdef Py_ssize_t i, min
        cdef list x = (<OrePolynomial_generic_dense>self)._coeffs
        cdef list y = (<OrePolynomial_generic_dense>right)._coeffs
        cdef Py_ssize_t dx = len(x), dy = len(y)
        cdef RingElement c
        if dx > dy:
            r = self._new_c([x[i] - y[i] for i in range(dy)] + x[dy:], self._parent, 0)
        elif dx < dy:
            r = self._new_c([x[i] - y[i] for i in range(dx)] + [-c for c in y[dx:]], self._parent, 0)
        else:
            r = self._new_c([x[i] - y[i] for i in range(dx)], self._parent, 1)
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
        c = self._new_c([-x for x in (<OrePolynomial_generic_dense>self)._coeffs],
                        self._parent, 0)
        return c

    def valuation(self):
        r"""
        Return the minimal degree of a nonzero monomial of ``self``.

        By convention, the zero Ore polynomial has valuation `+\infty`.

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

    cdef list _mul_list(self, list A):
        r"""
        Return the list of coefficients of the product of this
        Ore polynomial by that whose coefficients are given by ``A``.

        This is a helper function.
        """
        cdef list BA = [self.base_ring().zero()] * (len(self._coeffs) + len(A) - 1)
        cdef Morphism m = self._parent._morphism
        cdef d = self._parent._derivation
        cdef coeff
        cdef Py_ssize_t i, j
        for i in range(len(self._coeffs)):
            if i:
                lmul_gen(A, m, d)
            coeff = self._coeffs[i]
            if coeff:
                for j in range(len(A)):
                    BA[j] += coeff * A[j]
        return BA

    cpdef _lmul_(self, Element s):
        r"""
        Return the product ``self * right``.

        INPUT:

        - ``right`` -- an element of the base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: d*t  # indirect doctest
            t*d + 1
        """
        cdef coeffs = self._mul_list([s])
        return self._new_c(coeffs, self._parent, 1)

    cpdef _rmul_(self, Element s):
        r"""
        Return the product ``left * self``.

        INPUT:

        - ``left`` -- an element of the base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: t*(d + 1)  # indirect doctest
            t*d + t

        TESTS:

        We check that :issue:`32210` is fixed::

            sage: A.<t> = GF(5)[]
            sage: S.<X> = A['X', A.derivation()]
            sage: f = 0*X
            sage: f.degree()
            -1
        """
        return self._new_c([s * c for c in self._coeffs], self._parent, 1)

    cpdef _mul_(self, other):
        r"""
        Return the product ``self * right``.

        INPUT:

        - ``right`` -- an Ore polynomial in the same ring as ``self``

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: P = d^2 + t
            sage: Q = d^2 + (t + 1)*d
            sage: P * Q
            d^4 + (t + 1)*d^3 + (t + 2)*d^2 + (t^2 + t)*d
            sage: P * Q == Q * P
            False

        TESTS:

        We check associativity and distributivity::

            sage: U = A.random_element(degree=10)
            sage: V = A.random_element(degree=10)
            sage: W = A.random_element(degree=10)
            sage: U * (V * W) == (U * V) * W
            True
            sage: U * (V + W) == U*V + U*W
            True
            sage: (U + V) * W == U*W + V*W
            True
        """
        cdef coeffs = list((<OrePolynomial_generic_dense>other)._coeffs)
        if coeffs:
            coeffs = self._mul_list(coeffs)
        return self._new_c(coeffs, self._parent, 1)

    cdef _left_quo_rem(self, OrePolynomial other):
        r"""
        Return the quotient and remainder of the left Euclidean
        division of ``self`` by ``other`` (C implementation).
        """
        sig_check()
        cdef list A = list(self._coeffs)
        cdef Py_ssize_t degB = other.degree()
        cdef Morphism m = self._parent.twisting_morphism(-degB)
        cdef RingElement s, inv = other.leading_coefficient()
        cdef Py_ssize_t nb = len(A) - degB
        cdef list L
        cdef list quo = []
        if not inv.is_unit():
            raise ValueError("the leading coefficient of the divisor in not a unit")
        inv = inv.inverse_of_unit()
        for i in range(nb - 1, -1, -1):
            s = A[i + degB]
            if s:
                s = s * inv
                if m is not None:
                    s = m(s)
                L = (<OrePolynomial_generic_dense>other)._mul_list([s])
                for j in range(len(L) - 1):
                    A[i + j] -= L[j]
            quo.append(s)
        quo.reverse()
        return self._new_c(quo, self._parent), self._new_c(A[:degB], self._parent, 1)

    cdef _right_quo_rem(self, OrePolynomial other):
        r"""
        Return the quotient and remainder of the right Euclidean
        division of ``self`` by ``other`` (C implementation).
        """
        sig_check()
        cdef list A = list(self._coeffs)
        cdef Py_ssize_t i, j
        cdef Py_ssize_t degB = other.degree()
        cdef Py_ssize_t nb = len(A) - degB
        cdef Morphism m = self._parent._morphism
        cdef d = self._parent._derivation
        cdef RingElement s, inv = other.leading_coefficient()
        cdef list XnB = (<OrePolynomial_generic_dense>other)._coeffs
        cdef list XnBs = []
        cdef list quo = []
        if not inv.is_unit():
            raise ValueError("the leading coefficient of the divisor in not a unit")
        inv = inv.inverse_of_unit()
        for i in range(nb):
            if i:
                XnB = list(XnB)
                lmul_gen(XnB, m, d)
            XnBs.append(XnB)
        for i in range(nb - 1, -1, -1):
            s = A[i + degB]
            if s:
                if m:
                    s *= self._parent.twisting_morphism(i)(inv)
                else:
                    s *= inv
                XnB = XnBs[i]
                for j in range(len(XnB) - 1):
                    A[j] -= s * XnB[j]
            quo.append(s)
        quo.reverse()
        return self._new_c(quo, self._parent), self._new_c(A[:degB], self._parent, 1)

    cpdef list coefficients(self, sparse=True):
        r"""
        Return the coefficients of the monomials appearing in ``self``.

        If ``sparse=True`` (the default), return only the nonzero coefficients.
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

    def hilbert_shift(self, s, var=None):
        r"""
        Return this Ore polynomial with variable shifted by `s`,
        i.e. if this Ore polynomial is `P(x)`, return `P(x+s)`.

        INPUT:

        - ``s`` -- an element in the base ring

        - ``var`` -- string; the variable name

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]

            sage: L = d^3 + t*d^2
            sage: L.hilbert_shift(t)
            d^3 + 4*t*d^2 + (5*t^2 + 3)*d + 2*t^3 + 4*t
            sage: (d+t)^3 + t*(d+t)^2
            d^3 + 4*t*d^2 + (5*t^2 + 3)*d + 2*t^3 + 4*t

        One can specify another variable name::

            sage: L.hilbert_shift(t, var='x')
            x^3 + 4*t*x^2 + (5*t^2 + 3)*x + 2*t^3 + 4*t

        When the twisting morphism is not trivial, the output lies
        in a different Ore polynomial ring::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: P = x^2 + a*x + a^2
            sage: Q = P.hilbert_shift(a); Q
            x^2 + (2*a^2 + a + 4)*x + a^2 + 3*a + 4
            sage: Q.parent()
            Ore Polynomial Ring in x over
             Finite Field in a of size 5^3 twisted by a |--> a^5 and a*([a |--> a^5] - id)
            sage: Q.parent() is S
            False

        This behavior ensures that the Hilbert shift by a fixed element
        defines a homomorphism of rings::

            sage: # needs sage.rings.finite_rings
            sage: U = S.random_element(degree=5)
            sage: V = S.random_element(degree=5)
            sage: s = k.random_element()
            sage: (U+V).hilbert_shift(s) == U.hilbert_shift(s) + V.hilbert_shift(s)
            True
            sage: (U*V).hilbert_shift(s) == U.hilbert_shift(s) * V.hilbert_shift(s)
            True

        We check that shifting by an element and then by its opposite
        gives back the initial Ore polynomial::

            sage: # needs sage.rings.finite_rings
            sage: P = S.random_element(degree=10)
            sage: s = k.random_element()
            sage: P.hilbert_shift(s).hilbert_shift(-s) == P
            True
        """
        from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
        parent = self._parent
        k = parent.base_ring()
        morphism = parent._morphism
        derivation = parent._derivation
        if morphism is not None:
            if derivation is None:
                derivation = k.derivation(s, twist=morphism)
            else:
                derivation += k.derivation(s, twist=morphism)
        if var is None:
            var = parent.variable_name()
        if derivation is None:
            S = OrePolynomialRing(k, morphism, var)
        else:
            S = OrePolynomialRing(k, derivation, var)
        if not self:
            return S.zero()
        X = S.gen() + s
        Xi = S.one()
        ans = S(self[0])
        for i in range(1, self.degree() + 1):
            Xi = X * Xi
            ans += self[i] * Xi
        return ans

    def quotient_module(self, names=None):
        r"""
        Return the quotient ring `A/AP` as a module over `A`
        where `A` denotes the underlying Ore polynomial ring
        and `P` denotes this Ore polynomial.

        INPUT:

        - ``names`` (default: ``None``) -- a string or a list
          of string, the names of the vector of the canonical
          basis

        EXAMPLES::

        """
        from sage.modules.ore_module import OreModule
        coeffs = self.right_monic().list()
        f = companion_matrix(coeffs, format='bottom')
        return OreModule(f, self.parent(), names=names)


cdef class ConstantOrePolynomialSection(Map):
    r"""
    Representation of the canonical homomorphism from the constants of an Ore
    polynomial ring to the base ring.

    This class is necessary for automatic coercion from zero-degree Ore
    polynomial ring into the base ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.ore_polynomial_element import ConstantOrePolynomialSection
        sage: R.<t> = QQ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: m = ConstantOrePolynomialSection(S, R); m
        Generic map:
            From: Ore Polynomial Ring in x over Univariate Polynomial Ring in t
                  over Rational Field twisted by t |--> t + 1
            To:   Univariate Polynomial Ring in t over Rational Field
    """
    cpdef Element _call_(self, x):
        r"""
        Return the corresponding element of the base ring if ``self`` is a
        constant Ore polynomial. Otherwise, it fails.

        TESTS::

            sage: from sage.rings.polynomial.ore_polynomial_element import ConstantOrePolynomialSection
            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: m = ConstantOrePolynomialSection(S, R); m
            Generic map:
                From: Ore Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
                To:   Univariate Polynomial Ring in t over Rational Field
            sage: m(S([0,1])-S([0,1]))
            0
            sage: m(S([3,1])-S([0,1]))
            3
            sage: m(S([0,1])-S([0,t]))
            Traceback (most recent call last):
            ...
            TypeError: (-t + 1)*x is not a constant polynomial
        """
        if x.degree() <= 0:
            try:
                return <Element>(x.constant_coefficient())
            except AttributeError:
                return <Element>((<OrePolynomial>x).constant_coefficient())
        else:
            raise TypeError(f"{x} is not a constant polynomial")


cdef class OrePolynomialBaseringInjection(Morphism):
    r"""
    Representation of the canonical homomorphism from a ring `R` into an Ore
    polynomial ring over `R`.

    This class is necessary for automatic coercion from the base ring to the Ore
    polynomial ring.

    .. SEEALSO::

        :class:`~sage.rings.polynomial.polynomial_element.PolynomialBaseringInjection`

    EXAMPLES::

        sage: R.<t> = QQ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: S.coerce_map_from(S.base_ring()) #indirect doctest
        Ore Polynomial base injection morphism:
          From: Univariate Polynomial Ring in t over Rational Field
          To:   Ore Polynomial Ring in x over Univariate Polynomial Ring in t
                over Rational Field twisted by t |--> t + 1
    """
    def __init__(self, domain, codomain):
        r"""
        Construct a Skew Polynomial Basering Injection.

        INPUT:

        - ``domain`` -- a ring `R`; this will be the domain of the injection

        - ``codomain`` -- an Ore polynomial ring over ``domain``; this will be
          the codomain

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.polynomial.ore_polynomial_element import OrePolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: OrePolynomialBaseringInjection(k, k['x', Frob])
            Ore Polynomial base injection morphism:
              From: Finite Field in t of size 5^3
              To:   Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: R.<t> = QQ[]
            sage: OrePolynomialBaseringInjection(QQ, k['x', Frob])
            Traceback (most recent call last):
            ...
            AssertionError: the domain of the injection must be the base ring of the Ore polynomial ring
        """
        assert codomain.base_ring() is domain, \
            "the domain of the injection must be the base ring of the Ore polynomial ring"
        Morphism.__init__(self, Hom(domain, codomain))
        self._an_element = codomain.gen()
        self._repr_type_str = "Ore Polynomial base injection"
        self._new_constant_poly_ = self._an_element._new_constant_poly

    def an_element(self):
        r"""
        Return an element of the codomain of the ring homomorphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.polynomial.ore_polynomial_element import OrePolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: m = OrePolynomialBaseringInjection(k, k['x', Frob])
            sage: m.an_element()
            x
        """
        return self._an_element

    cpdef Element _call_(self, e):
        r"""
        Return the corresponding Ore polynomial to the element from the
        base ring according to ``self``.

        INPUT:

        - ``e`` -- element belonging to the base ring according to ``self``

        OUTPUT: the Ore polynomial corresponding to `e` according to ``self``

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.polynomial.ore_polynomial_element import OrePolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: m = OrePolynomialBaseringInjection(k, k['x', Frob])
            sage: m(4)
            4
            sage: parent(m(4))
            Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        try:
            return self._codomain._element_constructor_(e)
        except AttributeError:
            return self._codomain(e)

    def section(self):
        r"""
        Return the canonical homomorphism from the constants of an Ore
        polynomial ring to the base ring according to ``self``.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.polynomial.ore_polynomial_element import OrePolynomialBaseringInjection
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: m = OrePolynomialBaseringInjection(k, k['x', Frob])
            sage: m.section()
            Generic map:
              From: Ore Polynomial Ring in x over Finite Field in t of size 5^3
                    twisted by t |--> t^5
              To:   Finite Field in t of size 5^3
        """
        return ConstantOrePolynomialSection(self._codomain, self.domain())
