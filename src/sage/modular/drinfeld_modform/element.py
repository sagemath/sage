r"""
This module defines the elements of the class
:class:`~drinfeld_modular_forms.ring.DrinfeldModularForms`.

EXAMPLES::

    sage: q = 3
    sage: A = GF(q)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularForms(K, 2)  # rank 2
    sage: M.inject_variables()
    Defining g1, g2
    sage: g1.parent()
    Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

A *graded Drinfeld modular form* is a sum of modular forms having
potentially different weights::

    sage: F = g1*g2 + g2
    sage: F.is_drinfeld_modular_form()
    False
    sage: F.homogeneous_components()
    {8: g2, 10: g1*g2}
    sage: H = g1^4*g2^9 + T*g1^8*g2^8 + (T^2 - 1)*g1^28*g2^3
    sage: H.is_drinfeld_modular_form()
    True
    sage: H.weight()
    80

AUTHORS:

- David Ayotte (2022): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 DAVID AYOTTE <davidayotte94@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from sage.rings.polynomial.multi_polynomial import MPolynomial

class DrinfeldModularFormsElement(ModuleElement):
    r"""
    Element class of rings of Drinfeld modular forms.

    EXAMPLES::

        sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
        sage: M = DrinfeldModularForms(K, 2)
        sage: M.inject_variables()
        Defining g1, g2
        sage: (T^2 + 1)*(g1 + g1*g2)
        (T^2 + 1)*g1*g2 + (T^2 + 1)*g1
        sage: (g1).parent()
        Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
        sage: g2 in M
        True
    """
    def __init__(self, parent, polynomial):
        if not isinstance(polynomial, MPolynomial):
            raise TypeError("input must be a multivariate polynomial")
        if not parent.base_ring().has_coerce_map_from(polynomial.base_ring()):
            raise ValueError("unable to coerce base ring of the given "
                             "polynomial into Drinfeld modular form ring")
        poly = parent._poly_ring(polynomial)
        self._polynomial = poly

        ModuleElement.__init__(self, parent)

    def _repr_(self):
        r"""
        Return the string representation of self.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: (M.0)._repr_()
            'g1'
            sage: M.0 + M.1
            g2 + g1
        """
        return str(self._polynomial)

    def _add_(self, other):
        r"""
        Return the addition of self with other.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.0 + M.1  # indirect doctest
            g2 + g1
        """
        return self.__class__(self.parent(), self._polynomial + other._polynomial)

    def _mul_(self, other):
        r"""
        Return the multiplication of self with other.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.0*M.1  # indirect doctest
            g1*g2
            sage: M.0*(M.0 + M.1)
            g1*g2 + g1^2
            sage: (M.0 + M.1)*M.0
            g1*g2 + g1^2
        """
        return self.__class__(self.parent(), self._polynomial * other._polynomial)

    def _lmul_(self, c):
        r"""
        Return the scalar multiplication of self by `c`.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularForms(K, 2)
            sage: (T^2 + T + 2) * M.0  # indirect doctest
            (T^2 + T - 1)*g1
            sage: M.1 * (T^5 + T^2)
            (T^5 + T^2)*g2
            sage: 0 * M.1
            0
            sage: M.0 * 0
            0
        """
        return self.__class__(self.parent(), c * self._polynomial)

    def __neg__(self):
        r"""
        Return the negation of self.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: -M.0  # indirect doctest
            -g1
        """
        return self.__class__(self.parent(), -self._polynomial)

    def __bool__(self):
        r"""
        Return True whether self is nonzero.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: bool(M.0)
            True
        """
        return bool(self._polynomial)

    def _latex_(self):
        r"""
        Return the latex expression of self.

        TESTS::

            sage: A = GF(3)['T']; K.<T> = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.inject_variables()
            Defining g1, g2
            sage: latex(g1)
            g_{1}
            sage: latex(g2)
            g_{2}
            sage: latex(1/T*g1^5 + g2*g1)
            \frac{1}{T} g_{1}^{5} + g_{1} g_{2}
        """
        return self._polynomial._latex_()

    def _richcmp_(self, other, op):
        r"""
        Return the comparison of self with other.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.0 == M.1
            False
            sage: M.0 != M.1
            True
            sage: M.0 == M.0
            True
        """
        if op != op_EQ and op != op_NE:
            raise TypeError('invalid comparison between modular forms ring elements')
        return richcmp(self._polynomial, other._polynomial, op)

    def rank(self):
        r"""
        Return the rank of the graded Drinfeld form.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M2 = DrinfeldModularForms(K, 2)
            sage: (M2.0).rank()
            2
            sage: M5 = DrinfeldModularForms(K, 5)
            sage: (M5.0 + M5.3).rank()
            5
        """
        return self.parent()._rank

    def is_one(self):
        r"""
        Return ``True`` whether the given graded Drinfeld form is the
        multiplicative identity.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: u = M.one()
            sage: u.is_one()
            True
            sage: (M.0).is_one()
            False
        """
        return self._polynomial.is_one()

    def is_zero(self):
        r"""
        Return ``True`` whether the given graded Drinfeld form is the additive
        identity.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: z = M.zero()
            sage: z.is_zero()
            True
            sage: f = M.0
            sage: f.is_zero()
            False
            sage: (f - f).is_zero()
            True
            sage: (0 * M.0).is_zero()
            True
        """
        return not bool(self)

    def is_drinfeld_modular_form(self):
        r"""
        Return whether ``self`` is a Drinfeld modular form.

        We recall that elements of Drinfeld modular forms ring are not
        necessarily modular forms as they may have mixed weight components.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.inject_variables()
            Defining g1, g2
            sage: f = g1^5*g2^2  # homogeneous polynomial
            sage: f.is_drinfeld_modular_form()
            True
            sage: g = g1 + g2  # mixed weight components
            sage: g.is_drinfeld_modular_form()
            False
        """
        return self._polynomial.is_homogeneous()

    def homogeneous_components(self):
        r"""
        Return the homogeneous components of this graded Drinfeld
        modular form.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.inject_variables()
            Defining g1, g2
            sage: F = g1 + g1^2 + g1*g2^2 + g2^4
            sage: D = F.homogeneous_components(); D
            {2: g1, 4: g1^2, 18: g1*g2^2, 32: g2^4}
            sage: D[32]
            g2^4
        """
        elt_class = self.__class__
        M = self.parent()
        components = self._polynomial.homogeneous_components().items()
        return {k: elt_class(M, p) for k, p in components}

    def polynomial(self):
        r"""
        Return self as a multivariate polynomial over the generators of
        the ring.

        OUTPUT:

        A multivariate polynomial over the base ring.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.inject_variables()
            Defining g1, g2
            sage: P1 = g1.polynomial();
            sage: P2 = g2.polynomial();
            sage: P2^2 + P1^2 + P1
            g2^2 + g1^2 + g1
            sage: P1.parent()
            Multivariate Polynomial Ring in g1, g2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

        The degree of each variables corresponds to the weight of the
        generator::

            sage: P1.degree()
            2
            sage: P2.degree()
            8
        """
        return self._polynomial

    def type_m(self):  # Find a better name
        r"""
        Return the type of the graded Drinfeld form.

        Recall that the *type* is the integer `m` such that

        .. MATH::

            f(\gamma(w)) = \mathrm{det}(\gamma)^m j(\gamma, w)^k f(w).

        This method is only implemented when the rank is two.

        EXAMPLES::

            sage: A = GF(11)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2, has_type=True)
            sage: M.inject_variables()
            Defining g1, h
            sage: F = g1*h^9
            sage: F.type_m()
            9
            sage: (h^11).type_m()
            1
            sage: g1.type_m()
            0
        """
        if not self.is_drinfeld_modular_form():
            raise ValueError("self should be a Drinfeld modular form")
        if not self.parent()._has_type:
            return ZZ(0)
        q = self.base_ring().base_ring().cardinality()
        return self.polynomial().degrees()[-1]%(q-1)

    def weight(self):
        r"""
        Return the weight of self.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.inject_variables()
            Defining g1, g2
            sage: g1.weight()
            2
            sage: g2.weight()
            8
            sage: f = g1^5*g2^2
            sage: f.weight()
            26

        If the form is not modular, then the method returns an error::

            sage: f = g1 + g2
            sage: f.weight()
            Traceback (most recent call last):
            ...
            ValueError: the given ring element is not a Drinfeld modular form
        """
        if not self.is_drinfeld_modular_form():
            raise ValueError("the given ring element is not a Drinfeld modular form")
        return self._polynomial.degree()
