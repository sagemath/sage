r"""
Elements of Drinfeld modular forms rings

This module defines the elements of the class
:class:`~sage.modular.drinfeld_modform.ring.DrinfeldModularForms`.

AUTHORS:

- David Ayotte (2022): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 DAVID AYOTTE <da.ayotte@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial import MPolynomial


class DrinfeldModularFormsElement(ModuleElement):
    r"""
    Element class of rings of Drinfeld modular forms.

    Recall that a *graded Drinfeld form* is a sum of Drinfeld modular
    forms having potentially different weights:

    .. MATH::

        F = f_{k_1} + f_{k_2} + \cdots + f_{k_n}

    where `f_{k_i}` is a Drinfeld modular form of weight `k_i`. We also
    say that `f_{k_i}` is an *homogeneous component of weight* `k_i`. If
    `n=1`, then we say that `F` is *homogeneous of weight* `k_1`.

    EXAMPLES: use the ``inject_variable`` method of the parent to
    quickly assign variables names to the generators::

        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: M = DrinfeldModularForms(K, 2)
        sage: M.inject_variables()
        Defining g1, g2
        sage: g1 in M
        True
        sage: g2.parent()
        Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

    Next, via algebraic combination of the generator, we may create any
    element of the ring::

        sage: F = g1*g2 + g2
        sage: F
        g1*g2 + g2
        sage: F.is_homogeneous()
        False
        sage: F.homogeneous_components()
        {8: g2, 10: g1*g2}

    If the created form is homogeneous, we can ask for its weight in
    which case it will be a Drinfeld modular form::

        sage: H = g1^4*g2^9 + T*g1^8*g2^8 + (T^2 - 1)*g1^28*g2^3
        sage: H.is_homogeneous()
        True
        sage: H.weight()
        80

    You can also construct an element by simply passing a multivariate
    polynomial to the parent::

        sage: f1, f2 = polygens(K, 2, 'f1, f2')
        sage: M(f1)
        g1
        sage: M(f2)
        g2
        sage: M(T*f1 + f2^3 + T^2 + 1)
        g2^3 + T*g1 + (T^2 + 1)

    .. NOTE::

        This class should not be directly instanciated, instead create
        an instance of the parent
        :class:`~sage.modular.drinfeld_modform.ring.DrinfeldModularForms`
        and access its elements using the relevant methods.
    """
    def __init__(self, parent, polynomial):
        if not isinstance(polynomial, MPolynomial):
            raise TypeError("input must be a multivariate polynomial")
        if not parent.base_ring().has_coerce_map_from(polynomial.base_ring()):
            raise ValueError("unable to coerce base ring of the given "
                             "polynomial into Drinfeld modular form ring")
        poly = parent._poly_ring(polynomial)
        self._polynomial = poly

        super().__init__(parent)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

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
        Return the addition of ``self`` with ``other``.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.0 + M.1  # indirect doctest
            g2 + g1
        """
        return self.__class__(self.parent(), self._polynomial + other._polynomial)

    def _mul_(self, other):
        r"""
        Return the multiplication of ``self`` with ``other``.

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
        return self.__class__(self.parent(), self._polynomial*other._polynomial)

    def _lmul_(self, c):
        r"""
        Return the scalar multiplication of ``self`` by `c`.

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
        return self.__class__(self.parent(), c*self._polynomial)

    def __neg__(self):
        r"""
        Return the negation of ``self``.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: -M.0  # indirect doctest
            -g1
        """
        return self.__class__(self.parent(), -self._polynomial)

    def __bool__(self):
        r"""
        Return ``True`` whether ``self`` is nonzero.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: bool(M.0)
            True
        """
        return bool(self._polynomial)

    def _latex_(self):
        r"""
        Return the LaTeX expression of ``self``.

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
        Return the comparison of ``self`` with ``other``.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.0 == M.1
            False
            sage: M.0 != M.1
            True
            sage: M.0 == M.0
            True
            sage: M.0 < M.1
            Traceback (most recent call last):
            ...
            TypeError: '<' not supported between instances of 'DrinfeldModularForms_with_category.element_class' and 'DrinfeldModularForms_with_category.element_class'
        """
        if op != op_EQ and op != op_NE:
            return NotImplemented
        return richcmp(self._polynomial, other._polynomial, op)

    def rank(self):
        r"""
        Return the rank of this graded Drinfeld form.

        Note that the rank is independent of the chosen form and depends
        only on the parent.

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
        Return ``True`` whether this graded Drinfeld form is the
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
        Return ``True`` whether this graded Drinfeld form is the
        additive identity.

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

    def is_homogeneous(self):
        r"""
        Return whether the graded form is homogeneous in the weight.

        We recall that elements of Drinfeld modular forms ring are not
        necessarily modular forms as they may have mixed weight
        components.

        EXAMPLES::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2)
            sage: M.inject_variables()
            Defining g1, g2
            sage: f = g1^5*g2^2  # homogeneous polynomial
            sage: f.is_homogeneous()
            True
            sage: g = g1 + g2  # mixed weight components
            sage: g.is_homogeneous()
            False
        """
        return self._polynomial.is_homogeneous()

    def homogeneous_components(self):
        r"""
        Return the homogeneous components of this graded Drinfeld
        form.

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
        M = self.parent()
        components = self._polynomial.homogeneous_components().items()
        return {k: self.__class__(M, p) for k, p in components}

    def polynomial(self):
        r"""
        Return this graded Drinfeld forms as a multivariate polynomial
        over the generators of the ring.

        OUTPUT: a multivariate polynomial over the base ring

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

    def type(self):
        r"""
        Return the type of this graded Drinfeld form.

        Recall that the *type* is the integer `0 \leq m \leq q-1` such that

        .. MATH::

            f(\gamma(w)) = \mathrm{det}(\gamma)^m j(\gamma, w)^k f(w).

        EXAMPLES::

            sage: A = GF(11)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2, has_type=True)
            sage: M.inject_variables()
            Defining g1, h2
            sage: F = g1*h2^9
            sage: F.type()
            9
            sage: (h2^11).type()
            1
            sage: g1.type()
            0

        The type only makes sense when the form is homogeneous::

            sage: F = g1^4 + h2
            sage: F.type()
            Traceback (most recent call last):
            ...
            ValueError: the graded form is not homogeneous

        TESTS::

            sage: A = GF(2)['T']; K = Frac(A)
            sage: M = DrinfeldModularForms(K, 2, has_type=False)
            sage: (M.1).type()
            0
        """
        if not self.is_homogeneous():
            raise ValueError("the graded form is not homogeneous")
        if not self.parent()._has_type:
            return ZZ(0)
        q = self.base_ring().base_ring().cardinality()
        return self.polynomial().degrees()[-1] % (q-1)

    def weight(self):
        r"""
        Return the weight of this graded Drinfeld modular form.

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

        If the form is not homogeneous, then the method returns an
        error::

            sage: f = g1 + g2
            sage: f.weight()
            Traceback (most recent call last):
            ...
            ValueError: the graded form is not homogeneous
        """
        if not self.is_homogeneous():
            raise ValueError("the graded form is not homogeneous")
        return self._polynomial.degree()
