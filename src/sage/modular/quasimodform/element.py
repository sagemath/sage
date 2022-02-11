"""
Elements of quasimodular forms rings

AUTHORS:

- DAVID AYOTTE (2021-03-18): initial version
"""
# ****************************************************************************
#       Copyright (C) 2021 David Ayotte
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modular.modform.eis_series import eisenstein_series_qexp
from sage.modular.modform.element import GradedModularFormElement

from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from sage.rings.polynomial.polynomial_element import Polynomial

class QuasiModularFormsElement(ModuleElement):
    r"""
    A quasimodular forms ring element. Such an element is describbed by SageMath
    as a polynomial

    .. MATH::

        f_0 + f_1 E_2 + f_2 E_2^2 + \cdots + f_m E_2^m

    where each `f_i` a graded modular form element
    (see :class:`~sage.modular.modform.element.GradedModularFormElement`)

    EXAMPLES::

        sage: QM = QuasiModularForms(1)
        sage: QM.gens()
        [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
        1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        sage: QM.0 + QM.1
        2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
        sage: QM.0 * QM.1
        1 + 216*q - 3672*q^2 - 62496*q^3 - 322488*q^4 - 1121904*q^5 + O(q^6)
        sage: (QM.0)^2
        1 - 48*q + 432*q^2 + 3264*q^3 + 9456*q^4 + 21600*q^5 + O(q^6)
        sage: QM.0 == QM.1
        False

    Quasimodular forms ring element can be created via a polynomial in `E2` over the ring of modular forms::

        sage: E2 = QM.polygen()
        sage: E2.parent()
        Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
        sage: QM(E2)
        1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
        sage: M = QM.modular_forms_subring()
        sage: QM(M.0 * E2 + M.1 * E2^2)
        2 - 336*q + 4320*q^2 + 398400*q^3 - 3772992*q^4 - 89283168*q^5 + O(q^6)
    """
    def __init__(self, parent, polynomial):
        r"""
        INPUT:

        - ``parent`` - a quasimodular forms ring
        - ``polynomial`` - a polynomial `f_0 + f_1 E_2 + ... + f_n E_2^n` where
          each `f_i` are modular forms ring elements and `E_2` correspond to the
          weight 2 Eisenstein series

        OUTPUT:

        ``QuasiModularFormsElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.element_class(QM, 'E2')
            Traceback (most recent call last):
            ...
            TypeError: 'polynomial' argument should be of type 'Polynomial'
            sage: x = polygen(QQ)
            sage: QM.element_class(QM, x^2 + 1)
            Traceback (most recent call last):
            ...
            ValueError: at least one coefficient is not a 'GradedModularFormElement'
        """
        if not isinstance(polynomial, Polynomial):
            raise TypeError("'polynomial' argument should be of type 'Polynomial'")
        for f in polynomial.coefficients():
            if not isinstance(f, GradedModularFormElement):
                raise ValueError("at least one coefficient is not a 'GradedModularFormElement'")
        self._polynomial = polynomial
        ModuleElement.__init__(self, parent)

    def q_expansion(self, prec=6):
        r"""
        Computes the `q`-expansion of self to precision `prec`.

        An alias of this method is ``qexp``.

        EXAMPLES::

            sage: QM = QuasiModularForms()
            sage: E2 = QM.0
            sage: E2.q_expansion()
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: E2.q_expansion(prec=10)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 - 288*q^6 - 192*q^7 - 360*q^8 - 312*q^9 + O(q^10)
        """
        E2 = eisenstein_series_qexp(2, prec=prec, K=self.base_ring(), normalization='constant') #normalization -> to force integer coefficients
        coefficients = self._polynomial.coefficients(sparse=False)
        return sum(f.q_expansion(prec=prec)*E2**idx for idx, f in enumerate(coefficients))

    qexp = q_expansion # alias

    def _repr_(self):
        r"""
        String representation of self.

        TESTS::

            sage: QM = QuasiModularForms()
            sage: QM.0
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.1
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: QM.2
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
        """
        return str(self.q_expansion())

    def _richcmp_(self, other, op):
        r"""
        Compare self with other.

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.0 == QM.1
            False
            sage: QM.0 == QM.0
            True
            sage: QM.0 != QM.1
            True
            sage: QM.0 != QM.0
            False
            sage: QM.0 < QM.1
            Traceback (most recent call last):
            ...
            TypeError: invalid comparison between quasimodular forms ring elements
        """
        if op != op_EQ and op != op_NE:
            raise TypeError('invalid comparison between quasimodular forms ring elements')
        return richcmp(self._polynomial, other._polynomial, op)

    def _add_(self, other):
        r"""
        Addition of two ``QuasiModularFormElement``.

        INPUT:

        - ``other`` - ``QuasiModularFormElement``

        OUTPUT: a ``QuasiModularFormElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.0 + QM.1
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
            sage: QM.0 + (QM.1 + QM.2) == (QM.0 + QM.1) + QM.2
            True
            sage: QM = QuasiModularForms(5)
            sage: QM.0 + QM.1 + QM.2 + QM.3
            3 - 17*q - 54*q^2 - 62*q^3 - 98*q^4 + 137*q^5 + O(q^6)
        """
        return self.__class__(self.parent(), self._polynomial + other._polynomial)

    def __neg__(self):
        r"""
        The negation of ``self```

        TESTS::

            sage: -QuasiModularForms(1).0
            -1 + 24*q + 72*q^2 + 96*q^3 + 168*q^4 + 144*q^5 + O(q^6)
            sage: QuasiModularForms(1).0 - QuasiModularForms(1).0
            0
            sage: -QuasiModularForms(Gamma1(2)).2
            -1 - 240*q^2 - 2160*q^4 + O(q^6)
        """
        return self.__class__(self.parent(), -self._polynomial)

    def _mul_(self, other):
        r"""
        The multiplication of two ``QuasiModularFormElement``

        INPUT:

        - ``other`` - ``QuasiModularFormElement``

        OUTPUT: a ``QuasiModularFormElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.0 * QM.1
            1 + 216*q - 3672*q^2 - 62496*q^3 - 322488*q^4 - 1121904*q^5 + O(q^6)
            sage: (QM.0 * QM.1) * QM.2 == QM.0 * (QM.1 * QM.2)
            True
            sage: QM = QuasiModularForms(Gamma1(5))
            sage: QM.0 * QM.1 * QM.2
            q - 24*q^2 - 66*q^3 - 189*q^4 - 1917*q^5 + O(q^6)
        """
        return self.__class__(self.parent(), self._polynomial * other._polynomial)

    def _lmul_(self, c):
        r"""
        The left action of the base ring on self.

        INPUT:

        - ``other`` - ``QuasiModularFormElement``

        OUTPUT: a ``QuasiModularFormElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: (1/2) * QM.0
            1/2 - 12*q - 36*q^2 - 48*q^3 - 84*q^4 - 72*q^5 + O(q^6)
            sage: QM.0 * (3/2)
            3/2 - 36*q - 108*q^2 - 144*q^3 - 252*q^4 - 216*q^5 + O(q^6)
            sage: (5/2) * QuasiModularForms(Gamma0(7)).0 * (3/2)
            15/4 - 90*q - 270*q^2 - 360*q^3 - 630*q^4 - 540*q^5 + O(q^6)
        """
        return self.__class__(self.parent(), c * self._polynomial)

    def __bool__(self):
        r"""
        Return "True" if ``self`` is non-zero and "False" otherwise.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: bool(QM(0))
            False
            sage: bool(QM(1))
            True
            sage: bool(QM.0)
            True
        """
        return bool(self._polynomial)

    def is_zero(self):
        r"""
        Return "True" if the quasiform is 0 and "False" otherwise

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.zero().is_zero()
            True
            sage: QM(0).is_zero()
            True
            sage: QM(1/2).is_zero()
            False
            sage: (QM.0).is_zero()
            False
        """
        return not self

    def is_one(self):
        r"""
        Return "True" if the quasiform is 1 and "False" otherwise

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.one().is_one()
            True
            sage: QM(1).is_one()
            True
            sage: (QM.0).is_one()
            False
        """
        return self._polynomial.is_one()

    def is_graded_modular_form(self):
        r"""
        Return ``True`` if the given quasimodular form is a graded modular forms element
        and ``False`` otherwise.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0).is_graded_modular_form()
            False
            sage: (QM.1).is_graded_modular_form()
            True
            sage: (QM.1 + QM.0^2).is_graded_modular_form()
            False
            sage: (QM.1^2 + QM.2).is_graded_modular_form()
            True
            sage: QM = QuasiModularForms(Gamma0(6))
            sage: (QM.0).is_graded_modular_form()
            False
            sage: (QM.1 + QM.2 + QM.1 * QM.3).is_graded_modular_form()
            True
            sage: QM.zero().is_graded_modular_form()
            True

        .. NOTE::

            A graded modular form in SageMath is not necessarily a modular form
            as it can have mixed weight components. To check for modular forms
            only, see the method :meth:`is_modular_form`.
        """
        return self._polynomial.degree() <= 0

    def is_modular_form(self):
        r"""
        Return ``True`` if the given quasimodular form is a modular form and
        ``False`` otherwise.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0).is_modular_form()
            False
            sage: (QM.1).is_modular_form()
            True
            sage: (QM.1 + QM.2).is_modular_form() # mixed weight components
            False
            sage: QM.zero().is_modular_form()
            True
        """
        return self._polynomial.degree() <= 0 and self._polynomial[0].is_modular_form()

    def polynomial(self, names='E2, E4, E6'):
        r"""
        Return a multivariate polynomial `P(E_2, E_4, E_6)` corresponding to the
        given form where `E_2`, `E_4` and `E_6` are the generators of the
        quasimodular form ring given by
        :meth:`~sage.modular.quasiform.ring.QuasiModularForms.gens`.

        INPUT:

        - ``names`` (str, default: ``'E2, E4, E6'``) -- a list or tuple of names
          (strings), or a comma separated string. Correspond to the names of the
          variables;

        OUTPUT: A multivariate polynomial in the variables ``names``

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0 + QM.1).polynomial()
            E4 + E2
            sage: (1/2 + QM.0 + 2*QM.1^2 + QM.0*QM.2).polynomial()
            E2*E6 + 2*E4^2 + E2 + 1/2
        """
        P = self.parent().polynomial_ring(names)
        g0, g1 = self.parent().modular_forms_subring().polynomial_ring(names='x').gens()
        E2, E4, E6 = P.gens()
        return sum(f.to_polynomial().subs({g0:E4, g1:E6}) * E2 ** exp for exp, f in enumerate(self._polynomial.coefficients(sparse=False)))

    to_polynomial = polynomial # alias

    def weights_list(self):
        r"""
        Return the list of the weights of all the graded components of the given
        graded quasimodular form.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0).weights_list()
            [2]
            sage: (QM.0 + QM.1 + QM.2).weights_list()
            [2, 4, 6]
            sage: (QM.0 * QM.1 + QM.2).weights_list()
            [6]
            sage: QM(1/2).weights_list()
            [0]
        """
        return sorted(list(self.to_polynomial().homogeneous_components()))

    def is_homogeneous(self):
        r"""
        Return True if the graded quasimodular form is a homogeneous element,
        that is it lives in a unique graded components of the graded ring of
        self.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0).is_homogeneous()
            True
            sage: (QM.0 + QM.1).is_homogeneous()
            False
            sage: (QM.0 * QM.1 + QM.2).is_homogeneous()
            True
            sage: QM(1).is_homogeneous()
            True
            sage: (1 + QM.0).is_homogeneous()
            False
        """
        return len(self.weights_list()) == 1

    def weight(self):
        r"""
        Return the weight of the given quasimodular form.

        Note that the given form must be homogeneous.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0).weight()
            2
            sage: (QM.0 * QM.1 + QM.2).weight()
            6
            sage: QM(1/2).weight()
            0
            sage: (QM.0 + QM.1).weight()
            Traceback (most recent call last):
            ...
            ValueError: the given graded quasiform is not an homogeneous element
        """
        if self.is_homogeneous():
            return self.to_polynomial().degree()
        else:
            raise ValueError("the given graded quasiform is not an homogeneous element")

    def homogeneous_components(self):
        r"""
        Return a dictionary where the values are the homogeneous components of
        the given graded form and the keys are the weights of those components.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: (QM.0).homogeneous_components()
            {2: 1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)}
            sage: (QM.0 + QM.1 + QM.2).homogeneous_components()
            {2: 1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            4: 1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            6: 1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)}
            sage: (1 + QM.0).homogeneous_components()
            {0: 1, 2: 1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)}
        """
        QM = self.parent()
        poly_self = self.to_polynomial()
        pol_hom_comp = poly_self.homogeneous_components()
        return {k: QM.from_polynomial(pol) for k, pol in pol_hom_comp.items()}

    def serre_derivative(self):
        r"""
        Return the Serre derivative of the given quasimodular form.

        If the form is not homogeneous, then this method sums the Serre
        derivative of each homogeneous component.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: E2, E4, E6 = QM.gens()
            sage: DE2 = E2.serre_derivative(); DE2
            -1/6 - 16*q - 216*q^2 - 832*q^3 - 2248*q^4 - 4320*q^5 + O(q^6)
            sage: DE2 == (-E2^2 - E4)/12
            True
            sage: DE4 = E4.serre_derivative(); DE4
            -1/3 + 168*q + 5544*q^2 + 40992*q^3 + 177576*q^4 + 525168*q^5 + O(q^6)
            sage: DE4 == (-1/3) * E6
            True
            sage: DE6 = E6.serre_derivative(); DE6
            -1/2 - 240*q - 30960*q^2 - 525120*q^3 - 3963120*q^4 - 18750240*q^5 + O(q^6)
            sage: DE6 == (-1/2) * E4^2
            True

        The Serre derivative raises the weight of homogeneous elements by 2::

            sage: F = E6 + E4 * E2
            sage: F.weight()
            6
            sage: F.serre_derivative().weight()
            8
        """
        # initial variables:
        QM = self.parent()
        R = QM.base_ring()
        E2 = QM.gen(0)
        E4 = QM.gen(1)

        # compute the derivative of E2: q*dE2/dq
        E2deriv = R(12).inverse_of_unit() * (E2 ** 2 - E4)

        # sum the Serre derivative of each monomial of the form: f * E2^n
        # they are equal to:
        # [E2^n * serre_deriv(f)]  +  [n * f * E2^(n-1) * D(E2)]  -  [n/6 * f * E2^(n+1)]
        #   =      A               +              B               -           C
        der = QM.zero()
        u6 = R(6).inverse_of_unit()
        for n, f in enumerate(self._polynomial.coefficients(sparse=False)):
            if n == 0:
                der += QM(f.serre_derivative())
            else:
                A = (E2 ** n) * f.serre_derivative()
                B = R(n) * f * E2 ** (n - 1) * E2deriv
                C = R(n) * u6 * E2 ** (n + 1) * f
                der += QM(A + B - C)
        return der

    def derivative(self):
        r"""
        Return the derivative `q \frac{d}{dq}` of the given quasimodular form.

        If the form is not homogeneous, then this method sums the derivative of
        each homogeneous component.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: E2, E4, E6 = QM.gens()
            sage: dE2 = E2.derivative(); dE2
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 + O(q^6)
            sage: dE2 == (E2^2 - E4)/12 # Ramanujan identity
            True
            sage: dE4 = E4.derivative(); dE4
            240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + O(q^6)
            sage: dE4 == (E2 * E4 - E6)/3 # Ramanujan identity
            True
            sage: dE6 = E6.derivative(); dE6
            -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 + O(q^6)
            sage: dE6 == (E2 * E6 - E4^2)/2 # Ramanujan identity
            True

        Note that the derivative of a modular form is not necessarily a modular form::

            sage: dE4.is_modular_form()
            False
            sage: dE4.weight()
            6
        """
        QM = self.parent()
        E2 = QM.gen(0)
        R = self.base_ring()
        u = R(12).inverse_of_unit()
        hom_comp = self.homogeneous_components()

        return sum(f.serre_derivative() + R(k) * u * f * E2 for k, f in hom_comp.items())
