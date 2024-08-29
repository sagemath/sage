r"""
Graded quasimodular forms ring

Let `E_2` be the weight 2 Eisenstein series defined by

.. MATH::

    E_2(z) = 1 - \frac{2k}{B_k} \sum_{n=1}^{\infty} \sigma(n) q^n

where `\sigma` is the sum of divisors function and `q = \mathrm{exp}(2\pi i z)`
is the classical parameter at infinity, with `\mathrm{im}(z)>0`. This weight 2
Eisenstein series is not a modular form as it does not satisfy the
modularity condition:

.. MATH::

    z^2 E_2(-1/z) = E_2(z) + \frac{2k}{4\pi i B_k z}.

`E_2` is a quasimodular form of weight 2. General quasimodular forms of given
weight can also be defined. We denote by `QM` the graded ring of quasimodular
forms for the full modular group `\SL_2(\ZZ)`.

The SageMath implementation of the graded ring of quasimodular forms uses the
following isomorphism:

.. MATH::

    QM \cong M_* [E_2]

where `M_* \cong \CC[E_4, E_6]` is the graded ring of modular forms for
`\SL_2(\ZZ)`. (see :class:`sage.modular.modform.ring.ModularFormsRing`).

More generally, if `\Gamma \leq \SL_2(\ZZ)` is a congruence subgroup,
then the graded ring of quasimodular forms for `\Gamma` is given by
`M_*(\Gamma)[E_2]` where `M_*(\Gamma)` is the ring of modular forms for
`\Gamma`.

The SageMath implementation of the graded quasimodular forms ring allows
computation of a set of generators and perform usual arithmetic operations.

EXAMPLES::

    sage: QM = QuasiModularForms(1); QM
    Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
    sage: QM.category()
    Category of commutative graded algebras over Rational Field
    sage: QM.gens()
    [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
     1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
     1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
    sage: E2 = QM.0; E4 = QM.1; E6 = QM.2
    sage: E2 * E4 + E6
    2 - 288*q - 20304*q^2 - 185472*q^3 - 855216*q^4 - 2697408*q^5 + O(q^6)
    sage: E2.parent()
    Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field

The ``polygen`` method also return the weight-2 Eisenstein series as a
polynomial variable over the ring of modular forms::

    sage: QM = QuasiModularForms(1)
    sage: E2 = QM.polygen(); E2
    E2
    sage: E2.parent()
    Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field

An element of a ring of quasimodular forms can be created via a list of modular
forms or graded modular forms. The `i`-th index of the list will correspond to
the `i`-th coefficient of the polynomial in `E_2`::

    sage: QM = QuasiModularForms(1)
    sage: E2 = QM.0
    sage: Delta = CuspForms(1, 12).0
    sage: E4 = ModularForms(1, 4).0
    sage: F = QM([Delta, E4, Delta + E4]); F
    2 + 410*q - 12696*q^2 - 50424*q^3 + 1076264*q^4 + 10431996*q^5 + O(q^6)
    sage: F == Delta + E4 * E2 + (Delta + E4) * E2^2
    True

One may also create rings of quasimodular forms for certain congruence subgroups::

    sage: QM = QuasiModularForms(Gamma0(5)); QM
    Ring of Quasimodular Forms for Congruence Subgroup Gamma0(5) over Rational Field
    sage: QM.ngens()
    4

The first generator is the weight 2 Eisenstein series::

    sage: E2 = QM.0; E2
    1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)

The other generators correspond to the generators given by the method
:meth:`sage.modular.modform.ring.ModularFormsRing.gens`::

    sage: QM.gens()
    [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
     1 + 6*q + 18*q^2 + 24*q^3 + 42*q^4 + 6*q^5 + O(q^6),
     1 + 240*q^5 + O(q^6),
     q + 10*q^3 + 28*q^4 + 35*q^5 + O(q^6)]
    sage: QM.modular_forms_subring().gens()
    [1 + 6*q + 18*q^2 + 24*q^3 + 42*q^4 + 6*q^5 + O(q^6),
     1 + 240*q^5 + O(q^6),
     q + 10*q^3 + 28*q^4 + 35*q^5 + O(q^6)]

It is possible to convert a graded quasimodular form into a polynomial where
each variable corresponds to a generator of the ring::

    sage: QM = QuasiModularForms(1)
    sage: E2, E4, E6 = QM.gens()
    sage: F = E2*E4*E6 + E6^2; F
    2 - 1296*q + 91584*q^2 + 14591808*q^3 + 464670432*q^4 + 6160281120*q^5 + O(q^6)
    sage: p = F.polynomial('E2, E4, E6'); p
    E2*E4*E6 + E6^2
    sage: P = p.parent(); P
    Multivariate Polynomial Ring in E2, E4, E6 over Rational Field

The generators of the polynomial ring have degree equal to the weight of the
corresponding form::

    sage: P.inject_variables()
    Defining E2, E4, E6
    sage: E2.degree()
    2
    sage: E4.degree()
    4
    sage: E6.degree()
    6

This works also for congruence subgroup::

    sage: QM = QuasiModularForms(Gamma1(4))
    sage: QM.ngens()
    5
    sage: QM.polynomial_ring()
    Multivariate Polynomial Ring in E2, E2_0, E2_1, E3_0, E3_1 over Rational Field
    sage: (QM.0 + QM.1*QM.0^2 + QM.3 + QM.4^3).polynomial()
    E3_1^3 + E2^2*E2_0 + E3_0 + E2

One can also convert a multivariate polynomial into a quasimodular form::

    sage: QM.polynomial_ring().inject_variables()
    Defining E2, E2_0, E2_1, E3_0, E3_1
    sage: QM.from_polynomial(E3_1^3 + E2^2*E2_0 + E3_0 + E2)
    3 - 72*q + 396*q^2 + 2081*q^3 + 19752*q^4 + 98712*q^5 + O(q^6)

.. NOTE::

    - Currently, the only supported base ring is the Rational Field;
    - Spaces of quasimodular forms of fixed weight are not yet implemented.

REFERENCE:

See section 5.3 (page 58) of [Zag2008]_

AUTHORS:

- David Ayotte (2021-03-18): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 DAVID AYOTTE
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import product, chain

from sage.categories.graded_algebras import GradedAlgebras

from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0
from sage.modular.arithgroup.congroup_generic import CongruenceSubgroupBase
from sage.modular.modform.element import GradedModularFormElement, ModularFormElement
from sage.modular.modform.space import ModularFormsSpace
from sage.modular.modform.ring import ModularFormsRing

from sage.rings.integer import Integer
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.power_series_poly import PowerSeries_poly
from sage.rings.rational_field import QQ

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from .element import QuasiModularFormsElement


class QuasiModularForms(Parent, UniqueRepresentation):
    r"""
    The graded ring of quasimodular forms for the full modular group
    `\SL_2(\ZZ)`, with coefficients in a ring.

    EXAMPLES::

        sage: QM = QuasiModularForms(1); QM
        Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
        sage: QM.gens()
        [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
         1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
         1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]

    It is possible to access the weight 2 Eisenstein series::

        sage: QM.weight_2_eisenstein_series()
        1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)

    Currently, the only supported base ring is the rational numbers::

        sage: QuasiModularForms(1, GF(5))
        Traceback (most recent call last):
        ...
        NotImplementedError: base ring other than Q are not yet supported for quasimodular forms ring
    """
    Element = QuasiModularFormsElement

    def __init__(self, group=1, base_ring=QQ, name='E2'):
        r"""
        INPUT:

        - ``group`` -- (default: `\SL_2(\ZZ)`) a congruence subgroup of
          `\SL_2(\ZZ)`, or a positive integer `N` (interpreted as
          `\Gamma_0(N)`)

        - ``base_ring`` -- a base ring (default: `\QQ`); should be
          `\QQ`, `\ZZ`, or the integers mod `p` for some prime `p`

        - ``name`` -- string (default: ``'E2'``); a variable name corresponding to
          the weight 2 Eisenstein series

        TESTS:

            sage: M = QuasiModularForms(1)
            sage: M.group()
            Modular Group SL(2,Z)
            sage: M.base_ring()
            Rational Field
            sage: QuasiModularForms(Integers(5))
            Traceback (most recent call last):
            ...
            ValueError: Group (=Ring of integers modulo 5) should be a congruence subgroup

        ::

            sage: TestSuite(QuasiModularForms(1)).run()
            sage: TestSuite(QuasiModularForms(Gamma0(3))).run()
            sage: TestSuite(QuasiModularForms(Gamma1(3))).run()
        """
        if not isinstance(name, str):
            raise TypeError("`name` must be a string")
        # check if the group is SL2(Z)
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        elif not isinstance(group, CongruenceSubgroupBase):
            raise ValueError("Group (=%s) should be a congruence subgroup" % group)

        # Check if the base ring is the rational field
        if base_ring != QQ:
            raise NotImplementedError("base ring other than Q are not yet supported for quasimodular forms ring")

        self.__group = group
        self.__modular_forms_subring = ModularFormsRing(group, base_ring)
        self.__polynomial_subring = self.__modular_forms_subring[name]
        cat = GradedAlgebras(base_ring).Commutative()
        Parent.__init__(self, base=base_ring, category=cat)

    def group(self):
        r"""
        Return the congruence subgroup attached to the given quasimodular forms
        ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.group()
            Modular Group SL(2,Z)
            sage: QM.group() is SL2Z
            True
            sage: QuasiModularForms(3).group()
            Congruence Subgroup Gamma0(3)
            sage: QuasiModularForms(Gamma1(5)).group()
            Congruence Subgroup Gamma1(5)
        """
        return self.__group

    def modular_forms_subring(self):
        r"""
        Return the subring of modular forms of this ring of quasimodular forms.

        EXAMPLES::

            sage: QuasiModularForms(1).modular_forms_subring()
            Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
            sage: QuasiModularForms(5).modular_forms_subring()
            Ring of Modular Forms for Congruence Subgroup Gamma0(5) over Rational Field
        """
        return self.__modular_forms_subring

    def modular_forms_of_weight(self, weight):
        r"""
        Return the space of modular forms on this group of the given weight.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.modular_forms_of_weight(12)
            Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
            sage: QM = QuasiModularForms(Gamma1(3))
            sage: QM.modular_forms_of_weight(4)
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma1(3) of weight 4 over Rational Field
        """
        return self.__modular_forms_subring.modular_forms_of_weight(weight)

    def quasimodular_forms_of_weight(self, weight):
        r"""
        Return the space of quasimodular forms on this group of the given weight.

        INPUT:

        - ``weight`` -- integer

        OUTPUT: a quasimodular forms space of the given weight

        EXAMPLES::

            sage: QuasiModularForms(1).quasimodular_forms_of_weight(4)
            Traceback (most recent call last):
            ...
            NotImplementedError: spaces of quasimodular forms of fixed weight not yet implemented
        """
        raise NotImplementedError("spaces of quasimodular forms of fixed weight not yet implemented")

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: QuasiModularForms(1)._repr_()
            'Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field'
        """
        return "Ring of Quasimodular Forms for %s over %s" % (self.group(), self.base_ring())

    def _coerce_map_from_(self, M):
        r"""
        Code to make QuasiModularForms work well with coercion framework.

        TESTS::

            sage: E2 = QuasiModularForms(1).0
            sage: M = ModularFormsRing(1)
            sage: E2 + M.0
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
            sage: M.0 + E2
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
            sage: 1 + E2
            2 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: E2 + 1
            2 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: f = ModularForms(1, 12).0
            sage: E2 + f
            1 - 23*q - 96*q^2 + 156*q^3 - 1640*q^4 + 4686*q^5 + O(q^6)
            sage: f + E2
            1 - 23*q - 96*q^2 + 156*q^3 - 1640*q^4 + 4686*q^5 + O(q^6)
        """
        if isinstance(M, (ModularFormsRing, ModularFormsSpace)):
            if M.group() == self.group() and self.has_coerce_map_from(M.base_ring()):
                return True
        if self.base_ring().has_coerce_map_from(M):
            return True
        return False

    def _element_constructor_(self, datum):
        r"""
        The call method of ``self``.

        INPUT:

        - ``datum`` -- list; GradedModularFormElement, ModularFormElement,
          Polynomial, base ring element

        OUTPUT: QuasiModularFormElement

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: M = QM.modular_forms_subring()
            sage: m12 = QM.modular_forms_of_weight(12)
            sage: QM([M.0, M.1])
            2 - 288*q - 2448*q^2 + 319104*q^3 + 3681936*q^4 + 21775680*q^5 + O(q^6)
            sage: QM([m12.0, m12.1])
            1 + 49627/691*q + 132611664/691*q^2 + 8380115796/691*q^3 - 13290096200/691*q^4 - 4248043226454/691*q^5 + O(q^6)
            sage: QM([])
            Traceback (most recent call last):
            ...
            ValueError: the given list should be non-empty
            sage: QM(M.0)
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: QM(m12.0)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            sage: y = polygen(QQ)
            sage: QM(y)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM(1 + y + y^2)
            3 - 72*q + 360*q^2 + 3168*q^3 + 9288*q^4 + 21456*q^5 + O(q^6)
            sage: QM(1)
            1
            sage: QM(1/2)
            1/2
            sage: QM('E2')
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from <class 'str'> to Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
            sage: P.<q> = PowerSeriesRing(QQ)
            sage: QM(1 - 24 * q - 72 * q^2 - 96 * q^3 + O(q^4))
            Traceback (most recent call last):
            ...
            NotImplementedError: conversion from q-expansion not yet implemented
        """
        if isinstance(datum, list):
            if not datum:
                raise ValueError("the given list should be non-empty")
            for idx, f in enumerate(datum):
                if not isinstance(f, (GradedModularFormElement, ModularFormElement)):
                    raise ValueError("one list element is not a modular form")
                datum[idx] = self.__modular_forms_subring(f)  # to ensure that every form is a GradedModularFormElement
            datum = self.__polynomial_subring(datum)
        elif isinstance(datum, (GradedModularFormElement, ModularFormElement)):
            datum = self.__modular_forms_subring(datum)  # GradedModularFormElement
            datum = self.__polynomial_subring(datum)
        elif isinstance(datum, Polynomial):
            datum = self.__polynomial_subring(datum.coefficients(sparse=False))
        elif isinstance(datum, PowerSeries_poly):
            raise NotImplementedError("conversion from q-expansion not yet implemented")
        else:
            datum = self.__polynomial_subring.coerce(datum)
        return self.element_class(self, datum)

    def weight_2_eisenstein_series(self):
        r"""
        Return the weight 2 Eisenstein series.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: E2 = QM.weight_2_eisenstein_series(); E2
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: E2.parent()
            Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
        """
        return self(self.__polynomial_subring.gen())

    def gens(self):
        r"""
        Return a list of generators of the quasimodular forms ring.

        Note that the generators of the modular forms subring are the one given
        by the method :meth:`sage.modular.modform.ring.ModularFormsRing.gen_forms`

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.gens()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
            sage: QM.modular_forms_subring().gen_forms()
            [1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
            sage: QM = QuasiModularForms(5)
            sage: QM.gens()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 6*q + 18*q^2 + 24*q^3 + 42*q^4 + 6*q^5 + O(q^6),
            1 + 240*q^5 + O(q^6),
            q + 10*q^3 + 28*q^4 + 35*q^5 + O(q^6)]

        An alias of this method is ``generators``::

            sage: QuasiModularForms(1).generators()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        """
        gen_list = [self.weight_2_eisenstein_series()]
        for f in self.__modular_forms_subring.gen_forms():
            gen_list.append(self(f))
        return gen_list

    generators = gens  # alias

    def ngens(self):
        r"""
        Return the number of generators of the given graded quasimodular forms
        ring.

        EXAMPLES::

            sage: QuasiModularForms(1).ngens()
            3
        """
        return len(self.gens())

    def gen(self, n):
        r"""
        Return the `n`-th generator of the quasimodular forms ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.0
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.1
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: QM.2
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
            sage: QM = QuasiModularForms(5)
            sage: QM.0
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.1
            1 + 6*q + 18*q^2 + 24*q^3 + 42*q^4 + 6*q^5 + O(q^6)
            sage: QM.2
            1 + 240*q^5 + O(q^6)
            sage: QM.3
            q + 10*q^3 + 28*q^4 + 35*q^5 + O(q^6)
            sage: QM.4
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.gens()[n]

    def zero(self):
        r"""
        Return the zero element of this ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.zero()
            0
            sage: QM.zero().is_zero()
            True
        """
        return self.element_class(self, self.__polynomial_subring.zero())

    def one(self):
        r"""
        Return the one element of this ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.one()
            1
            sage: QM.one().is_one()
            True
        """
        return self.element_class(self, self.__polynomial_subring.one())

    def some_elements(self):
        r"""
        Return a list of generators of ``self``.

        EXAMPLES::

            sage: QuasiModularForms(1).some_elements()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        """
        return self.gens()

    def polygen(self):
        r"""
        Return the generator of this quasimodular form space as a polynomial
        ring over the modular form subring.

        Note that this generator correspond to the weight-2 Eisenstein series.
        The default name of this generator is ``E2``.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.polygen()
            E2
            sage: QuasiModularForms(1, name='X').polygen()
            X
            sage: QM.polygen().parent()
            Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
        """
        return self.__polynomial_subring.gen()

    def polynomial_ring(self, names=None):
        r"""
        Return a multivariate polynomial ring of which the quasimodular forms
        ring is a quotient.

        In the case of the full modular group, this ring is `R[E_2, E_4, E_6]`
        where `E_2`, `E_4` and `E_6` have degrees 2, 4 and 6 respectively.

        INPUT:

        - ``names``-- string (default: ``None``); list or tuple of names
          (strings), or a comma separated string. Defines the names for the
          generators of the multivariate polynomial ring. The default names are
          of the following form:

          - ``E2`` denotes the weight 2 Eisenstein series;

          - ``Ek_i`` and ``Sk_i`` denote the `i`-th basis element of the weight
            `k` Eisenstein subspace and cuspidal subspace respectively;

          - If the level is one, the default names are ``E2``, ``E4`` and
            ``E6``;

          - In any other cases, we use the letters ``Fk``, ``Gk``, ``Hk``, ...,
            ``FFk``, ``FGk``, ... to denote any generator of weight `k`.

        OUTPUT: a multivariate polynomial ring in the variables ``names``

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: P = QM.polynomial_ring(); P
            Multivariate Polynomial Ring in E2, E4, E6 over Rational Field
            sage: P.inject_variables()
            Defining E2, E4, E6
            sage: E2.degree()
            2
            sage: E4.degree()
            4
            sage: E6.degree()
            6

        Example when the level is not one::

            sage: QM = QuasiModularForms(Gamma0(29))
            sage: P_29 = QM.polynomial_ring()
            sage: P_29
            Multivariate Polynomial Ring in E2, F2, S2_0, S2_1, E4_0, F4, G4, H4 over Rational Field
            sage: P_29.inject_variables()
            Defining E2, F2, S2_0, S2_1, E4_0, F4, G4, H4
            sage: F2.degree()
            2
            sage: E4_0.degree()
            4

        The name ``Sk_i`` stands for the `i`-th basis element of the cuspidal subspace of weight `k`::

            sage: F2 = QM.from_polynomial(S2_0)
            sage: F2.qexp(10)
            q - q^4 - q^5 - q^6 + 2*q^7 - 2*q^8 - 2*q^9 + O(q^10)
            sage: CuspForms(Gamma0(29), 2).0.qexp(10)
            q - q^4 - q^5 - q^6 + 2*q^7 - 2*q^8 - 2*q^9 + O(q^10)
            sage: F2 == CuspForms(Gamma0(29), 2).0
            True

        The name ``Ek_i`` stands for the `i`-th basis element of the Eisenstein subspace of weight `k`::

            sage: F4 = QM.from_polynomial(E4_0)
            sage: F4.qexp(30)
            1 + 240*q^29 + O(q^30)
            sage: EisensteinForms(Gamma0(29), 4).0.qexp(30)
            1 + 240*q^29 + O(q^30)
            sage: F4 == EisensteinForms(Gamma0(29), 4).0
            True

        One may also choose the name of the variables::

            sage: QM = QuasiModularForms(1)
            sage: QM.polynomial_ring(names="P, Q, R")
            Multivariate Polynomial Ring in P, Q, R over Rational Field
        """
        gens = self.__modular_forms_subring.gen_forms()
        weights = [f.weight() for f in gens]
        gens = iter(gens)
        if names is None:
            if self.group() == Gamma0(1):
                names = ["E2", "E4", "E6"]
            else:
                names = ["E2"]
                letters = "FGHIJK"
                for unique_weight in set(weights):
                    same_weights = [k for k in weights if k == unique_weight]
                    # create all the names of the form:
                    #     F, G, H, I, J, K, FF, FG, FH,..., FFF, FFG,...
                    # the letters E and S are reserved for basis elements of the
                    # Eisenstein subspaces and cuspidal subspaces respectively.
                    iter_names = (product(letters, repeat=r)
                                  for r in range(1, len(same_weights)//len(letters) + 2))
                    iter_names = chain(*iter_names)
                    for k in same_weights:
                        form = next(gens)
                        Mk = self.__modular_forms_subring.modular_forms_of_weight(k)
                        if form.is_eisenstein():
                            Ek_basis = Mk.eisenstein_subspace().basis()
                            # check if form is a basis element of the Eisenstein subspace of weight k
                            try:
                                n = Ek_basis.index(form)
                                name = f"E{str(k)}_{str(n)}"
                            except ValueError:
                                name = "".join(next(iter_names)) + str(k)
                        elif form.is_cuspidal():
                            Sk_basis = Mk.cuspidal_subspace().basis()
                            # check if form is a basis element of the cuspidal subspace of weight k
                            try:
                                n = Sk_basis.index(form)
                                name = f"S{str(k)}_{str(n)}"
                            except ValueError:
                                name = "".join(next(iter_names)) + str(k)
                        else:
                            name = "".join(next(iter_names)) + str(k)
                        names.append(name)
        weights.insert(0, 2)  # add the weight 2 Eisenstein series
        return PolynomialRing(self.base_ring(), len(weights), names,
                              order=TermOrder('wdeglex', weights))

    def from_polynomial(self, polynomial):
        r"""
        Convert the given polynomial `P(x,\ldots, y)` to the graded quasiform
        `P(g_0, \ldots, g_n)` where the `g_i` are the generators given
        by :meth:`~sage.modular.quasimodform.ring.QuasiModularForms.gens`.

        INPUT:

        - ``polynomial`` -- a multivariate polynomial

        OUTPUT: the graded quasimodular forms `P(g_0, \ldots, g_n)`

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: P.<x, y, z> = QQ[]
            sage: QM.from_polynomial(x)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.from_polynomial(x) == QM.0
            True
            sage: QM.from_polynomial(y) == QM.1
            True
            sage: QM.from_polynomial(z) == QM.2
            True
            sage: QM.from_polynomial(x^2 + y + x*z + 1)
            4 - 336*q - 2016*q^2 + 322368*q^3 + 3691392*q^4 + 21797280*q^5 + O(q^6)
            sage: QM = QuasiModularForms(Gamma0(2))
            sage: P = QM.polynomial_ring()
            sage: P.inject_variables()
            Defining E2, E2_0, E4_0
            sage: QM.from_polynomial(E2)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.from_polynomial(E2 + E4_0*E2_0) == QM.0 + QM.2*QM.1
            True

        Naturally, the number of variable must not exceed the number of generators::

            sage: P = PolynomialRing(QQ, 'F', 4)
            sage: P.inject_variables()
            Defining F0, F1, F2, F3
            sage: QM.from_polynomial(F0 + F1 + F2 + F3)
            Traceback (most recent call last):
            ...
            ValueError: the number of variables (4) of the given polynomial cannot exceed the number of generators (3) of the quasimodular forms ring


        TESTS::

            sage: QuasiModularForms(1).from_polynomial('x')
            Traceback (most recent call last):
            ...
            TypeError: the input must be a polynomial
        """
        if not isinstance(polynomial, (MPolynomial, Polynomial)):
            raise TypeError('the input must be a polynomial')
        poly_parent = polynomial.parent()
        nb_var = poly_parent.ngens()
        if nb_var > self.ngens():
            raise ValueError("the number of variables (%s) of the given polynomial cannot exceed the number of generators (%s) of the quasimodular forms ring" % (nb_var, self.ngens()))
        gens_dict = {poly_parent.gen(i): self.gen(i) for i in range(nb_var)}
        return self(polynomial.subs(gens_dict))

    def basis_of_weight(self, weight):
        r"""
        Return a basis of elements generating the subspace of the given
        weight.

        INPUT:

        - ``weight`` -- integer; the weight of the subspace

        OUTPUT: list of quasimodular forms of the given weight

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.basis_of_weight(12)
            [q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6),
             1 + 65520/691*q + 134250480/691*q^2 + 11606736960/691*q^3 + 274945048560/691*q^4 + 3199218815520/691*q^5 + O(q^6),
             1 - 288*q - 129168*q^2 - 1927296*q^3 + 65152656*q^4 + 1535768640*q^5 + O(q^6),
             1 + 432*q + 39312*q^2 - 1711296*q^3 - 14159664*q^4 + 317412000*q^5 + O(q^6),
             1 - 576*q + 21168*q^2 + 308736*q^3 - 15034608*q^4 - 39208320*q^5 + O(q^6),
             1 + 144*q - 17712*q^2 + 524736*q^3 - 2279088*q^4 - 79760160*q^5 + O(q^6),
             1 - 144*q + 8208*q^2 - 225216*q^3 + 2634192*q^4 + 1488672*q^5 + O(q^6)]
            sage: QM = QuasiModularForms(Gamma1(3))
            sage: QM.basis_of_weight(3)
            [1 + 54*q^2 + 72*q^3 + 432*q^5 + O(q^6),
             q + 3*q^2 + 9*q^3 + 13*q^4 + 24*q^5 + O(q^6)]
            sage: QM.basis_of_weight(5)
            [1 - 90*q^2 - 240*q^3 - 3744*q^5 + O(q^6),
             q + 15*q^2 + 81*q^3 + 241*q^4 + 624*q^5 + O(q^6),
             1 - 24*q - 18*q^2 - 1320*q^3 - 5784*q^4 - 10080*q^5 + O(q^6),
             q - 21*q^2 - 135*q^3 - 515*q^4 - 1392*q^5 + O(q^6)]
        """
        basis = []
        E2 = self.weight_2_eisenstein_series()
        M = self.__modular_forms_subring
        E2_pow = self.one()
        for j in range(weight // 2):
            basis.extend(f * E2_pow
                         for f in M.modular_forms_of_weight(weight - 2*j).basis())
            E2_pow *= E2
        if not weight % 2:
            basis.append(E2_pow)
        return basis
