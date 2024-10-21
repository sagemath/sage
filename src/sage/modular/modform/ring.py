# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Graded rings of modular forms

This module contains functions to find generators for the graded ring of
modular forms of given level.

AUTHORS:

- William Stein (2007-08-24): first version
- David Ayotte (2021-06): implemented category and Parent/Element frameworks
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein
#                     2021 David Ayotte
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from random import shuffle

from sage.categories.graded_algebras import GradedAlgebras
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.misc.superseded import deprecated_function_alias
from sage.misc.verbose import verbose
from sage.modular.arithgroup.all import Gamma0, CongruenceSubgroupBase
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.power_series_poly import PowerSeries_poly
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_method, richcmp

from .constructor import ModularForms
from .element import ModularFormElement, GradedModularFormElement
from .space import ModularFormsSpace


def _span_of_forms_in_weight(forms, weight, prec, stop_dim=None, use_random=False):
    r"""
    Utility function. Given a nonempty list of pairs ``(k,f)``, where `k` is an
    integer and `f` is a power series, and a weight l, return all weight l
    forms obtained by multiplying together the given forms.

    INPUT:

    - ``forms`` -- list of pairs `(k, f)` with k an integer and f a power
      series (all over the same base ring)
    - ``weight`` -- integer
    - ``prec`` -- integer (less than or equal to the precision of all the
      forms in ``forms``); precision to use in power series computations
    - ``stop_dim`` -- integer; stop as soon as we have enough forms to span
      a submodule of this rank (a saturated one if the base ring is `\ZZ`).
      Ignored if ``use_random`` is ``False``.
    - ``use_random`` -- which algorithm to use. If ``True``, tries random products
      of the generators of the appropriate weight until a large enough
      submodule is found (determined by ``stop_dim``). If ``False``, just tries
      everything.

    Note that if the given forms do generate the whole space, then
    ``use_random=True`` will often be quicker (particularly if the weight is
    large); but if the forms don't generate, the randomized algorithm is no
    help and will actually be substantially slower, because it needs to do
    repeated echelon form calls to check if vectors are in a submodule, while
    the non-randomized algorithm just echelonizes one enormous matrix at the
    end.

    EXAMPLES::

        sage: import sage.modular.modform.ring as f
        sage: forms = [(4, 240*eisenstein_series_qexp(4,5)), (6,504*eisenstein_series_qexp(6,5))]
        sage: f._span_of_forms_in_weight(forms, 12, prec=5)
        Vector space of degree 5 and dimension 2 over Rational Field
        Basis matrix:
        [        1         0    196560  16773120 398034000]
        [        0         1       -24       252     -1472]
        sage: f._span_of_forms_in_weight(forms, 24, prec=5)
        Vector space of degree 5 and dimension 3 over Rational Field
        Basis matrix:
        [          1           0           0    52416000 39007332000]
        [          0           1           0      195660    12080128]
        [          0           0           1         -48        1080]
        sage: ModularForms(1, 24).q_echelon_basis(prec=5)
        [
        1 + 52416000*q^3 + 39007332000*q^4 + O(q^5),
        q + 195660*q^3 + 12080128*q^4 + O(q^5),
        q^2 - 48*q^3 + 1080*q^4 + O(q^5)
        ]

    Test the alternative randomized algorithm::

        sage: f._span_of_forms_in_weight(forms, 24, prec=5, use_random=True, stop_dim=3)
        Vector space of degree 5 and dimension 3 over Rational Field
        Basis matrix:
        [          1           0           0    52416000 39007332000]
        [          0           1           0      195660    12080128]
        [          0           0           1         -48        1080]
    """
    t = verbose('multiplying forms up to weight %s' % weight)
    # Algorithm: run through the monomials of the appropriate weight, and build
    # up the vector space they span.

    n = len(forms)
    R = forms[0][1].base_ring()
    V = R ** prec
    W = V.zero_submodule()
    shortforms = [f[1].truncate_powerseries(prec) for f in forms]

    # List of weights
    from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
    wts = list(WeightedIntegerVectors(weight, [f[0] for f in forms]))
    t = verbose("calculated weight list", t)
    N = len(wts)

    if use_random:
        if stop_dim is None:
            raise ValueError("stop_dim must be provided if use_random is True")
        shuffle(wts)

        for c in range(N):
            w = V(prod(shortforms[i]**wts[c][i] for i in range(n)).padded_list(prec))
            if w in W:
                continue
            W = V.span(list(W.gens()) + [w])
            if stop_dim and W.rank() == stop_dim:
                if R != ZZ or W.index_in_saturation() == 1:
                    verbose("Succeeded after %s of %s" % (c, N), t)
                    return W
        verbose("Nothing worked", t)
        return W
    else:
        G = [V(prod(forms[i][1]**c[i] for i in range(n)).padded_list(prec)) for c in wts]
        t = verbose('found %s candidates' % N, t)
        W = V.span(G)
        verbose('span has dimension %s' % W.rank(), t)
        return W


@richcmp_method
class ModularFormsRing(Parent):
    r"""
    The ring of modular forms (of weights 0 or at least 2) for a congruence
    subgroup of `\SL_2(\ZZ)`, with coefficients in a specified base ring.

    EXAMPLES::

        sage: ModularFormsRing(Gamma1(13))
        Ring of Modular Forms for Congruence Subgroup Gamma1(13) over Rational Field
        sage: m = ModularFormsRing(4); m
        Ring of Modular Forms for Congruence Subgroup Gamma0(4) over Rational Field
        sage: m.modular_forms_of_weight(2)
        Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(4) of weight 2 over Rational Field
        sage: m.modular_forms_of_weight(10)
        Modular Forms space of dimension 6 for Congruence Subgroup Gamma0(4) of weight 10 over Rational Field
        sage: m == loads(dumps(m))
        True
        sage: m.generators()
        [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10)),
         (2, q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10))]
        sage: m.q_expansion_basis(2,10)
        [1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10),
         q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10)]
        sage: m.q_expansion_basis(3,10)
        []
        sage: m.q_expansion_basis(10,10)
        [1 + 10560*q^6 + 3960*q^8 + O(q^10),
         q - 8056*q^7 - 30855*q^9 + O(q^10),
         q^2 - 796*q^6 - 8192*q^8 + O(q^10),
         q^3 + 66*q^7 + 832*q^9 + O(q^10),
         q^4 + 40*q^6 + 528*q^8 + O(q^10),
         q^5 + 20*q^7 + 190*q^9 + O(q^10)]

    Elements of modular forms ring can be initiated via multivariate polynomials (see :meth:`from_polynomial`)::

        sage: M = ModularFormsRing(1)
        sage: M.ngens()
        2
        sage: E4, E6 = polygens(QQ, 'E4, E6')
        sage: M(E4)
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
        sage: M(E6)
        1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
        sage: M((E4^3 - E6^2)/1728)
        q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
    """

    Element = GradedModularFormElement

    def __init__(self, group, base_ring=QQ):
        r"""
        INPUT:

        - ``group`` -- a congruence subgroup of `\SL_2(\ZZ)`, or a
          positive integer `N` (interpreted as `\Gamma_0(N)`)

        - ``base_ring`` -- ring (default: `\QQ`); a base ring, which
          should be `\QQ`, `\ZZ`, or the integers mod `p` for some prime
          `p`

        TESTS:

        Check that :issue:`15037` is fixed::

            sage: ModularFormsRing(3.4)
            Traceback (most recent call last):
            ...
            ValueError: group (=3.40000000000000) should be a congruence subgroup
            sage: ModularFormsRing(Gamma0(2), base_ring=PolynomialRing(ZZ, 'x'))
            Traceback (most recent call last):
            ...
            ValueError: base ring (=Univariate Polynomial Ring in x over Integer Ring) should be QQ, ZZ or a finite prime field

        ::

            sage: TestSuite(ModularFormsRing(1)).run()
            sage: TestSuite(ModularFormsRing(Gamma0(6))).run()
            sage: TestSuite(ModularFormsRing(Gamma1(4))).run()

        .. TODO::

            - Add graded modular forms over non-trivial Dirichlet character;
            - makes gen_forms returns modular forms over base rings other than `QQ`;
            - implement binary operations between two forms with different groups.
        """
        if isinstance(group, (int, Integer)):
            group = Gamma0(group)
        elif not isinstance(group, CongruenceSubgroupBase):
            raise ValueError("group (=%s) should be a congruence subgroup" % group)

        if base_ring != ZZ and not base_ring.is_field() and not base_ring.is_finite():
            raise ValueError("base ring (=%s) should be QQ, ZZ or a finite prime field" % base_ring)

        self.__group = group
        self.__cached_maxweight = ZZ(-1)
        self.__cached_gens = []
        self.__cached_cusp_maxweight = ZZ(-1)
        self.__cached_cusp_gens = []
        cat = GradedAlgebras(base_ring).Commutative()
        Parent.__init__(self, base=base_ring, category=cat)

    def change_ring(self, base_ring):
        r"""
        Return a ring of modular forms over a new base ring of the same
        congruence subgroup.

        INPUT:

        - ``base_ring`` -- a base ring, which should be `\QQ`, `\ZZ`, or
          the integers mod `p` for some prime `p`

        EXAMPLES::

            sage: M = ModularFormsRing(11); M
            Ring of Modular Forms for Congruence Subgroup Gamma0(11) over Rational Field
            sage: M.change_ring(Zmod(7))
            Ring of Modular Forms for Congruence Subgroup Gamma0(11) over Ring of integers modulo 7
            sage: M.change_ring(ZZ)
            Ring of Modular Forms for Congruence Subgroup Gamma0(11) over Integer Ring
        """
        return ModularFormsRing(self.group(), base_ring=base_ring)

    def some_elements(self):
        r"""
        Return some elements of this ring.

        EXAMPLES::

            sage: ModularFormsRing(1).some_elements()
            [1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
             1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        """
        return [self(f) for f in self.gen_forms()]

    def group(self):
        r"""
        Return the congruence subgroup of this ring of modular forms.

        EXAMPLES::

            sage: R = ModularFormsRing(Gamma1(13))
            sage: R.group() is Gamma1(13)
            True
        """
        return self.__group

    def gen(self, i):
        r"""
        Return the `i`-th generator of this ring.

        INPUT:

        - ``i`` -- integer

        OUTPUT: an instance of :class:`~sage.modular.modform.GradedModularFormElement`

        EXAMPLES::

            sage: M = ModularFormsRing(1)
            sage: E4 = M.0; E4 # indirect doctest
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: E6 = M.1; E6 # indirect doctest
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
        """
        if self.base_ring() is not QQ:
            raise NotImplementedError("the base ring of the given ring of modular form should be QQ")
        return self(self.gen_forms()[i])

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        EXAMPLES::

            sage: ModularFormsRing(1).ngens()
            2
            sage: ModularFormsRing(Gamma0(2)).ngens()
            2
            sage: ModularFormsRing(Gamma1(13)).ngens() # long time
            33

        .. WARNING::

            Computing the number of generators of a graded ring of modular form for a certain
            congruence subgroup can be very long.
        """
        return len(self.gen_forms())

    def polynomial_ring(self, names, gens=None):
        r"""
        Return a polynomial ring of which this ring of modular forms is
        a quotient.

        INPUT:

        - ``names`` -- a list or tuple of names (strings), or a comma
          separated string; consists in the names of the polynomial
          ring variables
        - ``gens`` -- list of modular forms generating this ring
          (default: ``None``); if ``gens`` is ``None`` then the list of
          generators returned by the method
          :meth:`~sage.modular.modform.find_generator.ModularFormsRing.gen_forms`
          is used instead. Note that we do not check if the list is
          indeed a generating set.

        OUTPUT: a multivariate polynomial ring in the variable
        ``names``. Each variable of the polynomial ring correspond to a
        generator given in the list ``gens`` (following the ordering of
        the list).

        EXAMPLES::

            sage: M = ModularFormsRing(1)
            sage: gens = M.gen_forms()
            sage: M.polynomial_ring('E4, E6', gens)
            Multivariate Polynomial Ring in E4, E6 over Rational Field
            sage: M = ModularFormsRing(Gamma0(8))
            sage: gens = M.gen_forms()
            sage: M.polynomial_ring('g', gens)
            Multivariate Polynomial Ring in g0, g1, g2 over Rational Field

        The degrees of the variables are the weights of the
        corresponding forms::

            sage: M = ModularFormsRing(1)
            sage: P.<E4, E6> = M.polynomial_ring()
            sage: E4.degree()
            4
            sage: E6.degree()
            6
            sage: (E4*E6).degree()
            10
        """
        if gens is None:
            gens = self.gen_forms()
        degs = [f.weight() for f in gens]
        return PolynomialRing(self.base_ring(), len(gens), names,
                              order=TermOrder('wdeglex', degs))

    def _generators_variables_dictionnary(self, poly_parent, gens):
        r"""
        Return a dictionary giving an association between polynomial
        ring generators and generators of modular forms ring.

        INPUT:

        - ``poly_parent`` -- a polynomial ring
        - ``gen`` -- list of generators of the modular forms ring

        TESTS::

            sage: M = ModularFormsRing(Gamma0(6))
            sage: P = QQ['x, y, z']
            sage: M._generators_variables_dictionnary(P, M.gen_forms())
            {z: q^2 - 2*q^3 + 3*q^4 + O(q^6),
             y: q + 5*q^3 - 2*q^4 + 6*q^5 + O(q^6),
             x: 1 + 24*q^3 + O(q^6)}
        """
        if poly_parent.base_ring() != self.base_ring():
            raise ValueError('the base ring of `poly_parent` must be the same as the base ring of the modular forms ring')
        nb_var = poly_parent.ngens()
        nb_gens = self.ngens()
        if nb_var != nb_gens:
            raise ValueError('the number of variables (%s) must be equal to'
                             ' the number of generators of the modular forms'
                             ' ring (%s)' % (nb_var, self.ngens()))
        return {poly_parent.gen(i): self(gens[i]) for i in range(0, nb_var)}

    def from_polynomial(self, polynomial, gens=None):
        r"""
        Return a graded modular form constructed by evaluating a given
        multivariate polynomial at a set of generators.

        INPUT:

        - ``polynomial`` -- a multivariate polynomial. The variables
          names of the polynomial should be different from ``'q'``. The
          number of variable of this polynomial should equal the number
          of given generators.
        - ``gens`` -- list of modular forms generating this ring
          (default: ``None``); if ``gens`` is ``None`` then the list of
          generators returned by the method
          :meth:`~sage.modular.modform.find_generator.ModularFormsRing.gen_forms`
          is used instead. Note that we do not check if the list is
          indeed a generating set.

        OUTPUT: a ``GradedModularFormElement`` given by the polynomial
        relation ``polynomial``

        EXAMPLES::

            sage: M = ModularFormsRing(1)
            sage: x,y = polygens(QQ, 'x,y')
            sage: M.from_polynomial(x^2+y^3)
            2 - 1032*q + 774072*q^2 - 77047584*q^3 - 11466304584*q^4 - 498052467504*q^5 + O(q^6)
            sage: M = ModularFormsRing(Gamma0(6))
            sage: M.ngens()
            3
            sage: x,y,z = polygens(QQ, 'x,y,z')
            sage: M.from_polynomial(x+y+z)
            1 + q + q^2 + 27*q^3 + q^4 + 6*q^5 + O(q^6)
            sage: M.0 + M.1 + M.2
            1 + q + q^2 + 27*q^3 + q^4 + 6*q^5 + O(q^6)
            sage: P = x.parent()
            sage: M.from_polynomial(P(1/2))
            1/2

        Note that the number of variables must be equal to the number of
        generators::

            sage: x, y = polygens(QQ, 'x, y')
            sage: M(x + y)
            Traceback (most recent call last):
            ...
            ValueError: the number of variables (2) must be equal to the number of generators of the modular forms ring (3)

        TESTS::

            sage: x,y = polygens(GF(7), 'x, y')
            sage: ModularFormsRing(1, GF(7))(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: conversion from polynomial is not implemented if the base ring is not Q

        ..TODO::

            * add conversion for symbolic expressions?
        """
        if not self.base_ring() == QQ: # this comes from the method gens_form
            raise NotImplementedError("conversion from polynomial is not implemented if the base ring is not Q")
        if not isinstance(polynomial, MPolynomial):
            raise TypeError('`polynomial` must be a multivariate polynomial')
        if gens is None:
            gens = self.gen_forms()
        dict = self._generators_variables_dictionnary(polynomial.parent(), gens)
        if polynomial.is_constant():
            return self(polynomial.constant_coefficient())
        return polynomial.substitute(dict)

    def _element_constructor_(self, forms_datum):
        r"""
        Return the graded modular form corresponding to the given data.

        INPUT:

        - ``forms_datum`` -- dictionary, list, ModularFormElement,
          GradedModularFormElement, RingElement, or Multivariate polynomial; try
          to coerce ``forms_datum`` into ``self``

        TESTS::

            sage: M = ModularFormsRing(1)
            sage: E4 = ModularForms(1,4).0; E6 = ModularForms(1,6).0
            sage: M([E4, E6])
            2 - 264*q - 14472*q^2 - 116256*q^3 - 515208*q^4 - 1545264*q^5 + O(q^6)
            sage: M([E4, E4])
            2 + 480*q + 4320*q^2 + 13440*q^3 + 35040*q^4 + 60480*q^5 + O(q^6)
            sage: M({4:E4, 6:E6})
            2 - 264*q - 14472*q^2 - 116256*q^3 - 515208*q^4 - 1545264*q^5 + O(q^6)
            sage: M(E4)
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: x,y = polygens(QQ, 'x,y')
            sage: M(x^2+x^3 + x*y + y^6)
            4 - 2088*q + 3816216*q^2 - 2296935072*q^3 + 720715388184*q^4 - 77528994304752*q^5 + O(q^6)
            sage: f = ModularForms(3, 10).0
            sage: M(f)
            Traceback (most recent call last):
            ...
            ValueError: the group (Congruence Subgroup Gamma0(3)) and/or the base ring (Rational Field) of the given modular form is not consistant with the base space: Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
            sage: M = ModularFormsRing(1, base_ring=ZZ)
            sage: M(ModularForms(1,4).0)
            Traceback (most recent call last):
            ...
            ValueError: the group (Modular Group SL(2,Z)) and/or the base ring (Rational Field) of the given modular form is not consistant with the base space: Ring of Modular Forms for Modular Group SL(2,Z) over Integer Ring
            sage: M('x')
            Traceback (most recent call last):
            ...
            TypeError: the defining data structure should be a single modular form, a ring element, a list of modular forms, a multivariate polynomial or a dictionary
            sage: P.<t> = PowerSeriesRing(QQ)
            sage: e = 1 + 240*t + 2160*t^2 + 6720*t^3 + 17520*t^4 + 30240*t^5 + O(t^6)
            sage: ModularFormsRing(1)(e)
            Traceback (most recent call last):
            ...
            NotImplementedError: conversion from q-expansion not yet implemented
        """
        if isinstance(forms_datum, (dict, list)):
            forms_dictionary = forms_datum
        elif isinstance(forms_datum, self.element_class):
            forms_dictionary = forms_datum._forms_dictionary
        elif isinstance(forms_datum, ModularFormElement):
            if self.group().is_subgroup(forms_datum.group()) and self.base_ring().has_coerce_map_from(forms_datum.base_ring()):
                forms_dictionary = {forms_datum.weight(): forms_datum}
            else:
                raise ValueError('the group (%s) and/or the base ring (%s) of the given modular form is not consistant with the base space: %s' % (forms_datum.group(), forms_datum.base_ring(), self))
        elif forms_datum in self.base_ring():
            forms_dictionary = {0: forms_datum}
        elif isinstance(forms_datum, MPolynomial):
            return self.from_polynomial(forms_datum)
        elif isinstance(forms_datum, PowerSeries_poly):
            raise NotImplementedError("conversion from q-expansion not yet implemented")
        else:
            raise TypeError('the defining data structure should be a single modular form, a ring element, a list of modular forms, a multivariate polynomial or a dictionary')
        return self.element_class(self, forms_dictionary)

    def zero(self):
        r"""
        Return the zero element of this ring.

        EXAMPLES::

            sage: M = ModularFormsRing(1)
            sage: zer = M.zero(); zer
            0
            sage: zer.is_zero()
            True
            sage: E4 = ModularForms(1,4).0; E4
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: E4 + zer
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: zer * E4
            0
            sage: E4 * zer
            0
        """
        return self.element_class(self, {})

    def one(self):
        r"""
        Return the one element of this ring.

        EXAMPLES::

            sage: M = ModularFormsRing(1)
            sage: u = M.one(); u
            1
            sage: u.is_one()
            True
            sage: u + u
            2
            sage: E4 = ModularForms(1,4).0; E4
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: E4 * u
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
        """
        return self.element_class(self, {0: self.base_ring().one()})

    def _coerce_map_from_(self, M):
        r"""
        Return ``True`` if there is a coercion map from ``M`` to this
        ring.

        TESTS::

            sage: M = ModularFormsRing(1)
            sage: E6 = ModularForms(1,6).0
            sage: D = ModularForms(1,12).1
            sage: M(D) + E6
            2 - 282744/691*q + 122757768/691*q^2 + 11521760544/691*q^3 + 274576933512/691*q^4 + 3198130142256/691*q^5 + O(q^6)
            sage: D + M(E6)
            2 - 282744/691*q + 122757768/691*q^2 + 11521760544/691*q^3 + 274576933512/691*q^4 + 3198130142256/691*q^5 + O(q^6)
            sage: 14 + M(D)
            15 + 65520/691*q + 134250480/691*q^2 + 11606736960/691*q^3 + 274945048560/691*q^4 + 3199218815520/691*q^5 + O(q^6)
            sage: M(D) + 53
            54 + 65520/691*q + 134250480/691*q^2 + 11606736960/691*q^3 + 274945048560/691*q^4 + 3199218815520/691*q^5 + O(q^6)
        """
        if isinstance(M, ModularFormsSpace):
            if M.group() == self.group() and self.has_coerce_map_from(M.base_ring()):
                return True
        if self.base_ring().has_coerce_map_from(M):
            return True
        return False

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        Rings are equal if and only if their groups and base rings are.

        EXAMPLES::

            sage: ModularFormsRing(3) == 3
            False
            sage: ModularFormsRing(Gamma0(3)) == ModularFormsRing(Gamma0(7))
            False
            sage: ModularFormsRing(Gamma0(3)) == ModularFormsRing(Gamma0(3))
            True
        """
        if not isinstance(other, ModularFormsRing):
            return NotImplemented

        return richcmp((self.group(), self.base_ring()),
                       (other.group(), other.base_ring()), op)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ModularFormsRing(Gamma0(13))._repr_()
            'Ring of Modular Forms for Congruence Subgroup Gamma0(13) over Rational Field'
            sage: ModularFormsRing(Gamma1(13), base_ring=ZZ)._repr_()
            'Ring of Modular Forms for Congruence Subgroup Gamma1(13) over Integer Ring'
        """
        return "Ring of Modular Forms for %s over %s" % (self.group(), self.base_ring())

    def modular_forms_of_weight(self, weight):
        """
        Return the space of modular forms of the given weight and the
        same congruence subgroup.

        EXAMPLES::

            sage: R = ModularFormsRing(13)
            sage: R.modular_forms_of_weight(10)
            Modular Forms space of dimension 11 for Congruence Subgroup Gamma0(13) of weight 10 over Rational Field
            sage: ModularFormsRing(Gamma1(13)).modular_forms_of_weight(3)
            Modular Forms space of dimension 20 for Congruence Subgroup Gamma1(13) of weight 3 over Rational Field
        """
        return ModularForms(self.group(), weight)

    def generators(self, maxweight=8, prec=10, start_gens=[], start_weight=2):
        r"""
        Return a list of generator of this ring as a list of pairs
        `(k, f)` where `k` is an integer and `f` is a univariate power
        series in `q` corresponding to the `q`-expansion of a modular
        form of weight `k`.

        More precisely, if `R` is the base ring of self, then this
        function calculates a set of modular forms which generate the
        `R`-algebra of all modular forms of weight up to ``maxweight``
        with coefficients in `R`.

        INPUT:

        - ``maxweight`` -- integer (default: 8); check up to this weight
          for generators

        - ``prec`` -- integer (default: 10); return `q`-expansions to
          this precision

        - ``start_gens`` -- list (default: ``[]``); list of pairs
          `(k, f)`, or triples `(k, f, F)`, where:

          - `k` is an integer,
          - `f` is the `q`-expansion of a modular form of weight `k`,
            as a power series over the base ring of self,
          - `F` (if provided) is a modular form object corresponding to F.

          If this list is nonempty, we find a minimal generating set containing
          these forms. If `F` is not supplied, then `f` needs to have
          sufficiently large precision (an error will be raised if this is not
          the case); otherwise, more terms will be calculated from the modular
          form object `F`.

        - ``start_weight`` -- integer (default: 2); calculate the graded
          subalgebra of forms of weight at least ``start_weight``

        OUTPUT:

        a list of pairs (k, f), where f is the `q`-expansion to precision
        ``prec`` of a modular form of weight k.

        .. SEEALSO::

            :meth:`gen_forms`, which does exactly the same thing, but returns
            Sage modular form objects rather than bare power series, and keeps
            track of a lifting to characteristic 0 when the base ring is a
            finite field.

        .. NOTE::

            If called with the default values of ``start_gens`` (an empty list)
            and ``start_weight`` (2), the values will be cached for re-use on
            subsequent calls to this function. (This cache is shared with
            :meth:`gen_forms`). If called with non-default values for these
            parameters, caching will be disabled.

        EXAMPLES::

            sage: ModularFormsRing(SL2Z).generators()
            [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + 82560*q^7 + 140400*q^8 + 181680*q^9 + O(q^10)),
             (6, 1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 - 8471232*q^7 - 17047800*q^8 - 29883672*q^9 + O(q^10))]
            sage: s = ModularFormsRing(SL2Z).generators(maxweight=5, prec=3); s
            [(4, 1 + 240*q + 2160*q^2 + O(q^3))]
            sage: s[0][1].parent()
            Power Series Ring in q over Rational Field

            sage: ModularFormsRing(1).generators(prec=4)
            [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + O(q^4)),
             (6, 1 - 504*q - 16632*q^2 - 122976*q^3 + O(q^4))]
            sage: ModularFormsRing(2).generators(prec=12)
            [(2, 1 + 24*q + 24*q^2 + 96*q^3 + 24*q^4 + 144*q^5 + 96*q^6 + 192*q^7 + 24*q^8 + 312*q^9 + 144*q^10 + 288*q^11 + O(q^12)),
             (4, 1 + 240*q^2 + 2160*q^4 + 6720*q^6 + 17520*q^8 + 30240*q^10 + O(q^12))]
            sage: ModularFormsRing(4).generators(maxweight=2, prec=20)
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + 144*q^10 + 96*q^12 + 192*q^14 + 24*q^16 + 312*q^18 + O(q^20)),
             (2, q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + 12*q^11 + 14*q^13 + 24*q^15 + 18*q^17 + 20*q^19 + O(q^20))]

        Here we see that for ``\Gamma_0(11)`` taking a basis of forms in weights 2
        and 4 is enough to generate everything up to weight 12 (and probably
        everything else).::

            sage: v = ModularFormsRing(11).generators(maxweight=12)
            sage: len(v)
            3
            sage: [k for k, _ in v]
            [2, 2, 4]
            sage: from sage.modular.dims import dimension_modular_forms
            sage: dimension_modular_forms(11,2)
            2
            sage: dimension_modular_forms(11,4)
            4

        For congruence subgroups not containing -1, we miss out some forms since we
        can't calculate weight 1 forms at present, but we can still find generators
        for the ring of forms of weight `\ge 2`::

            sage: ModularFormsRing(Gamma1(4)).generators(prec=10, maxweight=10)
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10)),
             (2, q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10)),
             (3, 1 + 12*q^2 + 64*q^3 + 60*q^4 + 160*q^6 + 384*q^7 + 252*q^8 + O(q^10)),
             (3, q + 4*q^2 + 8*q^3 + 16*q^4 + 26*q^5 + 32*q^6 + 48*q^7 + 64*q^8 + 73*q^9 + O(q^10))]

        Using different base rings will change the generators::

            sage: ModularFormsRing(Gamma0(13)).generators(maxweight=12, prec=4)
            [(2, 1 + 2*q + 6*q^2 + 8*q^3 + O(q^4)),
             (4, 1 + O(q^4)), (4, q + O(q^4)),
             (4, q^2 + O(q^4)), (4, q^3 + O(q^4)),
             (6, 1 + O(q^4)),
             (6, q + O(q^4))]
            sage: ModularFormsRing(Gamma0(13),base_ring=ZZ).generators(maxweight=12, prec=4)
            [(2, 1 + 2*q + 6*q^2 + 8*q^3 + O(q^4)),
             (4, q + 4*q^2 + 10*q^3 + O(q^4)),
             (4, 2*q^2 + 5*q^3 + O(q^4)),
             (4, q^2 + O(q^4)),
             (4, -2*q^3 + O(q^4)),
             (6, O(q^4)),
             (6, O(q^4)),
             (12, O(q^4))]
            sage: [k for k,f in ModularFormsRing(1, QQ).generators(maxweight=12)]
            [4, 6]
            sage: [k for k,f in ModularFormsRing(1, ZZ).generators(maxweight=12)]
            [4, 6, 12]
            sage: [k for k,f in ModularFormsRing(1, Zmod(5)).generators(maxweight=12)]
            [4, 6]
            sage: [k for k,f in ModularFormsRing(1, Zmod(2)).generators(maxweight=12)]
            [4, 6, 12]

        An example where ``start_gens`` are specified::

            sage: M = ModularForms(11, 2); f = (M.0 + M.1).qexp(8)
            sage: ModularFormsRing(11).generators(start_gens = [(2, f)])
            Traceback (most recent call last):
            ...
            ValueError: Requested precision cannot be higher than precision of approximate starting generators!
            sage: f = (M.0 + M.1).qexp(10); f
            1 + 17/5*q + 26/5*q^2 + 43/5*q^3 + 94/5*q^4 + 77/5*q^5 + 154/5*q^6 + 86/5*q^7 + 36*q^8 + 146/5*q^9 + O(q^10)
            sage: ModularFormsRing(11).generators(start_gens = [(2, f)])
            [(2, 1 + 17/5*q + 26/5*q^2 + 43/5*q^3 + 94/5*q^4 + 77/5*q^5 + 154/5*q^6 + 86/5*q^7 + 36*q^8 + 146/5*q^9 + O(q^10)),
             (2, 1 + 12*q^2 + 12*q^3 + 12*q^4 + 12*q^5 + 24*q^6 + 24*q^7 + 36*q^8 + 36*q^9 + O(q^10)),
             (4, 1 + O(q^10))]
        """
        sgs = []
        for x in start_gens:
            if len(x) == 2:
                if x[1].prec() < prec:
                    raise ValueError("Requested precision cannot be higher"
                                     " than precision of approximate starting generators!")
                sgs.append((x[0], x[1], None))
            else:
                sgs.append(x)

        G = self._find_generators(maxweight, tuple(sgs), start_weight)

        ret = []
        # Returned generators may be a funny mixture of precisions if start_gens has been used.
        for k, f, F in G:
            if f.prec() < prec:
                f = F.qexp(prec).change_ring(self.base_ring())
            else:
                f = f.truncate_powerseries(prec)
            ret.append((k, f))

        return ret

    def gen_forms(self, maxweight=8, start_gens=[], start_weight=2):
        r"""
        Return a list of modular forms generating this ring (as an algebra
        over the appropriate base ring).

        This method differs from :meth:`generators` only in that it returns
        graded modular form objects, rather than bare `q`-expansions.

        INPUT:

        - ``maxweight`` -- integer (default: 8); calculate forms
          generating all forms up to this weight

        - ``start_gens`` -- list (default: ``[]``); a list of
          modular forms. If this list is nonempty, we find a minimal
          generating set containing these forms.

        - ``start_weight`` -- integer (default: 2); calculate the graded
          subalgebra of forms of weight at least ``start_weight``

        .. NOTE::

            If called with the default values of ``start_gens`` (an empty list)
            and ``start_weight`` (2), the values will be cached for re-use on
            subsequent calls to this function. (This cache is shared with
            :meth:`generators`). If called with non-default values for these
            parameters, caching will be disabled.

        EXAMPLES::

            sage: A = ModularFormsRing(Gamma0(11), Zmod(5)).gen_forms(); A
            [1 + 12*q^2 + 12*q^3 + 12*q^4 + 12*q^5 + O(q^6),
             q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6),
             q - 9*q^4 - 10*q^5 + O(q^6)]
            sage: A[0].parent()
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
        """
        sgs = tuple((F.weight(), None, F) for F in start_gens)
        G = self._find_generators(maxweight, sgs, start_weight)
        return [F for k, f, F in G]

    gens = gen_forms

    def _find_generators(self, maxweight, start_gens, start_weight):
        r"""
        Returns a list of triples `(k, f, F)` where `F` is a modular
        form of weight `k` and `f` is its `q`-expansion coerced into the
        base ring of self.

        For internal use. This function is called by :meth:`generators`
        and :meth:`gen_forms`.

        INPUT:

        - ``maxweight`` -- maximum weight to try
        - ``start_weight`` -- minimum weight to try
        - ``start_gens`` -- a sequence of tuples of the form `(k, f, F)`, where
          `F` is a modular form of weight `k` and `f` is its `q`-expansion
          coerced into ``self.base_ring()``. Either (but not both) of `f` and `F`
          may be ``None``.

        OUTPUT: list of tuples, formatted as with ``start_gens``

        EXAMPLES::

            sage: R = ModularFormsRing(Gamma1(4))
            sage: R._find_generators(8, (), 2)
            [(2, 1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^9), 1 + 24*q^2 + 24*q^4 + O(q^6)),
             (2, q + 4*q^3 + 6*q^5 + 8*q^7 + O(q^9), q + 4*q^3 + 6*q^5 + O(q^6)),
             (3, 1 + 12*q^2 + 64*q^3 + 60*q^4 + 160*q^6 + 384*q^7 + 252*q^8 + O(q^9), 1 + 12*q^2 + 64*q^3 + 60*q^4 + O(q^6)),
             (3, q + 4*q^2 + 8*q^3 + 16*q^4 + 26*q^5 + 32*q^6 + 48*q^7 + 64*q^8 + O(q^9), q + 4*q^2 + 8*q^3 + 16*q^4 + 26*q^5 + O(q^6))]
        """
        default_params = (start_gens == () and start_weight == 2)

        if default_params and self.__cached_maxweight != -1:
            verbose("Already know generators up to weight %s -- using those" % self.__cached_maxweight)

            if self.__cached_maxweight >= maxweight:
                return [(k, f, F) for k, f, F in self.__cached_gens if k <= maxweight]

            start_gens = self.__cached_gens
            start_weight = self.__cached_maxweight + 1

        if self.group().is_even():
            increment = 2
        else:
            increment = 1

        working_prec = self.modular_forms_of_weight(maxweight).sturm_bound()

        # parse the list of start gens
        G = []
        for x in start_gens:
            k, f, F = x
            if F is None and f.prec() < working_prec:
                raise ValueError("Need start gens to precision at least %s" % working_prec)
            elif f is None or f.prec() < working_prec:
                f = F.qexp(working_prec).change_ring(self.base_ring())
            G.append((k, f, F))

        k = start_weight
        if increment == 2 and (k % 2) == 1:
            k += 1

        while k <= maxweight:

            if self.modular_forms_of_weight(k).dimension() == 0:
                k += increment
                continue

            verbose('Looking at k = %s' % k)
            M = self.modular_forms_of_weight(k)

            # 1. Multiply together all forms in G that give an element
            #    of M.
            if G:
                F = _span_of_forms_in_weight(G, k, M.sturm_bound(), None, False)
            else:
                F = (self.base_ring() ** M.sturm_bound()).zero_submodule()

            # 2. If the dimension of the span of the result is equal
            #    to the dimension of M, increment k.
            if F.rank() == M.dimension():
                if self.base_ring().is_field() or F.index_in_saturation() == 1:
                    # TODO: Do something clever if the submodule's of the right
                    # rank but not saturated -- avoid triggering needless
                    # modular symbol computations.
                    verbose('Nothing new in weight %s' % k)
                    k += increment
                    continue

            # 3. If the dimension is less, compute a basis for G, and
            #    try adding basis elements of M into G.

            verbose("Known generators span a subspace of dimension %s of space of dimension %s" % (F.dimension(), M.dimension()))
            if self.base_ring() == ZZ:
                verbose("saturation index is %s" % F.index_in_saturation())

            t = verbose("Computing more modular forms at weight %s" % k)
            kprec = M.sturm_bound()
            if self.base_ring() == QQ:
                B = M.q_echelon_basis(working_prec)
            else:
                B = M.q_integral_basis(working_prec)
            t = verbose("done computing forms", t)
            V = F.ambient_module().submodule_with_basis([f.padded_list(kprec) for f in B])
            Q = V / F
            for q in Q.gens():
                try:
                    qc = V.coordinates(Q.lift(q))
                except AttributeError:
                    # work around a silly free module bug
                    qc = V.coordinates(q.lift())
                qcZZ = [ZZ(_) for _ in qc] # lift to ZZ so we can define F
                f = sum([B[i] * qcZZ[i] for i in range(len(B))])
                F = M(f)
                G.append((k, f.change_ring(self.base_ring()), F))

            verbose('added %s new generators' % Q.ngens(), t)
            k += increment

        if default_params:
            self.__cached_maxweight = maxweight
            self.__cached_gens = G

        return G

    @cached_method
    def q_expansion_basis(self, weight, prec=None, use_random=True):
        r"""
        Return a basis of `q`-expansions for the space of modular forms
        of the given weight for this group, calculated using the ring
        generators given by ``find_generators``.

        INPUT:

        - ``weight`` -- the weight
        - ``prec`` -- integer (default: ``None``); power series
          precision. If ``None``, the precision defaults to the Sturm
          bound for the requested level and weight.
        - ``use_random`` -- boolean (default: ``True``); whether or not to
          use a randomized algorithm when building up the space of forms
          at the given weight from known generators of small weight.

        EXAMPLES::

            sage: m = ModularFormsRing(Gamma0(4))
            sage: m.q_expansion_basis(2,10)
            [1 + 24*q^2 + 24*q^4 + 96*q^6 + 24*q^8 + O(q^10),
             q + 4*q^3 + 6*q^5 + 8*q^7 + 13*q^9 + O(q^10)]
            sage: m.q_expansion_basis(3,10)
            []

            sage: X = ModularFormsRing(SL2Z)
            sage: X.q_expansion_basis(12, 10)
            [1 + 196560*q^2 + 16773120*q^3 + 398034000*q^4 + 4629381120*q^5 + 34417656000*q^6 + 187489935360*q^7 + 814879774800*q^8 + 2975551488000*q^9 + O(q^10),
             q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 + O(q^10)]

        We calculate a basis of a massive modular forms space, in two ways.
        Using this module is about twice as fast as Sage's generic code. ::

            sage: A = ModularFormsRing(11).q_expansion_basis(30, prec=40) # long time (5s)
            sage: B = ModularForms(Gamma0(11), 30).q_echelon_basis(prec=40) # long time (9s)
            sage: A == B # long time
            True

        Check that absurdly small values of ``prec`` don't mess things up::

            sage: ModularFormsRing(11).q_expansion_basis(10, prec=5)
            [1 + O(q^5), q + O(q^5), q^2 + O(q^5), q^3 + O(q^5),
             q^4 + O(q^5), O(q^5), O(q^5), O(q^5), O(q^5), O(q^5)]
        """
        d = self.modular_forms_of_weight(weight).dimension()
        if d == 0:
            return []

        if prec is None:
            prec = self.modular_forms_of_weight(weight).sturm_bound()

        working_prec = max(prec, self.modular_forms_of_weight(weight).sturm_bound())

        gen_weight = min(6, weight)

        while True:
            verbose("Trying to generate the %s-dimensional space at weight %s using generators of weight up to %s" % (d, weight, gen_weight))
            G = self.generators(maxweight=gen_weight, prec=working_prec)
            V = _span_of_forms_in_weight(G, weight, prec=working_prec, use_random=use_random, stop_dim=d)
            if V.rank() == d and (self.base_ring().is_field() or V.index_in_saturation() == 1):
                break
            else:
                gen_weight += 1
                verbose("Need more generators: trying again with generators of weight up to %s" % gen_weight)

        R = G[0][1].parent()
        return [R(list(x), prec=prec) for x in V.gens()]

    def cuspidal_ideal_generators(self, maxweight=8, prec=None):
        r"""
        Return a set of generators for the ideal of cuspidal forms in
        this ring, as a module over the whole ring.

        EXAMPLES::

            sage: ModularFormsRing(Gamma0(3)).cuspidal_ideal_generators(maxweight=12)
            [(6, q - 6*q^2 + 9*q^3 + 4*q^4 + O(q^5), q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 + O(q^6))]
            sage: [k for k,f,F in ModularFormsRing(13, base_ring=ZZ).cuspidal_ideal_generators(maxweight=14)]
            [4, 4, 4, 6, 6, 12]
        """
        working_prec = self.modular_forms_of_weight(maxweight).sturm_bound()

        if self.__cached_cusp_maxweight > -1:
            k = self.__cached_cusp_maxweight + 1
            verbose("Already calculated cusp gens up to weight %s -- using those" % (k-1))

            # we may need to increase the precision of the cached cusp
            # generators
            G = []
            for j,f,F in self.__cached_cusp_gens:
                if f.prec() >= working_prec:
                    f = F.qexp(working_prec).change_ring(self.base_ring())
                G.append((j, f, F))
        else:
            k = 2
            G = []

        while k <= maxweight:
            t = verbose(f"Looking for cusp generators in weight {k}")

            kprec = self.modular_forms_of_weight(k).sturm_bound()

            flist = []

            for (j, f, F) in G:
                for g in self.q_expansion_basis(k - j, prec=kprec):
                    flist.append(g*f)
            A = self.base_ring() ** kprec
            W = A.span([A(f.padded_list(kprec)) for f in flist])

            S = self.modular_forms_of_weight(k).cuspidal_submodule()
            if (W.rank() == S.dimension()
                and (self.base_ring().is_field() or W.index_in_saturation() == 1)):
                    verbose("Nothing new in weight %s" % k, t)
                    k += 1
                    continue

            t = verbose("Known cusp generators span a submodule of dimension %s of space of dimension %s" % (W.rank(), S.dimension()), t)

            B = S.q_integral_basis(prec=working_prec)
            V = A.span([A(f.change_ring(self.base_ring()).padded_list(kprec)) for f in B])
            Q = V/W

            for q in Q.gens():
                try:
                    qc = V.coordinates(Q.lift(q))
                except AttributeError:
                    # work around a silly free module bug
                    qc = V.coordinates(q.lift())
                qcZZ = [ZZ(_) for _ in qc] # lift to ZZ so we can define F
                f = sum([B[i] * qcZZ[i] for i in range(len(B))])
                F = S(f)
                G.append((k, f.change_ring(self.base_ring()), F))

            verbose('added %s new generators' % Q.ngens(), t)
            k += 1

        self.__cached_cusp_maxweight = maxweight
        self.__cached_cusp_gens = G

        if prec is None:
            return G
        elif prec <= working_prec:
            return [(k, f.truncate_powerseries(prec), F) for k,f,F in G]
        else:
            # user wants increased precision, so we may as well cache that
            Gnew = [(k, F.qexp(prec).change_ring(self.base_ring()), F) for k, f, F in G]
            self.__cached_cusp_gens = Gnew
            return Gnew

    def cuspidal_submodule_q_expansion_basis(self, weight, prec=None):
        r"""
        Return a basis of `q`-expansions for the space of cusp forms of
        weight ``weight`` for this group.

        INPUT:

        - ``weight`` -- the weight
        - ``prec`` -- integer (default: ``None``) precision of
          `q`-expansions to return

        ALGORITHM: Uses the method :meth:`cuspidal_ideal_generators` to
        calculate generators of the ideal of cusp forms inside this ring. Then
        multiply these up to weight ``weight`` using the generators of the
        whole modular form space returned by :meth:`q_expansion_basis`.

        EXAMPLES::

            sage: R = ModularFormsRing(Gamma0(3))
            sage: R.cuspidal_submodule_q_expansion_basis(20)
            [q - 8532*q^6 - 88442*q^7 + O(q^8), q^2 + 207*q^6 + 24516*q^7 + O(q^8),
             q^3 + 456*q^6 + O(q^8), q^4 - 135*q^6 - 926*q^7 + O(q^8), q^5 + 18*q^6 + 135*q^7 + O(q^8)]

        We compute a basis of a space of very large weight, quickly (using this
        module) and slowly (using modular symbols), and verify that the answers
        are the same. ::

            sage: A = R.cuspidal_submodule_q_expansion_basis(80, prec=30)  # long time (1s on sage.math, 2013)
            sage: B = R.modular_forms_of_weight(80).cuspidal_submodule().q_expansion_basis(prec=30)  # long time (19s on sage.math, 2013)
            sage: A == B # long time
            True
        """
        d = self.modular_forms_of_weight(weight).cuspidal_submodule().dimension()
        if d == 0:
            return []

        minprec = self.modular_forms_of_weight(weight).sturm_bound()
        if prec is None:
            prec = working_prec = minprec
        else:
            working_prec = max(prec, minprec)

        gen_weight = min(6, weight)

        while True:
            verbose("Trying to generate the %s-dimensional cuspidal submodule at weight %s using generators of weight up to %s" % (d, weight, gen_weight))
            G = self.cuspidal_ideal_generators(maxweight=gen_weight, prec=working_prec)

            flist = []
            for (j, f, F) in G:
                for g in self.q_expansion_basis(weight - j, prec=working_prec):
                    flist.append(g*f)

            A = self.base_ring() ** working_prec
            W = A.span([A(f.padded_list(working_prec)) for f in flist])
            if W.rank() == d and (self.base_ring().is_field() or W.index_in_saturation() == 1):
                break
            else:
                gen_weight += 1
                verbose("Need more generators: trying again with generators of weight up to %s" % gen_weight)

        R = G[0][1].parent()
        return [R(list(x), prec=prec) for x in W.gens()]


# Deprecated functions
find_generators = deprecated_function_alias(31559, ModularFormsRing.generators)
basis_for_modform_space = deprecated_function_alias(31559, ModularFormsRing.q_expansion_basis)
