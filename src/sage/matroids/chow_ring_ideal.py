r"""
Chow ring ideals of matroids

AUTHORS:

- Shriya M
"""
# ****************************************************************************
#       Copyright (C) 2024 Shriya M <25shriya at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.combinat.posets.posets import Poset
from itertools import product


class ChowRingIdeal(MPolynomialIdeal):
    def matroid(self):
        r"""
        Return the matroid of the given Chow ring ideal.

        EXAMPLES::

            sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False)
            sage: ch.defining_ideal().matroid()
            U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
            {3: {{0, 1, 2, 3, 4, 5}}}
        """
        return self._matroid

    def _lattice_flats(self):
        r"""
        Return the ranks and chains of lattice of flats of the matroid.

        EXAMPLES::

            sage: ch = matroids.catalog.NonFano().chow_ring(QQ, True, 'atom-free')
            sage: ch.defining_ideal()._lattice_flats()
            ({frozenset({'a'}): 1, frozenset({'b'}): 1, frozenset({'c'}): 1,
            frozenset({'d'}): 1, frozenset({'e'}): 1, frozenset({'f'}): 1,
            frozenset({'g'}): 1, frozenset({'d', 'e'}): 2,
            frozenset({'d', 'f'}): 2, frozenset({'e', 'f'}): 2,
            frozenset({'a', 'b', 'f'}): 2, frozenset({'a', 'c', 'e'}): 2,
            frozenset({'a', 'd', 'g'}): 2, frozenset({'b', 'c', 'd'}): 2,
            frozenset({'b', 'e', 'g'}): 2, frozenset({'c', 'f', 'g'}): 2,
            frozenset({'a', 'b', 'c', 'd', 'e', 'f', 'g'}): 3},
            Set of chains of Finite poset containing 17 elements)
        """
        F = self._matroid.lattice_of_flats()
        H = F.hasse_diagram()
        H.delete_vertex(self._matroid.flats(0)[0])  # remove the empty flat
        lattice_flats = Poset(H)
        flats = list(lattice_flats)
        flats.sort(key=lambda X: (len(X), sorted(X)))
        ranks = {F: self._matroid.rank(F) for F in flats}
        chains = lattice_flats.chains()  # Only chains
        return (ranks, chains)

    def flats_to_generator_dict(self):
        r"""
        Return the corresponding generators of flats/groundset elements of
        Chow ring ideal.

        EXAMPLES::

            sage: ch = matroids.Uniform(4, 6).chow_ring(QQ, True, 'atom-free')
            sage: ch.defining_ideal().flats_to_generator_dict()
            {frozenset({0}): A0, frozenset({1}): A1, frozenset({2}): A2,
             frozenset({3}): A3, frozenset({4}): A4, frozenset({5}): A5,
             frozenset({0, 1}): A01, frozenset({0, 2}): A02,
             frozenset({0, 3}): A03, frozenset({0, 4}): A04,
             frozenset({0, 5}): A05, frozenset({1, 2}): A12,
             frozenset({1, 3}): A13, frozenset({1, 4}): A14,
             frozenset({1, 5}): A15, frozenset({2, 3}): A23,
             frozenset({2, 4}): A24, frozenset({2, 5}): A25,
             frozenset({3, 4}): A34, frozenset({3, 5}): A35,
             frozenset({4, 5}): A45, frozenset({0, 1, 2}): A012,
             frozenset({0, 1, 3}): A013, frozenset({0, 1, 4}): A014,
             frozenset({0, 1, 5}): A015, frozenset({0, 2, 3}): A023,
             frozenset({0, 2, 4}): A024, frozenset({0, 2, 5}): A025,
             frozenset({0, 3, 4}): A034, frozenset({0, 3, 5}): A035,
             frozenset({0, 4, 5}): A045, frozenset({1, 2, 3}): A123,
             frozenset({1, 2, 4}): A124, frozenset({1, 2, 5}): A125,
             frozenset({1, 3, 4}): A134, frozenset({1, 3, 5}): A135,
             frozenset({1, 4, 5}): A145, frozenset({2, 3, 4}): A234,
             frozenset({2, 3, 5}): A235, frozenset({2, 4, 5}): A245,
             frozenset({3, 4, 5}): A345,
             frozenset({0, 1, 2, 3, 4, 5}): A012345}
        """
        return dict(self._flats_generator)


class ChowRingIdeal_nonaug_fy(ChowRingIdeal):
    r"""
    The Chow ring ideal of a matroid `M` in Feitchner-Yuzvinsky presentation.

    The *Chow ring ideal* for a matroid `M` in Feitchner-Yuzvinsky presentation
    is defined as the ideal `(I_M + J_M)` of the polynomial ring

    .. MATH::

        R[x_{F_1}, \ldots, x_{F_k}],

    where

    - `F_1, \ldots, F_k` are the non-empty flats of `M`,
    - `I_M` is the Stanley-Reisner ideal, i.e., it is generated
      by products `x_{F_1}, \ldots, x_{F_t}` for subsets `\{F_1, \ldots, F_t\}`
      of flats that are not chains, and
    - `J_M` is the ideal generated by all linear forms

      .. MATH::

          \sum_{a \in F} x_F

      for all atoms `a` in the lattice of flats.

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring

    REFERENCES:

    - [ANR2023]_

    EXAMPLES:

    Chow ring ideal of uniform matroid of rank 3 on 6 elements::

        sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False, 'fy')
        sage: ch.defining_ideal()
        Chow ring ideal of U(3, 6): Matroid of rank 3 on 6 elements with
         circuit-closures {3: {{0, 1, 2, 3, 4, 5}}} - non augmented
        sage: ch = matroids.catalog.Fano().chow_ring(QQ, False, 'fy')
        sage: ch.defining_ideal()
        Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements,
         type (3, 0) - non augmented in Feitchner-Yuzvinsky presentation
    """
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: I = matroids.catalog.Fano().chow_ring(QQ, False, 'fy').defining_ideal()
            sage: TestSuite(I).run(skip="_test_category")
        """
        self._matroid = M
        flats = [X for i in range(1, self._matroid.rank() + 1)
                 for X in self._matroid.flats(i)]
        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]
        try:
            poly_ring = PolynomialRing(R, names)  # self.ring
        except ValueError:  # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(flats))
        gens = poly_ring.gens()
        self._flats_generator = dict(zip(flats, gens))
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))

    def _gens_constructor(self, poly_ring):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: ch = matroids.catalog.NonFano().chow_ring(QQ, False)
            sage: sorted(ch.defining_ideal()._gens_constructor(ch.defining_ideal().ring()))
            [Ag + Aadg + Abeg + Acfg + Aabcdefg,
             Af + Aabf + Acfg + Adf + Aef + Aabcdefg,
             Ae + Aace + Abeg + Ade + Aef + Aabcdefg,
             Ad + Aadg + Abcd + Ade + Adf + Aabcdefg,
             Ac + Aace + Abcd + Acfg + Aabcdefg,
             Ab + Aabf + Abcd + Abeg + Aabcdefg,
             Aa + Aabf + Aace + Aadg + Aabcdefg,
             Adf*Aef, Ade*Aef, Acfg*Aef, Abeg*Aef, Abcd*Aef, Aadg*Aef,
             Aace*Aef, Aabf*Aef, Ag*Aef, Ad*Aef, Ac*Aef, Ab*Aef, Aa*Aef,
             Ade*Adf, Acfg*Adf, Abeg*Adf, Abcd*Adf, Aadg*Adf, Aace*Adf,
             Aabf*Adf, Ag*Adf, Ae*Adf, Ac*Adf, Ab*Adf, Aa*Adf, Acfg*Ade,
             Abeg*Ade, Abcd*Ade, Aadg*Ade, Aace*Ade, Aabf*Ade, Ag*Ade, Af*Ade,
             Ac*Ade, Ab*Ade, Aa*Ade, Abeg*Acfg, Abcd*Acfg, Aadg*Acfg,
             Aace*Acfg, Aabf*Acfg, Ae*Acfg, Ad*Acfg, Ab*Acfg, Aa*Acfg,
             Abcd*Abeg, Aadg*Abeg, Aace*Abeg, Aabf*Abeg, Af*Abeg, Ad*Abeg,
             Ac*Abeg, Aa*Abeg, Aadg*Abcd, Aace*Abcd, Aabf*Abcd, Ag*Abcd,
             Af*Abcd, Ae*Abcd, Aa*Abcd, Aace*Aadg, Aabf*Aadg, Af*Aadg, Ae*Aadg,
             Ac*Aadg, Ab*Aadg, Aabf*Aace, Ag*Aace, Af*Aace, Ad*Aace, Ab*Aace,
             Ag*Aabf, Ae*Aabf, Ad*Aabf, Ac*Aabf, Af*Ag, Ae*Ag, Ad*Ag, Ac*Ag,
             Ab*Ag, Aa*Ag, Ae*Af, Ad*Af, Ac*Af, Ab*Af, Aa*Af, Ad*Ae, Ac*Ae,
             Ab*Ae, Aa*Ae, Ac*Ad, Ab*Ad, Aa*Ad, Ab*Ac, Aa*Ac, Aa*Ab]
        """
        flats = list(self._flats_generator)
        lattice_flats = Poset((flats, lambda x, y: x <= y))
        I = []
        subsets = lattice_flats.antichains().elements_of_depth_iterator(2)
        for subset in subsets:
            term = poly_ring.one()
            for el in subset:
                term *= self._flats_generator[el]
            I.append(term)  # Stanley-Reisner Ideal
        atoms = self._matroid.lattice_of_flats().atoms()
        atoms_gen = {a: poly_ring.zero() for a in atoms}
        for F in flats:
            for a in atoms:
                if a.issubset(F):
                    atoms_gen[a] += self._flats_generator[F]
        J = list(atoms_gen.values())  # Linear generators
        return I + J

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ch = matroids.catalog.Fano().chow_ring(QQ, False)
            sage: ch.defining_ideal()
            Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements,
            type (3, 0) - non augmented in Feitchner-Yuzvinsky presentation
        """
        return "Chow ring ideal of {} - non augmented in Feitchner-Yuzvinksy presentation".format(self._matroid)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: M1 = Matroid(groundset='abcd', bases=['ab','ad', 'bc'])
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch.defining_ideal()._latex_()
            '(I_{\\text{\\texttt{Matroid{ }of{ }rank{ }2{ }on{ }4{ }elements{ }with{ }3{ }bases}}} + J_{\\text{\\texttt{Matroid{ }of{ }rank{ }2{ }on{ }4{ }elements{ }with{ }3{ }bases}}}'
        """
        from sage.misc.latex import latex
        return '(I_{{{M}}} + J_{{{M}}}'.format(M=latex(self._matroid))

    def groebner_basis(self, algorithm='', *args, **kwargs):
        r"""
        Return a Groebner basis of ``self``.

        EXAMPLES::

            sage: ch = Matroid(groundset='abc', bases=['ab', 'ac']).chow_ring(QQ, False)
            sage: ch.defining_ideal().groebner_basis()
            [Aa*Abc, Aa + Aabc, Abc + Aabc, Aa*Aabc, Abc*Aabc, Aabc^2]
            sage: ch.defining_ideal().groebner_basis().is_groebner()
            True
            sage: ch.defining_ideal().hilbert_series() == ch.defining_ideal().gens().ideal().hilbert_series()
            True

        Another example would be the Groebner basis of the Chow ring ideal of
        the matroid of the length 3 cycle graph::

            sage: ch = Matroid(graphs.CycleGraph(3)).chow_ring(QQ, False)
            sage: ch.defining_ideal().groebner_basis()
            [A0*A1, A0*A2, A1*A2, A0 + A3, A1 + A3, A2 + A3, A0*A3, A1*A3, A2*A3, A3^2]
            sage: ch.defining_ideal().groebner_basis().is_groebner()
            True
            sage: ch.defining_ideal().hilbert_series() == ch.defining_ideal().gens().ideal().hilbert_series()
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().groebner_basis(algorithm=algorithm, *args, **kwargs)
        flats = sorted(list(self._flats_generator), key=len)
        ranks = {F: self._matroid.rank(F) for F in flats}
        gb = []
        R = self.ring()
        flats_gen = self._flats_generator
        lattice_flats = Poset((flats, lambda x, y: x <= y))
        antichains = lattice_flats.antichains().elements_of_depth_iterator(2)
        for subset in antichains:  # Taking antichains of size 2
            term = R.one()
            for x in subset:
                term *= flats_gen[x]
            gb.append(term)
        for F in flats:  # Reduced groebner basis by computing the sum first and then the product
            term = R.zero()
            for G in lattice_flats.order_filter([F]):
                term += flats_gen[G]
            for G in lattice_flats.order_ideal([F]):
                if G != F:
                    gb.append(flats_gen[G]*(term) ** (ranks[F] - ranks[G]))

            gb.append(term ** ranks[F])

        return PolynomialSequence(R, [gb])

    def normal_basis(self, algorithm='', *args, **kwargs):
        r"""
        Return the monomial basis of the quotient ring of this ideal.

        EXAMPLES::

            sage: ch = matroids.Z(3).chow_ring(QQ, False)
            sage: I = ch.defining_ideal()
            sage: I.normal_basis()
            [1, Ax2x3y1, Ax1x3y2, Ax1x2y3, Ay1y2y3, Atx1y1, Atx2y2, Atx3y3, Atx1x2x3y1y2y3, Atx1x2x3y1y2y3^2]
            sage: set(I.gens().ideal().normal_basis()) == set(I.normal_basis())
            True
            sage: ch = matroids.AG(2,3).chow_ring(QQ, False)
            sage: I = ch.defining_ideal()
            sage: I.normal_basis()
            [1, A012, A236, A046, A156, A345, A247, A057, A137, A258, A678,
            A038, A148, A012345678, A012345678^2]
            sage: set(I.gens().ideal().normal_basis()) == set(I.normal_basis())
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().normal_basis(algorithm=algorithm, *args, **kwargs)
        R = self.ring()
        flats_gen = self._flats_generator
        monomial_basis = []
        ranks, chains = self._lattice_flats()
        for subset in chains:
            max_powers = []
            k = len(subset)
            if (k == 0):
                monomial_basis.append(R.one())
            elif not ((k == 1) & (ranks[subset[0]] == 1)):
                for i in range(k):
                    if i == 0:
                        max_powers.append(ranks[subset[i]])
                    else:
                        max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                for combination in product(*(range(1, p) for p in max_powers)):
                    expression = R.one()
                    for val, c in zip(subset, combination):
                        expression *= flats_gen[val] ** c
                    monomial_basis.append(expression)
        return PolynomialSequence(R, [monomial_basis])

class ChowRingIdeal_nonaug_af(ChowRingIdeal):
    r"""
    The Chow ring ideal of a matroid `M` in atom-free presentation.

    The *Chow ring ideal* for a matroid `M` in atom-free presentation
    is defined as the ideal `(I_M + J_M + K_M)` of the polynomial ring

    .. MATH::

        R[x_{F_1}, \ldots, x_{F_k}],

    where

    - `F_1, \ldots, F_k` are flats of `M` of at least rank 2,
    - `I_M` is the ideal generated by products `x_{F_i} x_{F_j}` for
      incomparable flats `F_i` and `F_j`,
    - `J_M` is the ideal generated by all linear forms

      .. MATH::

          x_F \sum_{F' \superseteq F \vee i} x_{F'}

      for all flats `F` and `i` in `E - F`, where `E` is the groundset of
      the matroid, and
    - `K_M` is the ideal generated by all quadratic forms

      .. MATH::

          \sum_{F \superseteq i \vee j} x_F^2 + \sum_{F' \superset F \superseteq i \vee j} 2 x_F x_{F'}

      for all elements `i, j` in `E`, such that `i \neq j`.

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring

    REFERENCES:

    - [MM2022]_

    EXAMPLES:

    Chow ring ideal of uniform matroid of rank 3 on 6 elements::

        sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False, 'atom-free')
        sage: ch.defining_ideal()
        Chow ring ideal of U(3, 6): Matroid of rank 3 on 6 elements with
        circuit-closures {3: {{0, 1, 2, 3, 4, 5}}} - non augmented in the
        atom-free presentation
        sage: ch = matroids.catalog.Fano().chow_ring(QQ, False)
        sage: ch.defining_ideal()
        Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements,
         type (3, 0) - non augmented in the atom-free presentation
    """
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: I = matroids.catalog.Fano().chow_ring(QQ, False, 'atom-free').defining_ideal()
            sage: TestSuite(I).run(skip="_test_category")
        """
        self._matroid = M
        flats = [X for i in range(2, self._matroid.rank() + 1)
                 for X in self._matroid.flats(i)]
        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]
        try:
            poly_ring = PolynomialRing(R, names)  # self.ring
        except ValueError:  # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(flats))
        gens = poly_ring.gens()
        self._flats_generator = dict(zip(flats, gens))
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))

    def _gens_constructor(self, poly_ring):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: ch = matroids.catalog.NonFano().chow_ring(QQ, False)
            sage: sorted(ch.defining_ideal()._gens_constructor(ch.defining_ideal().ring()))
            [0, Adef*Aabcdefg, Adef*Aabcdefg, Adef*Aabcdefg, 4*Adef*Aabcdefg,
             Acfg*Aabcdefg, Acfg*Aabcdefg, Acfg*Aabcdefg, 4*Acfg*Aabcdefg,
             Abeg*Aabcdefg, Abeg*Aabcdefg, Abeg*Aabcdefg, 4*Abeg*Aabcdefg,
             Abcd*Aabcdefg, Abcd*Aabcdefg, Abcd*Aabcdefg, 4*Abcd*Aabcdefg,
             Aadg*Aabcdefg, Aadg*Aabcdefg, Aadg*Aabcdefg,
             Aadg*Aabcdefg + Abeg*Aabcdefg + Acfg*Aabcdefg,
             Aadg*Aabcdefg + Abcd*Aabcdefg + Adef*Aabcdefg,
             4*Aadg*Aabcdefg, Aace*Aabcdefg, Aace*Aabcdefg, Aace*Aabcdefg,
             Aace*Aabcdefg + Abeg*Aabcdefg + Adef*Aabcdefg,
             Aace*Aabcdefg + Abcd*Aabcdefg + Acfg*Aabcdefg,
             4*Aace*Aabcdefg, Aabf*Aabcdefg, Aabf*Aabcdefg, Aabf*Aabcdefg,
             Aabf*Aabcdefg + Acfg*Aabcdefg + Adef*Aabcdefg,
             Aabf*Aabcdefg + Abcd*Aabcdefg + Abeg*Aabcdefg,
             Aabf*Aabcdefg + Aace*Aabcdefg + Aadg*Aabcdefg, 4*Aabf*Aabcdefg,
             Adef^2 + Aabcdefg^2, Adef^2 + Aabcdefg^2, Adef^2 + Aabcdefg^2,
             Acfg*Adef, Abeg*Adef, Abcd*Adef, Aadg*Adef, Aace*Adef, Aabf*Adef,
             Acfg^2 + Aabcdefg^2, Acfg^2 + Aabcdefg^2, Acfg^2 + Aabcdefg^2,
             Abeg*Acfg, Abcd*Acfg, Aadg*Acfg, Aace*Acfg, Aabf*Acfg,
             Abeg^2 + Aabcdefg^2, Abeg^2 + Aabcdefg^2, Abeg^2 + Aabcdefg^2,
             Abcd*Abeg, Aadg*Abeg, Aace*Abeg, Aabf*Abeg, Abcd^2 + Aabcdefg^2,
             Abcd^2 + Aabcdefg^2, Abcd^2 + Aabcdefg^2, Aadg*Abcd, Aace*Abcd,
             Aabf*Abcd, Aadg^2 + Aabcdefg^2, Aadg^2 + Aabcdefg^2,
             Aadg^2 + Aabcdefg^2, Aadg^2 + Abeg^2 + Acfg^2 + Aabcdefg^2, 
             Aadg^2 + Abcd^2 + Adef^2 + Aabcdefg^2, Aace*Aadg, Aabf*Aadg,
             Aace^2 + Aabcdefg^2, Aace^2 + Aabcdefg^2, Aace^2 + Aabcdefg^2,
             Aace^2 + Abeg^2 + Adef^2 + Aabcdefg^2,
             Aace^2 + Abcd^2 + Acfg^2 + Aabcdefg^2, Aabf*Aace,
             Aabf^2 + Aabcdefg^2, Aabf^2 + Aabcdefg^2,
             Aabf^2 + Aabcdefg^2, Aabf^2 + Acfg^2 + Adef^2 + Aabcdefg^2,
             Aabf^2 + Abcd^2 + Abeg^2 + Aabcdefg^2,
             Aabf^2 + Aace^2 + Aadg^2 + Aabcdefg^2]
        """
        E = list(self._matroid.groundset())
        flats = list(self._flats_generator)
        lattice_flats = Poset((flats, lambda x, y: x <= y))
        I = []
        flats_gen = self._flats_generator
        subsets = lattice_flats.antichains().elements_of_depth_iterator(2)
        for subset in subsets:
            I.append(flats_gen[subset[0]] * flats_gen[subset[1]])  # Stanley-Reisner Ideal
        J = []
        for F in flats:
            term = poly_ring.zero()
            for i in [el for el in E if el not in F]:
                for G in lattice_flats.order_filter([F]):
                    if G >= F.union(frozenset({i})):
                        term += flats_gen[G]
            J.append(flats_gen[F] * term)
        K = []
        for x, i in enumerate(E):
            for j in E[x:]:
                term1 = poly_ring.zero()
                term2 = poly_ring.zero()
                for F in flats:
                    if F >= frozenset({i}).union(frozenset({j})):
                        term1 += flats_gen[F] ** 2
                        for G in lattice_flats.order_filter([F]):
                            if G != F:
                                term2 += flats_gen[F]*flats_gen[G]
                K.extend([term1, term2])
        return I + J + K

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ch = matroids.catalog.Fano().chow_ring(QQ, False, 'atom-free')
            sage: ch.defining_ideal()
            Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements,
            type (3, 0) - non augmented in the atom-free presentation
        """
        return "Chow ring ideal of {} - non augmented in the atom-free presentation".format(self._matroid)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: M1 = Matroid(groundset='abcd', bases=['ab','ad', 'bc'])
            sage: ch = M1.chow_ring(QQ, False, 'atom-free')
            sage: ch.defining_ideal()._latex_()
            '(I_{\\text{\\texttt{Matroid{ }of{ }rank{ }2{ }on{ }4{ }elements{ }with{ }3{ }bases}}} + J_{\\text{\\texttt{Matroid{ }of{ }rank{ }2{ }on{ }4{ }elements{ }with{ }3{ }bases}}} + K_{\\text{\\texttt{Matroid{ }of{ }rank{ }2{ }on{ }4{ }elements{ }with{ }3{ }bases}}}'
        """
        from sage.misc.latex import latex
        return '(I_{{{M}}} + J_{{{M}}} + K_{{{M}}}'.format(M=latex(self._matroid))

    def groebner_basis(self, algorithm = '', *args, **kwargs):
        r"""
        Return a Groebner basis of ``self``.

        EXAMPLES::

            sage: ch = matroids.catalog.Fano().chow_ring(QQ, False, 'atom-free')
            sage: ch.defining_ideal().groebner_basis()
            Polynomial Sequence with 36 Polynomials in 8 Variables
            sage: ch.defining_ideal().groebner_basis().is_groebner()
            True
            sage: ch.defining_ideal().hilbert_series() == ch.defining_ideal().gens().ideal().hilbert_series()
            True

        Another example would be the Groebner basis of the Chow ring ideal of
        the matroid of the length 3 cycle graph::

            sage: ch = Matroid(graphs.CycleGraph(3)).chow_ring(QQ, False, 'atom-free')
            sage: ch.defining_ideal().groebner_basis()
            [A^2]
            sage: ch.defining_ideal().groebner_basis().is_groebner()
            True
            sage: ch.defining_ideal().hilbert_series() == ch.defining_ideal().gens().ideal().hilbert_series()
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().groebner_basis(algorithm=algorithm, *args, **kwargs)
        flats = list(self._flats_generator)
        lattice_flats = Poset((flats, lambda x, y: x <= y))
        gb = []
        poly_ring = self.ring()
        flats_gen = self._flats_generator
        subsets = lattice_flats.antichains().elements_of_depth_iterator(2)
        for subset in subsets:
            gb.append(flats_gen[subset[0]] * flats_gen[subset[1]])
        subsets = lattice_flats.chains().elements_of_depth_iterator(2)
        ranks, chains = self._lattice_flats()
        for subset in subsets:
            term = poly_ring.zero()
            for F in lattice_flats.order_filter([subset[1]]):
                term += flats_gen[F]
            gb.append(flats_gen[subset[0]] * (term ** (ranks[subset[1]] - ranks[subset[0]])))
        for F in flats:
            term = poly_ring.zero()
            for G in lattice_flats.order_filter([F]):
                term += flats_gen[G]
            gb.append(term ** ranks[F])
        return PolynomialSequence(poly_ring, [gb])

    def normal_basis(self, algorithm = '', *args, **kwargs):
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().normal_basis(algorithm=algorithm, *args, **kwargs)
        R = self.ring()
        flats_gen = self._flats_generator
        monomial_basis = []
        ranks, chains = self._lattice_flats()
        for subset in chains:
            max_powers = []
            k = len(subset)
            if not subset:
                monomial_basis.append(R.one())
            else:
                for i in range(k):
                    if i == 0:
                        max_powers.append(ranks[subset[i]])
                    else:
                        max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                ranges = [range(1, p) for p in max_powers]
                first_rank = ranks[subset[k-1]]
                for combination in product(*(r for r in ranges)):
                    # Generating combinations for all powers from 1 to max_powers
                    if sum(combination) <= first_rank:
                        expression = R.one()
                        for val, c in zip(subset, combination):
                            expression *= flats_gen[val] ** c
                        monomial_basis.append(expression)
        return PolynomialSequence(R, [monomial_basis])

class ChowRingIdeal_nonaug_sp(ChowRingIdeal):
    def __init__(self, M, R):
        self._matroid = M
        flats = [X for i in range(1, self._matroid.rank() + 1)
                 for X in self._matroid.flats(i)]
        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]
        try:
            poly_ring = PolynomialRing(R, names)  # self.ring
        except ValueError:  # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(flats))
        gens = poly_ring.gens()
        self._flats_generator = dict(zip(flats, gens))
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))

    def _gens_constructor(self, poly_ring):
        flats = list(self._flats_generator)
        lattice_flats = Poset((flats, lambda x, y: x <= y))
        flats_gen = self._flats_generator
        ranks, chains = self._lattice_flats()
        atoms = [F for F in flats if ranks[F] == 1]
        I = [flats_gen[a] for a in atoms]
        J = []
        subsets = lattice_flats.antichains().elements_of_depth_iterator(2)
        for subset in subsets:
            F = subset[0]
            G = subset[1]
            term1 = poly_ring.zero()
            term2 = poly_ring.zero()
            for H in lattice_flats.order_filter([F]):
                term1 += lattice_flats.moebius_function(F, H) * flats_gen[H]
            for H in lattice_flats.order_filter([G]):
                term2 += lattice_flats.moebius_function(G, H) * flats_gen[H]
            J.append(term1 * term2)
        return I + J

    def _repr_(self):
        return "Chow ring ideal of {} - non augmented in simplicial presentation".format(self._matroid)

    def _latex_(self):
        from sage.misc.latex import latex
        return '(I_{{{M}}} + J_{{{M}}}'.format(M=latex(self._matroid))

    def groebner_basis(self, algorithm = '', *args, **kwargs):
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().groebner_basis(algorithm=algorithm, *args, **kwargs)
        flats = list(self._flats_generator)
        lattice_flats = Poset((flats, lambda x, y: x <= y))
        gb = []
        poly_ring = self.ring()
        flats_gen = self._flats_generator
        ranks, chains = self._lattice_flats()
        subsets = lattice_flats.antichains().elements_of_depth_iterator(2)
        for subset in subsets:
            F = subset[0]
            G = subset[1]
            term1 = poly_ring.zero()
            term2 = poly_ring.zero()
            for H in lattice_flats.order_filter([F]):
                term1 += lattice_flats.moebius_function(F, H) * flats_gen[H]
            for H in lattice_flats.order_filter([G]):
                term2 += lattice_flats.moebius_function(G, H) * flats_gen[H]
            gb.append(term1 * term2)
        subsets = lattice_flats.chains().elements_of_depth_iterator(2)
        for subset in subsets:
            F = subset[0]
            G = subset[1]
            term = poly_ring.zero()
            for H in lattice_flats.order_filter([F]):
                term += lattice_flats.moebius_function(F, H) * flats_gen[H]
            gb.append(term * (flats_gen[G] ** (ranks[G] - ranks[F])))

        for F in flats:
            gb.append(flats_gen[F] ** ranks[F])
        return PolynomialSequence(poly_ring, [gb])
            
    def normal_basis(self, algorithm = '', *args, **kwargs):
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().normal_basis(algorithm=algorithm, *args, **kwargs)
        r = self._matroid.rank() - 1
        R = self.ring()
        flats_gen = self._flats_generator
        monomial_basis = []
        ranks, chains = self._lattice_flats()
        for subset in chains:
            if not subset:
                monomial_basis.append(R.one())
            else:
                k = len(subset)
                max_powers = []
                max_powers.append(ranks[subset[0]])
                for i in range(1, k):
                    max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                ranges = [range(1, p) for p in max_powers]
                ranges[0] = range(1, max_powers[0] + 1)
                for combination in product(*(r for r in ranges)):
                    # generating combinations for all powers up to max_powers
                    expression = R.one()
                    for val, c in zip(subset, combination):
                        expression *= flats_gen[val] ** c
                    if sum(combination) <= r:
                        monomial_basis.append(expression)
        return PolynomialSequence(R, [monomial_basis])


class AugmentedChowRingIdeal_fy(ChowRingIdeal):
    r"""
    The augmented Chow ring ideal of matroid `M` over ring `R` in
    the Feitchner-Yuzvinsky presentation.

    The augmented Chow ring ideal for a matroid `M` is defined as the ideal
    `(I_M + J_M)` of the following polynomial ring

    .. MATH::

        R[y_{e_1}, \ldots, y_{e_n}, x_{F_1}, \ldots, x_{F_k}],

    where

    - `F_1, \ldots, F_k` are the proper flats of `M`,
    - `e_1, \ldots, e_n` are `n` elements of groundset of `M`,
    - `J_M` is the ideal generated by all quadratic monomials `x_{F} x_{F'}`,
      where `F` and `F'` are incomparable elements in the lattice of
      flats and `y_{i} x_F` for all flats `F` and `i \in E \setminus F` and
    - `I_M` is the ideal generated by all linear forms

      .. MATH::

          y_i - \sum_{i \notin F} x_F

      for all `i \in E`.

    The augmented Chow ring ideal in the Feitchner-Yuzvinsky presentation
    for a simple matroid `M` is defined as the ideal `I_{FY}(M)` of the
    following polynomial ring

    .. MATH::

        R[y_{e_1}, \ldots, y_{e_n}, y_{F_1 \cup e}, \ldots, y_{F_k \cup e}],

    where `F_1, \ldots, F_k` are the flats of `M`, `e_1, \ldots, e_n` are
    `n` elements of groundset of `M`, and `I_{FY}(M)` is the ideal generated by

    - all quadratic monomials `y_{F \cup e} y_{F' \cup e}`, for incomparable
      elements `F` and `F'` in the lattice of flats,

    - `y_{i} y_{F \cup e}` for all flats `F` and all `i \in E \setminus F`

    - for all `i \in E`

      .. MATH::

          y_i + \sum_{i \in F} y_{F \cup e}

    - and

      .. MATH::

          \sum_{F} y_{F \cup e}.

    Setting `x_F = y_{F \cup e}` and using the last linear
    form to eliminate `x_E` recovers the usual presentation of
    augmented Chow ring of `M`.

    REFERENCES:

    - [MM2022]_

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring

    EXAMPLES:

    Augmented Chow ring ideal of Wheel matroid of rank 3::

        sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'fy')
        sage: ch.defining_ideal()
        Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 6
        elements with 16 bases of Feitchner-Yuzvinsky presentation
    """
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: I = matroids.Wheel(3).chow_ring(QQ, True, 'fy').defining_ideal()
            sage: TestSuite(I).run(skip="_test_category")
        """
        self._matroid = M
        self._flats = [X for i in range(self._matroid.rank() + 1)
                       for X in self._matroid.flats(i)]
        E = list(self._matroid.groundset())
        self._flats_generator = dict()
        try:
            names_groundset = ['A{}'.format(''.join(str(x))) for x in E]
            names_flats = ['B{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self._flats]
            poly_ring = PolynomialRing(R, names_groundset + names_flats)  # self.ring()
        except ValueError:  # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(E) + len(self._flats))
        for i, x in enumerate(E):
            self._flats_generator[x] = poly_ring.gens()[i]
        for i, F in enumerate(self._flats):
            self._flats_generator[F] = poly_ring.gens()[len(E) + i]
        self._flats_containing = {x: [] for x in E}
        for F in self._flats:
            for x in F:
                self._flats_containing[x].append(F)
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))

    def _gens_constructor(self, poly_ring):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'fy')
            sage: sorted(ch.defining_ideal()._gens_constructor(ch.defining_ideal().ring()))
            [B + B0 + B1 + B2 + B3 + B4 + B5 + B013 + B025 + B04 + B124 + B15 + B23 + B345 + B012345,
             B + B0 + B1 + B2 + B3 + B4 + B5 + B013 + B025 + B04 + B124 + B15 + B23 + B345 + B012345,
             B + B0 + B1 + B2 + B3 + B4 + B5 + B013 + B025 + B04 + B124 + B15 + B23 + B345 + B012345,
             B + B0 + B1 + B2 + B3 + B4 + B5 + B013 + B025 + B04 + B124 + B15 + B23 + B345 + B012345,
             B + B0 + B1 + B2 + B3 + B4 + B5 + B013 + B025 + B04 + B124 + B15 + B23 + B345 + B012345,
             B + B0 + B1 + B2 + B3 + B4 + B5 + B013 + B025 + B04 + B124 + B15 + B23 + B345 + B012345,
             A5 + B5 + B025 + B15 + B345 + B012345,
             A4 + B4 + B04 + B124 + B345 + B012345,
             A3 + B3 + B013 + B23 + B345 + B012345,
             A2 + B2 + B025 + B124 + B23 + B012345,
             A1 + B1 + B013 + B124 + B15 + B012345,
             A0 + B0 + B013 + B025 + B04 + B012345,
             B23*B345, B15*B345, B124*B345, B04*B345, B025*B345, B013*B345,
             B2*B345, B1*B345, B0*B345, A2*B345, A1*B345, A0*B345, B15*B23,
             B124*B23, B04*B23, B025*B23, B013*B23, B5*B23, B4*B23, B1*B23,
             B0*B23, A5*B23, A4*B23, A1*B23, A0*B23, B124*B15, B04*B15,
             B025*B15, B013*B15, B4*B15, B3*B15, B2*B15, B0*B15, A4*B15,
             A3*B15, A2*B15, A0*B15, B04*B124, B025*B124, B013*B124, B5*B124,
             B3*B124, B0*B124, A5*B124, A3*B124, A0*B124, B025*B04, B013*B04,
             B5*B04, B3*B04, B2*B04, B1*B04, A5*B04, A3*B04, A2*B04, A1*B04,
             B013*B025, B4*B025, B3*B025, B1*B025, A4*B025, A3*B025, A1*B025,
             B5*B013, B4*B013, B2*B013, A5*B013, A4*B013, A2*B013, B4*B5,
             B3*B5, B2*B5, B1*B5, B0*B5, A4*B5, A3*B5, A2*B5, A1*B5, A0*B5,
             B3*B4, B2*B4, B1*B4, B0*B4, A5*B4, A3*B4, A2*B4, A1*B4, A0*B4,
             B2*B3, B1*B3, B0*B3, A5*B3, A4*B3, A2*B3, A1*B3, A0*B3, B1*B2,
             B0*B2, A5*B2, A4*B2, A3*B2, A1*B2, A0*B2, B0*B1, A5*B1, A4*B1,
             A3*B1, A2*B1, A0*B1, A5*B0, A4*B0, A3*B0, A2*B0, A1*B0, A5*B,
             A4*B, A3*B, A2*B, A1*B, A0*B]
        """
        E = list(self._matroid.groundset())
        Q = []
        L = []
        lattice_flats = Poset((self._flats, lambda x, y: x <= y))
        antichains = lattice_flats.antichains().elements_of_depth_iterator(2)
        for F, G in antichains:
            Q.append(self._flats_generator[F] * self._flats_generator[G])  # Quadratic generators

        for x in E:
            term = poly_ring.zero()
            term1 = poly_ring.zero()
            for F in self._flats:
                term1 += self._flats_generator[F]
                if F not in self._flats_containing[x]:
                    Q.append(self._flats_generator[x] * self._flats_generator[F])
                else:
                    term += self._flats_generator[F]
            L.append(self._flats_generator[x] + term)  # Linear generators
            L.append(term1)
        return Q + L

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'fy')
            sage: ch.defining_ideal()
            Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on
            6 elements with 16 bases of Feitchner-Yuzvinsky presentation
        """
        return "Augmented Chow ring ideal of {} of Feitchner-Yuzvinsky presentation".format(self._matroid)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: M1 = Matroid(graphs.CycleGraph(3))
            sage: ch = M1.chow_ring(QQ, True, 'fy')
            sage: ch.defining_ideal()._latex_()
            'I_{FY}(\\text{\\texttt{Graphic{ }matroid{ }of{ }rank{ }2{ }on{ }3{ }elements}})'
        """
        from sage.misc.latex import latex
        return 'I_{{FY}}({})'.format((latex(self._matroid)))

    def groebner_basis(self, algorithm='', *args, **kwargs):
        r"""
        Return the Groebner basis of ``self``.

        EXAMPLES::

            sage: ch = matroids.catalog.NonFano().chow_ring(QQ, True, 'fy')
            sage: ch.defining_ideal().groebner_basis(algorithm='')
            Polynomial Sequence with 178 Polynomials in 25 Variables
            sage: ch.defining_ideal().groebner_basis(algorithm='').is_groebner()
            True
            sage: ch.defining_ideal().hilbert_series() == ch.defining_ideal().gens().ideal().hilbert_series()
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().groebner_basis(algorithm=algorithm, *args, **kwargs)
        gb = []  # reduced groebner basis with two eliminated cases
        E = list(self._matroid.groundset())
        poly_ring = self.ring()
        reln = lambda x,y: x <= y
        lattice_flats = Poset((self._flats, reln))
        antichains = lattice_flats.antichains().elements_of_depth_iterator(2)
        for F, G in antichains:
            gb.append(self._flats_generator[F] * self._flats_generator[G]) # non-nested flats
        for i in E:
            term = poly_ring.zero()
            for H in self._flats_containing[i]:
                term += self._flats_generator[H]
            if term != poly_ring.zero():
                gb.append(self._flats_generator[i] + term)  # 5.7 (MM2022)

        for F in self._flats:
            term1 = poly_ring.zero()
            for H in lattice_flats.order_filter([F]):
                term1 += self._flats_generator[H]
            if term1 != poly_ring.zero():
                gb.append(term1**(self._matroid.rank(F) + 1))  # 5.6 (MM2022)
                order_ideal_modified = lattice_flats.order_ideal([F])
                order_ideal_modified.remove(F)
                for G in order_ideal_modified:  # nested flats
                    gb.append(self._flats_generator[G]*term1**(self._matroid.rank(F) - self._matroid.rank(G)))

        return PolynomialSequence(poly_ring, [gb])

    def normal_basis(self, algorithm='', *args, **kwargs):
        r"""
        Return the monomial basis of the quotient ring of this ideal.

        EXAMPLES::

            sage: ch = matroids.Uniform(2, 5).chow_ring(QQ, True, 'fy')
            sage: I = ch.defining_ideal()
            sage: I.normal_basis()
            [1, B0, B1, B2, B3, B4, B01234, B01234^2]
            sage: set(I.gens().ideal().normal_basis()) == set(I.normal_basis())
            True
            sage: ch = matroids.catalog.Fano().chow_ring(QQ, True, 'fy')
            sage: I = ch.defining_ideal()
            sage: I.normal_basis()
            Polynomial Sequence with 32 Polynomials in 15 Variables
            sage: set(I.gens().ideal().normal_basis()) == set(I.normal_basis())
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().normal_basis(algorithm=algorithm, *args, **kwargs)
        R = self.ring()
        flats_gen = self._flats_generator
        monomial_basis = []
        ranks, chains = self._lattice_flats()
        for subset in chains:
            if not subset:
                monomial_basis.append(R.one())
            else:
                k = len(subset)
                max_powers = []
                max_powers.append(ranks[subset[0]])
                for i in range(1, k):
                    max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                ranges = [range(1, p) for p in max_powers]
                ranges[0] = range(1, max_powers[0] + 1)
                for combination in product(*(r for r in ranges)):
                    # generating combinations for all powers up to max_powers
                    expression = R.one()
                    for val, c in zip(subset, combination):
                        expression *= flats_gen[val] ** c
                    monomial_basis.append(expression)
        return PolynomialSequence(R, [monomial_basis])


class AugmentedChowRingIdeal_atom_free(ChowRingIdeal):
    r"""
    The augmented Chow ring ideal for a matroid `M` over ring `R` in the
    atom-free presentation.

    The augmented Chow ring ideal in the atom-free presentation for a matroid
    `M` is defined as the ideal `I_{af}(M)` of the polynomial ring:

    .. MATH::

        R[x_{F_1}, \ldots, x_{F_k}],

    where `F_1, \ldots, F_k` are the non-empty flats of `M` and `I_{af}(M)` is
    the ideal generated by

    - all quadratic monomials `x_{F} x_{F'}` for all incomparable elements
      `F` and `F'` in the lattice of flats,

    - for all flats `F` and `i \in E \setminus F`

      .. MATH::

        x_F \sum_{i \in F'} x_{F'}

    - and for all `i \in E`

      .. MATH::

        \sum_{i \in F'} (x_{F'})^2.

    REFERENCES:

    - [MM2022]_

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring

    EXAMPLES:

    Augmented Chow ring ideal of Wheel matroid of rank 3::

        sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
        sage: ch.defining_ideal()
        Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 6
        elements with 16 bases in the atom-free presentation
    """
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: I = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free').defining_ideal()
            sage: TestSuite(I).run(skip="_test_category")
        """
        self._matroid = M
        self._flats = [X for i in range(1, self._matroid.rank() + 1)
                       for X in self._matroid.flats(i)]
        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self._flats]
        try:
            poly_ring = PolynomialRing(R, names)  # self.ring
        except ValueError:  # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(self._flats))
        gens = poly_ring.gens()
        self._flats_generator = dict(zip(self._flats, gens))
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))

    def _gens_constructor(self, poly_ring):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: M1 = Matroid(graphs.CycleGraph(3))
            sage: ch = M1.chow_ring(QQ, True, 'atom-free')
            sage: sorted(ch.defining_ideal()._gens_constructor(ch.defining_ideal().ring()))
            [A2^2 + 2*A2*A3 + A3^2, A1*A2, A1*A2 + A2*A3, A1*A2 + A1*A3, A0*A2,
             A0*A2 + A2*A3, A0*A2 + A0*A3, A1^2 + 2*A1*A3 + A3^2, A0*A1,
             A0*A1 + A1*A3, A0*A1 + A0*A3, A0^2 + 2*A0*A3 + A3^2]
        """
        E = list(self._matroid.groundset())
        Q = []  # Quadratic generators
        flats_containing = {x: [] for x in E}
        for F in self._flats:
            for x in F:
                flats_containing[x].append(F)
        reln = lambda x,y: x <= y
        lattice_flats = Poset((self._flats, reln))
        antichains = lattice_flats.antichains().elements_of_depth_iterator(2)
        for F, G in antichains:
            Q.append(self._flats_generator[F] * self._flats_generator[G])
        for F in self._flats:
            for x in E:  # generators for every set of flats containing element
                term = poly_ring.zero()
                for H in flats_containing[x]:
                    term += self._flats_generator[H]
                if term**2 not in Q:
                    Q.append(term**2)

                if F not in flats_containing[x]:  # generators for every set of flats not containing element
                    Q.append(self._flats_generator[F]*term)
        return Q

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
            sage: ch.defining_ideal()
            Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on
            6 elements with 16 bases in the atom-free presentation
        """
        return "Augmented Chow ring ideal of {} in the atom-free presentation".format(self._matroid)

    def _latex_(self):
        r"""
        Return the LaTeX output of the ring and generators of ``self``.

        EXAMPLES::

            sage: M1 = Matroid(graphs.CycleGraph(3))
            sage: ch = M1.chow_ring(QQ, True, 'atom-free')
            sage: ch.defining_ideal()._latex_()
            'I_{af}(\\text{\\texttt{Graphic{ }matroid{ }of{ }rank{ }2{ }on{ }3{ }elements}})'
        """
        from sage.misc.latex import latex
        return 'I_{{af}}({})'.format(latex(self._matroid))

    def groebner_basis(self, algorithm='', *args, **kwargs):
        """
        Return the Groebner basis of ``self``.

        EXAMPLES::

            sage: M1 = matroids.Uniform(3, 6)
            sage: ch = M1.chow_ring(QQ, True, 'atom-free')
            sage: ch.defining_ideal().groebner_basis(algorithm='')
            Polynomial Sequence with 253 Polynomials in 22 Variables
            sage: ch.defining_ideal().groebner_basis(algorithm='').is_groebner()
            True
            sage: ch.defining_ideal().hilbert_series() == ch.defining_ideal().gens().ideal().hilbert_series()
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().groebner_basis(algorithm=algorithm, *args, **kwargs)
        gb = []
        poly_ring = self.ring()
        lattice_flats = Poset((self._flats, lambda x, y: x <= y))
        antichains = lattice_flats.antichains().elements_of_depth_iterator(2)
        for F, G in antichains:
            gb.append(self._flats_generator[F]*self._flats_generator[G])
        for F in self._flats:
            term = poly_ring.zero()
            for H in lattice_flats.order_filter([F]):
                term += self._flats_generator[H]
            if term != poly_ring.zero():
                order_ideal_modified = lattice_flats.order_ideal([F])
                order_ideal_modified.remove(F)
                for G in order_ideal_modified:
                    gb.append(self._flats_generator[G] * (term ** (self._matroid.rank(F) - self._matroid.rank(G))))
            gb.append(term ** (self._matroid.rank(F) + 1))

        return PolynomialSequence(poly_ring, [gb])

    def normal_basis(self, algorithm='', *args, **kwargs):
        r"""
        Return the monomial basis of the quotient ring of this ideal.

        EXAMPLES::

            sage: ch = Matroid(graphs.CycleGraph(3)).chow_ring(QQ, True, 'atom-free')
            sage: I = ch.defining_ideal()
            sage: I.normal_basis()
            [1, A0, A1, A2, A3, A3^2]
            sage: set(I.gens().ideal().normal_basis()) == set(I.normal_basis())
            True
            sage: ch = matroids.Wheel(3).chow_ring(QQ, True, 'atom-free')
            sage: I = ch.defining_ideal()
            sage: I.normal_basis()
            Polynomial Sequence with 30 Polynomials in 14 Variables
            sage: set(I.gens().ideal().normal_basis()) == set(I.normal_basis())
            True
        """
        if algorithm == '':
            algorithm = 'constructed'
        if algorithm != 'constructed':
            return super().normal_basis(algorithm=algorithm, *args, **kwargs)
        R = self.ring()
        flats_gen = self._flats_generator
        monomial_basis = []
        ranks, chains = self._lattice_flats()
        for subset in chains:
            max_powers = []
            k = len(subset)
            if not subset:
                monomial_basis.append(R.one())
            else:
                for i in range(k):
                    if i == 0:
                        max_powers.append(ranks[subset[i]])
                    else:
                        max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                ranges = [range(1, p) for p in max_powers]
                ranges[0] = range(1, max_powers[0] + 1)
                first_rank = ranks[subset[k-1]] + 1
                for combination in product(*(r for r in ranges)):
                    # Generating combinations for all powers from 1 to max_powers
                    if sum(combination) <= first_rank:
                        expression = R.one()
                        for val, c in zip(subset, combination):
                            expression *= flats_gen[val] ** c
                        monomial_basis.append(expression)
        return PolynomialSequence(R, [monomial_basis])
