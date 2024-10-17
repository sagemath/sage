r"""
Chow rings of matroids

AUTHORS:

- Shriya M
"""

from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug, AugmentedChowRingIdeal_fy, AugmentedChowRingIdeal_atom_free
from sage.rings.quotient_ring import QuotientRing_generic
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.commutative_rings import CommutativeRings
from itertools import product, combinations
from sage.combinat.posets.posets import Poset

class ChowRing(QuotientRing_generic):
    r"""
    The Chow ring of a matroid.

    The *Chow ring of a matroid* `M` is defined as the quotient ring

    .. MATH::

        A^*(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / (Q_M + L_M),

    where `(Q_M + L_M)` is the :class:`Chow ring ideal
    <sage.matroids.chow_ring_ideal.ChowRingIdeal_nonaug>` of matroid `M`.

    The *augmented Chow ring of matroid* `M` in the Feitchner-Yuzvinsky presentation
    is the quotient ring

    .. MATH::

        A(M)_R := R[y_{e_1}, \ldots, y_{e_n}, x_{F_1}, \ldots, x_{F_k}] / (I_M + J_M),

    where `(I_M + J_M)` is the :class:`augmented Chow ring ideal
    <sage.matroids.chow_ring_ideal.AugmentedChowRingIdeal_fy>` of matroid `M`
    in Feitchner-Yuzvinsky presentation.

    The *augmented Chow ring of atom-free presentation* is the quotient ring

    .. MATH::

        A(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / I_M^{af},

    where `I_M^{af}` is the :class:`augmented Chow ring ideal
    <sage.matroids.chow_ring_ideal.AugmentedChowRingIdeal_atom_free>`
    of matroid `M` in the atom-free presentation.

    .. SEEALSO::

        :mod:`sage.matroids.chow_ring_ideal

    INPUT:

    - ``M`` -- matroid
    - ``R`` -- commutative ring
    - ``augmented`` -- boolean; when ``True``, this is the augmented
      Chow ring and if ``False``, this is the non-augmented Chow ring
    - ``presentation`` -- string (default: ``None``); one of the following:

      * ``"fy"`` - the Feitchner-Yuzvinsky presentation
      * ``"atom-free"`` - the atom-free presentation

    REFERENCES:

    - [FY2004]_
    - [AHK2015]_

    EXAMPLES::

        sage: M1 = matroids.catalog.P8pp() #more examples
        sage: ch = M1.chow_ring(QQ, False)
        sage: ch
        Chow ring of P8'': Matroid of rank 4 on 8 elements with 8 nonspanning circuits
    """
    def __init__(self, R, M, augmented, presentation=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, False)
            sage: TestSuite(ch).run()
        """
        self._matroid = M
        self._augmented = augmented
        self._presentation = presentation
        if augmented is True:
            if presentation == 'fy':
                self._ideal = AugmentedChowRingIdeal_fy(M, R)
            elif presentation == 'atom-free':
                self._ideal = AugmentedChowRingIdeal_atom_free(M, R)
        else:
            self._ideal = ChowRingIdeal_nonaug(M, R)
        C = CommutativeRings().Quotients() & GradedAlgebrasWithBasis(R).FiniteDimensional()
        QuotientRing_generic.__init__(self, R=self._ideal.ring(),
                                      I=self._ideal,
                                      names=self._ideal.ring().variable_names(),
                                      category=C)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: M1 = matroids.catalog.Fano()
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch
            Chow ring of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
        """
        if self._augmented is True:
            if self._presentation == 'fy':
                return "Augmented Chow ring of {} in Feitchner-Yuzvinsky presentation".format(self._matroid)
            elif self._presentation == 'atom-free':
                return "Augmented Chow ring of {} in atom-free presentation".format(self._matroid)
        return "Chow ring of {}".format(self._matroid)

    def _latex_(self):
        r"""
        Return the LaTeX output of the polynomial ring and Chow ring ideal.

        EXAMPLES::

            sage: from sage.matroids.graphic_matroid import GraphicMatroid

            sage: M1 = matroids.Uniform(2,5)
            sage: ch = M1.chow_ring(QQ, False)
            sage: ch._latex_()
            '\\Bold{Q}[A_{0}, A_{1}, A_{2}, A_{3}, A_{4}, A_{01234}] / \\left(A_{0} A_{1}, A_{0} A_{2}, A_{0} A_{3}, A_{0} A_{4}, A_{1} A_{2}, A_{1} A_{3}, A_{1} A_{4}, A_{2} A_{3}, A_{2} A_{4}, A_{3} A_{4}, A_{0} + A_{01234}, A_{1} + A_{01234}, A_{2} + A_{01234}, A_{3} + A_{01234}, A_{4} + A_{01234}\\right)\\Bold{Q}[A_{0}, A_{1}, A_{2}, A_{3}, A_{4}, A_{01234}]'
        """
        from sage.misc.latex import latex
        return "{} / {}".format(latex(self._ideal.ring()), latex(self._ideal))

    def _coerce_map_from_base_ring(self):
        r"""
        Disable the coercion from the base ring from the category.

        TESTS::

            sage: ch = matroids.Wheel(3).chow_ring(QQ, False)
            sage: ch._coerce_map_from_base_ring() is None
            True
        """
        return None  # don't need anything special

    def basis(self):
        r"""
        Return the monomial basis of the given Chow ring.

        EXAMPLES::

            sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, True, 'fy')
            sage: ch.basis()
            [B0*B01, 1]
            sage: len(ch.defining_ideal().normal_basis()) == ch.basis().cardinality()
            True
            sage: ch = matroids.catalog.NonFano().chow_ring(QQ, False)
            [A0*A013, 1]
            sage: len(ch.defining_ideal().normal_basis()) == ch.basis().cardinality()
            True

        """
        flats = [X for i in range(1, self._matroid.rank() + 1)
                 for X in self._matroid.flats(i)] #Non empty flats
        flats_gen = self._ideal.flats_generator()
        R = self._ideal.ring()
        flats = sorted(flats, key=lambda X: (len(X), sorted(X)))
        ranks = {F: self._matroid.rank(F) for F in flats}
        monomial_basis = []
        reln = lambda x,y: x <= y
        lattice_flats = Poset((flats, reln))
        chains = lattice_flats.chains() #Only chains
        if self._augmented:
            if self._presentation == 'fy':
                for subset in chains:
                    k = len(subset)
                    if k == 0:
                        monomial_basis.append(R.one())
                    else:
                        max_powers = []
                        max_powers.append(ranks[subset[0]])
                        for i in range(1, k):
                            max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                        ranges = [range(1, p) for p in max_powers]
                        ranges[0] = range(1, max_powers[0] + 1)
                        for combination in product(*(r for r in ranges)):
                            #Generating combinations for all powers up to max_powers
                            expression = R.one()
                            for i in range(k):
                                expression *= flats_gen[subset[i]]**combination[i]
                            monomial_basis.append(expression)

            elif self._presentation == 'atom-free':
                for subset in chains:
                    max_powers = []
                    k = len(subset)
                    if subset == []:
                        monomial_basis.append(R.one())
                    else:
                        for i in range(k):
                            if i == 0:
                                max_powers.append(ranks[subset[i]])
                            else:
                                max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                        first_rank = ranks[subset[k-1]] + 1
                        last_rank = ranks[subset[0]]
                        for combination in product(*(range(1, p) for p in max_powers)):
                            #Generating combinations for all powers from 1 to max_powers
                            if sum(combination) == first_rank:
                                expression = R.one()
                                for i in range(k):
                                    expression *= flats_gen[subset[i]]**combination[i]
                                monomial_basis.append(expression)
                        max_powers[0] = 0
                        for combination in product(*(range(1, p) for p in max_powers)):
                            #Generating all combinations including 0 power and max_power for first flat
                            if sum(combination) <= first_rank:
                                expression = R.one()
                                for i in range(len(combination)):
                                    expression *= flats_gen[subset[i]]**combination[i]
                                monomial_basis.append(expression)
                            else:
                                expression *= flats_gen[subset[0]]**last_rank
                                if sum(combination) + last_rank == first_rank:
                                    for i in range(k):
                                        expression *= flats_gen[subset[i]]**combination[i]
                                monomial_basis.append(expression)

        else:
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
                        for i in range(k):
                            expression *= flats_gen[subset[i]]**combination[i]
                        monomial_basis.append(expression)

        from sage.sets.family import Family
        return Family([self.element_class(self, mon, reduce=False) for mon in monomial_basis])

    class Element(QuotientRing_generic.Element):
        def to_vector(self, order=None):
            r"""
            Return ``self`` as a (dense) free module vector.

            EXAMPLES::

                sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False)
                sage: v = ch.an_element(); v
                A0
                sage: v.to_vector()
                (0, 0, 0, 1, 0, 0, 0)
            """
            P = self.parent()
            B = P.basis()
            FM = P._dense_free_module()
            f = self.lift()
            return FM([f.monomial_coefficient(b.lift()) for b in B])

        _vector_ = to_vector

        def monomial_coefficients(self, copy=None):
            r"""
            Return the monomial coefficients of ``self``.

            EXAMPLES::

                sage: ch = matroids.catalog.NonFano().chow_ring(QQ, True, 'fy')
                sage: v = ch.an_element(); v
                Ag
                sage: v.monomial_coefficients()
                {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}
            """
            B = self.parent().basis()
            f = self.lift()
            return {i: f.monomial_coefficient(b.lift()) for i, b in enumerate(B)}

        def degree(self):
            r"""
            Return the degree of ``self``.

            EXAMPLES::

                sage: ch = matroids.Uniform(3, 6).chow_ring(QQ, False)
                sage: for b in ch.basis():
                ....:     print(b, b.degree())
                1 0
                A1 1
                A3 1
                A0 1
                A2 1
                A4 1
                A5 1
                sage: v = sum(ch.basis())
                sage: v.degree()
                1
            """
            return self.lift().degree()

        def homogeneous_degree(self):
            r"""
            Return the (homogeneous) degree of ``self`` if homogeneous
            otherwise raise an error.

            EXAMPLES::

                sage: ch = matroids.catalog.Fano().chow_ring(QQ, True, 'fy')
                sage: for b in ch.basis():
                ....:     print(b, b.homogeneous_degree())
                1 0
                Bb 1
                Ba 1
                Bc 1
                Bd 1
                Be 1
                Bf 1
                Bg 1
                sage: v = sum(ch.basis()); v
                Ba + Bb + Bc + Bd + Be + Bf + Bg + 1
                sage: v.homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            TESTS::

                sage: ch = matroids.Wheel(3).chow_ring(QQ, False)
                sage: ch.zero().homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if not self:
                raise ValueError("the zero element does not have a well-defined degree")
            f = self.lift()
            if not f.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return f.degree()