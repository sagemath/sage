r"""
Chow rings of matroids

AUTHORS:

- Shriya M
"""

from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug, AugmentedChowRingIdeal_fy, AugmentedChowRingIdeal_atom_free
from sage.rings.quotient_ring import QuotientRing_generic
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.commutative_rings import CommutativeRings
from sage.misc.misc_c import prod
from itertools import product

class ChowRing(QuotientRing_generic):
    r"""
    The class of Chow ring, a multi-polynomial quotient ring.

    INPUT:


    - `M` -- a matroid
    - `R` -- a commutative ring
    - ``augmented`` -- a Boolean value. When ``True``, it returns the augmented
        Chow ring. If ``False``, it returns the Chow ring
    - ``presentation`` -- a string literal. Takes in the presentation of
        augmented Chow ring (`fy` or `atom-free`). Default value `None`

    OUTPUT: Chow ring of matroid `M`

    These are the classes of Chow rings for matroids. It also takes in
    a parameter boolean ``augmented`` which creates the augmented Chow
    ring if given ``True``.

    The Chow ring of a matroid is defined as the quotient ring:

    .. MATH::

        A^*(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / (Q_M + L_M)

        where

    ..MATH::

        (Q_M + L_M)

        is the Chow ring ideal of matroid `M`.

    The augmented Chow ring of matroid `M` of Feitchner-Yuzvinsky presentation
    is the quotient ring:

    ..MATH::

        A(M)_R := R[y_{e_1}, \ldots, y_{e_n}, x_{F_1}, \ldots, x_{F_k}] / (I_M + J_M)

        where

    ..MATH::

        (I_M + J_M)

        is the augmented Chow ring ideal of matroid `M` of Feitchner-Yuzvinsky
        presentation.

    The augmented Chow ring of atom-free presentation is the quotient ring:

    ..MATH::

        A(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / I_{af}M

        where

        I_{af}M
        
        is the augmented Chow ring ideal of matroid `M` of atom-free presentation.

    .. SEEALSO::

        :mod: sage.matroids.chow_ring_ideal

    REFERENCES:

        - [FY2004]_
        - [AHK2015]_

    EXAMPLES::

        sage: M1 = matroids.catalog.P8pp()
        sage: ch = M1.chow_ring(QQ, False)
        sage: ch
        Chow ring of P8'': Matroid of rank 4 on 8 elements with 8 nonspanning circuits
    """
    def __init__(self, R, M, augmented, presentation=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: I = matroids.Wheel(3).chow_ring(QQ, False)
            sage: TestSuite(I).run()
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
            self._ideal = ChowRingIdeal_nonaug(M, R) #check method to get ring
        C = CommutativeRings().Quotients() & GradedAlgebrasWithBasis(R).FiniteDimensional()
        QuotientRing_generic.__init__(self, R=self._ideal.ring(), I=self._ideal, names=self._ideal.ring().variable_names(), category=C)

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
                return "Augmented Chow ring of {} of Feitchner-Yuzvinsky presentation".format(self._matroid)
            elif self._presentation == 'atom-free':
                return "Augmented Chow ring of {} of atom-free presentation".format(self._matroid)
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        r"""
        Return the LaTeX output of the polynomial ring and Chow ring ideal.

        EXAMPLES::

            sage: from sage.matroids.graphic_matroid import GraphicMatroid
            
            sage: M1 = GraphicMatroid(graphs.CycleGraph(3))
            sage: ch = M1.chow_ring(QQ, True, 'fy')
            sage: ch._latex_()
            \Bold{Q}[A_{0}, A_{1}, A_{2}, B_{0}, B_{1}, B_{2}]/\left(B_{0}^{2},
            B_{0} B_{1}, B_{0} B_{2}, B_{0} B_{1}, B_{1}^{2}, B_{1} B_{2},
            B_{0} B_{2}, B_{1} B_{2}, B_{2}^{2}, A_{0} - B_{1} - B_{2},
            A_{1} - B_{0} - B_{2}, A_{2} - B_{0} - B_{1}\right)
            \Bold{Q}[A_{0}, A_{1}, A_{2}, B_{0}, B_{1}, B_{2}]

        """
        import sage.misc.latex as latex
        return "%s/%s" % (latex.latex(self._ideal.ring()), latex.latex(self._ideal))

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
            sage: ch = matroids.Wheel(3).chow_ring(ZZ, False)
            [A0*A013, 1]
        """
        flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        flats.append(frozenset())
        maximum_rank = max(self._matroid.rank(F) for F in flats)
        flats_gen = self._ideal.flats_generator()
        flats = sorted(flats, key=lambda X: (len(X), sorted(X)))
        ranks = [self._matroid.rank(F) for F in flats]
        monomial_basis = []
        if self._augmented is True:
            if self._presentation == 'fy':
                for i in range(maximum_rank):
                    term = self._ideal.ring().one()
                    for j in range(len(flats)):
                        if j == 0:
                            if i <= ranks[0]:
                                if flats[j] in flats_gen:
                                    term *= flats_gen[flats[j]]**(i + 1)
                        else:
                            if i < ranks[j] - ranks[j-1]:
                                if flats[j] in flats_gen:
                                    term *= flats_gen[flats[j]]**(i + 1)
                    monomial_basis.append(term)
                
            elif self._presentation == 'atom-free': #all double equals need spacing
                first_rank = self._matroid.rank(flats[len(flats) - 1])
                print(first_rank)
                for i in range(maximum_rank):
                    pow = []
                    for j in range(1, len(flats) - 1):
                        if i < ranks[j] - ranks[j-1]:
                            pow.append((j , i))
                    if sum(p[1] for p in pow) == first_rank + 1:
                        term = prod(flats_gen[flats[p[0]]] ** p[1] for p in pow) 
                        monomial_basis.append(term)

        else:
            print(ranks)
            rank_diff = []
            for i in range(len(flats)):
                max_pow = 0
                for j in range(i):
                    if flats[j] < flats[i]:
                        max_pow += ranks[j]
                rank_diff.append(ranks[i] - max_pow)
            print(rank_diff)
            def generate_terms(current_term, index, max_powers, f_dict):
                if index == len(max_powers):
                    term = self._ideal.ring().one()
                    term *= f_dict[i+1]**(current_term[i] for i in range(len(current_term)))
                    yield term
                    return
                
                for power in range(max_powers[index]):
                    current_term[index] = power
                    generate_terms(current_term, index + 1, max_powers, f_dict)

            current_term = [0]*len(flats)
            monomial_basis = list(generate_terms(current_term, 0, rank_diff, flats_gen))
            print(monomial_basis)
            

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
                sage: v.to_vector() #Error in output!

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
                0
                sage: v.monomial_coefficients() #error in output!

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
                A0*A01 2
                1 0
                sage: v = sum(ch.basis())
                sage: v.degree()
                2
            """
            return self.lift().degree()

        def homogeneous_degree(self):
            r"""
            Return the (homogeneous) degree of ``self`` if homogeneous
            otherwise raise an error.

            EXAMPLES::

                ch = matroids.catalog.Fano().chow_ring(QQ, True, 'fy')
                sage: for b in ch.basis():
                ....:     print(b, b.homogeneous_degree())
                Ba*Babf 2
                Ba^2 2
                sage: v = sum(ch.basis()); v
                Ba^2 + Ba*Babf
                sage: v.homogeneous_degree() #error - basis() type is wrong!
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            TESTS::
              
                sage: ch = matroids.Wheel(3).chow_ring(QQ, False)
                sage: ch.zero().homogeneous_degree() #error!
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