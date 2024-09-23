r"""
Chow rings of matroids

AUTHORS:

- Shriya M
"""

from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug, AugmentedChowRingIdeal_fy, AugmentedChowRingIdeal_atom_free
from sage.rings.quotient_ring import QuotientRing_generic
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.commutative_rings import CommutativeRings
from sage.combinat.posets.posets import Poset
from itertools import product, combinations

class ChowRing(QuotientRing_generic):
    r"""
    The Chow ring of a matroid.

    The *Chow ring of a matroid* `M` is defined as the quotient ring:

    .. MATH::

        A^*(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / (Q_M + L_M)

    where `(Q_M + L_M)` is the Chow ring ideal of matroid `M`.

    The *augmented Chow ring of matroid* `M` in the Feitchner-Yuzvinsky presentation
    is the quotient ring:

    ..MATH::

        A(M)_R := R[y_{e_1}, \ldots, y_{e_n}, x_{F_1}, \ldots, x_{F_k}] / (I_M + J_M)

    where `(I_M + J_M)` is the augmented Chow ring ideal of matroid `M`
    of Feitchner-Yuzvinsky presentation.

    The *augmented Chow ring of atom-free presentation* is the quotient ring:

    ..MATH::

        A(M)_R := R[x_{F_1}, \ldots, x_{F_k}] / I_{af}M

    where `I_{af}M` is the augmented Chow ring ideal of matroid `M` of atom-free presentation.

    .. SEEALSO::

        :mod:`sage.matroids.chow_ring_ideal`
    
    INPUT:

    - `M` -- a matroid
    - `R` -- a commutative ring
    - ``augmented`` -- boolean; when ``True``, this is the augmented
        Chow ring and if ``False``, this is the non-augmented Chow ring
    - ``presentation`` -- string (default: ``None``); one of the following:
    
      * ``"fy"`` - Feitchner-Yuzvinsky presentation*
      * ``"atom-free" - Atom-free presentation*


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
        flats_gen = self._ideal.flats_generator()
        R = self._ideal.ring()
        flats = sorted(flats, key=lambda X: (len(X), sorted(X)))
        ranks = {F:self._matroid.rank(F) for F in flats}
        monomial_basis = []
        if self._augmented is True:
            if self._presentation == 'fy':
                flats.remove(frozenset()) #Non empty proper flats
                max_powers = []
                max_powers[0] = ranks[flats[0]]
                for i in range(1, len(flats)):
                    max_powers = ranks[flats[i]] - ranks[flats[i-1]]
                for combination in product(*(range(p) for p in max_powers)): #Generating combinations for all powers up to max_powers
                        expression = R.one()
                        for i in range(len(flats)):
                            expression *= flats_gen[subset[i]]**combination[i] 
                        monomial_basis.append(expression)
                        if combination[0] == 0: #Generating combinations for all powers including first max_powers
                            expression *= flats_gen[subset[0]]**max_powers[0]
                            monomial_basis.append(expression)
                
            elif self._presentation == 'atom-free':
                subsets = []
        # Generate all subsets of the frozenset using combinations
            for r in range(len(flats) + 1):  # r is the size of the subset
                subsets.extend(list(subset) for subset in combinations(flats, r))
            for subset in subsets:
                flag = True
                sorted_list = sorted(subset, key=len)
                for i in range (len(sorted_list)): #Taking only chains
                    if (i != 0) & (len(sorted_list[i]) == len(sorted_list[i-1])):
                        flag = False
                        break
                if flag is True: #For every chain
                    max_powers = []
                    x_dict = dict()
                    k = len(subset)
                    for i in range(k-1):
                        if i == 0:
                            max_powers.append(ranks[subset[i]])
                            x_dict[subset[i]] = flats_gen[subset[i]]
                        else:
                            max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                            x_dict[subset[i]] = flats_gen[subset[i]]
                    x_dict[subset[k-1]] = flats_gen[subset[k-1]]
                    max_powers[k-1] = ranks[subset[k-1]]
                    first_rank = ranks[subset[0]] + 1
                    last_rank = ranks[subset[k-1]]
                    for combination in product(*(range(1, p) for p in max_powers)): #Generating combinations for all powers from 1 to max_powers
                        expression = R.one()
                        if sum(combination) == first_rank:
                            for i in range(k):
                                expression *= x_dict[subset[i]]**combination[i] 
                            monomial_basis.append(expression)
                    max_powers.remove(last_rank)
                    for combination in product(*(range(1, p) for p in max_powers)): #Generating all combinations including 0 power and max_power for first flat
                        expression = R.one()
                        if sum(combination) == first_rank:
                            for i in range(len(combination)):
                                expression *= x_dict[subset[i]]**combination[i] 
                            monomial_basis.append(expression)
                        else:
                            expression *= x_dict[subset[k-1]]**last_rank
                            if sum(combination) + last_rank == first_rank:
                                for i in range(k):
                                    expression *= x_dict[subset[i]]**combination[i] 
                            monomial_basis.append(expression)

        else:
            flats.remove(frozenset()) #Non empty proper flats
            subsets = []
        # Generate all subsets of the frozenset using combinations
            for r in range(len(flats) + 1):  # r is the size of the subset
                subsets.extend(list(subset) for subset in combinations(flats, r))
            for subset in subsets:
                flag = True
                sorted_list = sorted(subset, key=len)
                for i in range (len(sorted_list)):
                    if (i != 0) & (len(sorted_list[i]) == len(sorted_list[i-1])):
                        flag = False
                        break

                if flag is True:
                    max_powers = []
                    x_dict = dict()
                    for i in range(len(subset)):
                        if i == 0:
                            max_powers.append(ranks[subset[i]])
                            x_dict[subset[i]] = flats_gen[subset[i]]
                        else:
                            max_powers.append(ranks[subset[i]] - ranks[subset[i-1]])
                            x_dict[subset[i]] = flats_gen[subset[i]]
                    k = len(subset)
                    for combination in product(*(range(p) for p in max_powers)):
                        expression = R.one()
                        for i in range(k):
                            expression *= x_dict[subset[i]]**combination[i] 
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
                sage: v.monomial_coefficients()

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