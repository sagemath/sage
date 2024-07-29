r"""
Chow rings of matroids

AUTHORS:

- Shriya M

These are the classes of Chow rings for matroids. It also takes in
a parameter boolean ``augmented`` which creates the augmented Chow
ring if given ``True``.


REFERENCES

- :arxiv:`2309.14312`
- :arxiv:`2111.00393`
"""

from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug, AugmentedChowRingIdeal_fy, AugmentedChowRingIdeal_atom_free
from sage.rings.quotient_ring import QuotientRing_nc
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.sets.set import Set
from sage.combinat.posets.posets import Poset


class ChowRing(QuotientRing_nc):
    r"""
    The class of Chow ring, a multi-polynomial quotient ring.
    Base class - ``QuotientRing_nc``.

    INPUT:


    - `M` -- a matroid.
    - `R` -- a ring.
    - ``augmented`` -- a Boolean value. When ``True``, it returns the augmented Chow
      ring. If ``False``, it returns the Chow ring

    OUTPUT: Chow ring of matroid `M`.

    EXAMPLES::

        sage: from sage.matroids.chow_ring import ChowRing

        sage: M1 = matroids.catalog.P8pp()
        sage: ch = ChowRing(M=M1, R=QQ, augmented=False)
        sage: ch
        Chow ring of P8'': Matroid of rank 4 on 8 elements with 8 nonspanning circuits
    """
    def __init__(self, R, M, augmented, presentation=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring import ChowRing

            sage: I = ChowRing(M=matroids.Wheel(3), R=QQ, augmented=False)
            sage: TestSuite(I).run()
        """
        self._matroid = M
        self._augmented = augmented
        self._presentation = presentation
        if augmented:
            if presentation=='fy':
                self._ideal = AugmentedChowRingIdeal_fy(M, R)
            elif presentation=='atom-free':
                self._ideal = AugmentedChowRingIdeal_atom_free(M, R)
        else:
            self._ideal = ChowRingIdeal_nonaug(M, R) #check method to get ring
        QuotientRing_nc.__init__(self, R=self._ideal.ring(), I=self._ideal, names=self._ideal.ring().variable_names(), category=GradedAlgebrasWithBasis(R))

    def _repr_(self):
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        import sage.misc.latex as latex
        return "%s/%s" % (latex.latex(self.poly_ring), latex.latex(self._ideal))

    def _coerce_map_from_base_ring(self):
        r"""
        Disable the coercion from the base ring from the category.

        TESTS::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: GP = SGA.garsia_procesi_module([2, 2])
            sage: GP._coerce_map_from_base_ring() is None
            True
        """
        return None  # don't need anything special

    def basis(self):
        flats = [Set(X) for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        flats.append(Set([]))
        maximum_rank = max([self._matroid.rank(F) for F in flats])
        func = lambda A,B: A.issubset(B)
        flats = Poset((flats, func), cover_relations=True).directed_subsets('down') #DEBUG
        monomial_basis = []
        if self._augmented:
            if self._presentation=='fy':
                for i in range(maximum_rank):
                    term = self._ideal.ring().one()
                    for j in range(len(flats)):
                        if j == 0:
                            if i > self._matroid.rank(flats[0]):
                                term *= self._flats_generator[flats[j]]**(0)
                            else:
                                term *= self._flats_generator[flats[j]]**(i + 1)
                        else:
                            if i >= (self._matroid.rank(flats[j]) - self._matroid.rank(flats[j-1])):
                                term *= self._flats_generator[flats[j]]**(0)
                            else:
                                term *= self._flats_generator[flats[j]]**(i + 1)
                    monomial_basis.append(term)
                
            elif self._presentation=='atom-free':
                for i in range(maximum_rank):
                    term = self._ideal.ring().one()
                    pow = []
                    for j in range(1, len(flats) - 1):
                        if i >= (self._matroid.rank(flats[j-1]) - self._matroid.rank(flats[j])):
                            pow.append((j , i))
                    if sum(pow) == self._matroid.rank(flats[0]):
                        for p in pow:
                            term *= self._flats_generator[flats[p[0]]]**(p[1])
                    monomial_basis.append(term)
        

        else:
            for i in range(maximum_rank):
                    term = self._ideal.ring().one()
                    for j in range(len(flats)):
                        if i > (self._matroid.rank(flats[j]) - self._matroid.rank(flats[j-1]) - 1):
                                term *= self._flats_generator[flats[j]]**(0)
                        else:
                            term *= self._flats_generator[flats[j]]**(i + 1)
                    monomial_basis.append(term)


        m_basis = PolynomialSequence([monomial_basis], self._ideal.ring())
        return m_basis
                        
                



                    
                        







    class Element(QuotientRing_nc.Element):
        def to_vector(self, order=None):
            r"""
            Return ``self`` as a (dense) free module vector.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP22 = SGA.garsia_procesi_module([2, 2])
                sage: v = GP22.an_element(); v
                -gp1 - gp2 - gp3
                sage: v.to_vector()
                (0, 0, 2, 2, 2, 0)
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

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP31 = SGA.garsia_procesi_module([3, 1])
                sage: v = GP31.an_element(); v
                -gp1 - gp2 - gp3
                sage: v.monomial_coefficients()
                {0: 2, 1: 2, 2: 2, 3: 0}
            """
            B = self.parent().basis()
            f = self.lift()
            return {i: f.monomial_coefficient(b.lift()) for i, b in enumerate(B)}

        def degree(self):
            r"""
            Return the degree of ``self``.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP22 = SGA.garsia_procesi_module([2, 2])
                sage: for b in GP22.basis():
                ....:     print(b, b.degree())
                gp2*gp3 2
                gp1*gp3 2
                gp3 1
                gp2 1
                gp1 1
                1 0
                sage: v = sum(GP22.basis())
                sage: v.degree()
                2
            """
            return self.lift().degree()

        def homogeneous_degree(self):
            r"""
            Return the (homogeneous) degree of ``self`` if homogeneous
            otherwise raise an error.

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(GF(2), 4)
                sage: GP31 = SGA.garsia_procesi_module([3, 1])
                sage: for b in GP31.basis():
                ....:     print(b, b.homogeneous_degree())
                gp3 1
                gp2 1
                gp1 1
                1 0
                sage: v = sum(GP31.basis()); v
                gp1 + gp2 + gp3 + 1
                sage: v.homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            TESTS::

                sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
                sage: GP4 = SGA.garsia_procesi_module([4])
                sage: GP4.zero().homogeneous_degree()
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