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
from functools import cmp_to_key
from sage.misc.misc_c import prod

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
        r"""
        EXAMPLE::

        sage: from sage.matroids.chow_ring import ChowRing

        sage: M1 = matroids.catalog.Fano()
        sage: ch = ChowRing(M=M1, R=QQ, augmented=False)
        sage: ch
        Chow ring of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
        """
        if self._augmented:
            if self._presentation=='fy':
                return "Augmented Chow ring of {} of Feitchner-Yuzvinsky presentation".format(self._matroid)
            elif self._presentation=='atom-free':
                return "Augmented Chow ring of {} of atom-free presentation".format(self._matroid)
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        r"""
        Return the LaTeX output of the polynomial ring and Chow ring ideal.

        EXAMPLE::

            sage: from sage.matroids.chow_ring import ChowRing
            sage: from sage.matroids.graphic_matroid import GraphicMatroid
            
            sage: M1 = GraphicMatroid(graphs.CycleGraph(3))
            sage: ch = ChowRing(M=M1, R=QQ, augmented=True, presentation='fy')
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

            sage: from sage.matroids.chow_ring import ChowRing

            sage: ch = ChowRing(M=matroids.Wheel(3), R=QQ, augmented=False)
            sage: ch._coerce_map_from_base_ring() is None
            True
        """
        return None  # don't need anything special

    def basis(self):
        r"""
        Return the monomial basis of the given Chow ring.
        
        EXAMPLES::

            sage: from sage.matroids.chow_ring import ChowRing

            sage: ch = ChowRing(M=matroids.Uniform(3, 6), R=QQ, augmented=True, presentation='fy')
            sage: ch.basis()
            [B0*B01, 1]
            sage: ch = ChowRing(M=matroids.Wheel(3), R=ZZ, augmented=False)
            [A0*A013, 1]
        """
        flats = [Set(X) for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        flats.append(Set([]))
        maximum_rank = max([self._matroid.rank(F) for F in flats])
        flats_gen = self._ideal.flats_generator()
        def func(A, B):
            if A.issubset(B):
                return -1
            elif B.issubset(A):
                return 1
            else:
                return 0
        flats = sorted(flats, key=cmp_to_key(func))
        ranks = [self._matroid.rank(F) for F in flats]
        monomial_basis = []
        if self._augmented:
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
            for i in range(maximum_rank):
                    term = self._ideal.ring().one()
                    for j in range(len(flats)):
                        if i > ranks[j] - ranks[j-1] - 1:
                                if flats[j] in list(flats_gen):
                                    term *= flats_gen[flats[j]]**(0)
                        else:
                            if flats[j] in list(flats_gen):
                                term *= flats_gen[flats[j]]**(i + 1)
                    monomial_basis.append(term)


        m_basis = PolynomialSequence(self._ideal.ring(), monomial_basis)
        return m_basis

    class Element(QuotientRing_nc.Element):
        def to_vector(self, order=None):
            r"""
            Return ``self`` as a (dense) free module vector.

            EXAMPLES::

                sage: ch = ChowRing(M=matroids.Uniform(3, 6), R=QQ, augmented=False)
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

                sage: ch = ChowRing(M=matroids.catalog.NonFano(), R=QQ, augmented=True, presentation='fy')
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

                sage: ch = ChowRing(M=matroids.Uniform(3, 6), R=QQ, augmented=False)
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

                ch = ChowRing(M=matroids.catalog.Fano(), R=QQ, augmented=True, presentation='fy')
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
              
                sage: ch = ChowRing(M=matroids.Wheel(3), R=QQ, augmented=False)
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