r"""
Chow ring ideals of matroids

AUTHORS:

- Shriya M

These are the classes of Chow ring ideals for matroids. There are three classes
created - Chow ring ideal and augmented Chow ring ideal. The augmented
Chow ring ideal has two different presentations implemented - the Feitchner-
Yuzvinsky presentation and atom-free presentation. Both classes have 
``grobner_basis()`` methods as well, as an explicit Groebner basis is known 
in each case.

REFERENCES

- :arxiv:`2309.14312`
- :arxiv:`2111.00393`
"""

from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.misc.abstract_method import abstract_method

class ChowRingIdeal(MPolynomialIdeal):
    @abstract_method
    def _gens_constructor():
        pass

    def matroid(self):
        r"""
        Return the matroid of the given Chow ring ideal.

        EXAMPLE::
            
            sage: ch = ChowRingIdeal_nonaug(M=matroids.Uniform(3,6), R=QQ)
            sage: ch.matroid()
            U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
            {3: {{0, 1, 2, 3, 4, 5}}}
        """
        M = self._matroid
        return M
    
    def flats_generator(self):
        r"""
        Return the variables of every corresponding flat/groundset element
        of the matroid.

        EXAMPLE::

            sage: ch = ChowRingIdeal_nonaug(M=matroids.catalog.Fano(), R=QQ)
            sage: ch.flats_generator() #WHERE IS OUTPUT?
        """
        return dict(self._flats_generator)
    

class ChowRingIdeal_nonaug(ChowRingIdeal): 
    r"""
    The Chow ring ideal. 

    INPUT:

    - `M` -- a matroid
    - `R` -- a ring

    OUTPUT: Chow ring ideal of matroid `M`

    EXAMPLES:

    Chow ring ideal of uniform matroid of rank 3 on 6 elements::

        sage: from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug

        sage: ch = ChowRingIdeal_nonaug(M=matroids.Uniform(3,6), R=QQ)
        sage: ch
        Chow ring ideal of U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
        {3: {{0, 1, 2, 3, 4, 5}}}
        sage: ch = ChowRingIdeal_nonaug(M=matroids.catalog.Fano(), R=QQ)
        sage: ch
        Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
    """
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug

            sage: I = ChowRingIdeal_nonaug(M=matroids.catalog.Fano(), R=QQ)
            sage: TestSuite(I).run(skip="_test_category")
        """
        self._matroid = M
        flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]
        try:
            poly_ring = PolynomialRing(R, names) #self.ring
        except ValueError: # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(self.flats))
        gens = poly_ring.gens()
        self._flats_generator = dict(zip(flats, gens))
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))

    def _gens_constructor(self, poly_ring):
        r"""
        Returns the generators of the Chow ring ideal.
        
        EXAMPLE::
            sage: ch = ChowRingIdeal_nonaug(M=matroids.catalog.NonFano(), R=QQ)
            sage: ch._gens_constructor() #WHERE IS OUTPUT?

        """
        E = list(self._matroid.groundset())
        flats = list(self._flats_generator.keys())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(flats):
            for x in F:
                flats_containing[x].append(i)
        gens = poly_ring.gens()
        Q = [gens[i] * gens[i+j+1] for i,F in enumerate(flats)
            for j,G in enumerate(flats[i+1:]) if not (F < G or G < F)]
        L = [sum(gens[i] for i in flats_containing[x])
            - sum(gens[i] for i in flats_containing[y])
            for j,x in enumerate(E) for y in E[j+1:]]
        return Q + L

    def _repr_(self):
        r"""
        EXAMPLE::

            sage: ch = ChowRingIdeal_nonaug(M=matroids.catalog.Fano(), R=QQ)
            sage: ch
            Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
        """
        return "Chow ring ideal of {}".format(self._matroid)

    def groebner_basis(self):
        r"""
        Returns the Groebner basis of the Chow ring ideal.
        Return type - ``PolynomialSequence``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring_ideal import ChowRingIdeal_nonaug
            sage: from sage.matroids.basis_matroid import BasisMatroid

            sage: ch = ChowRingIdeal_nonaug(M=BasisMatroid(groundset='abc', bases=['ab', 'ac']), R=QQ)
            sage: ch.groebner_basis()
            [Aa^2, Aa*Abc, Aa*Abc, Abc^2]
            sage: ch.groebner_basis().is_groebner()
            True
        
        Another example would be the Groebner basis of the Chow ring ideal of
        the Non-Fano matroid::

            sage: ch = ChowRingIdeal_nonaug(M=matroids.catalog.NonFano(), R=QQ)
            sage: ch.groebner_basis()
            Polynomial Sequence with 232 Polynomials in 16 Variables
        """

        flats = list(self._flats_generator)
        gb = list()
        R = self.ring()
        for F in flats:
            for G in flats: 
                if not (F < G or G < F):
                    gb.append(self._flats_generator[F]*self._flats_generator[G])
                else:
                    term = R.zero()
                    for H in flats:
                        if H > G:
                            term += self._flats_generator[H]
                    if Set(F).is_empty():
                        gb.append(term**self._matroid.rank(G)) 
                    elif F < G:
                        gb.append(term**(self._matroid.rank(G)-self._matroid.rank(F)))

        g_basis = PolynomialSequence(R, [gb])
        return g_basis
    
class AugmentedChowRingIdeal_fy(ChowRingIdeal):
    r"""
        The class of augmented Chow ring ideal of Feitchner-Yuzvinsky
        presentation, a multi-polynomial ideal.
        Base class - ``MPolynomialIdeal``.

        INPUT:


        - `M` -- a matroid.
        - `R` -- a ring.

        OUTPUT: augmented Chow ring ideal of matroid `M` of Feitchner-Yuzvinsky
        presentation.

        EXAMPLES::

        Augmented Chow ring ideal of Wheel matroid of rank 3::

            sage: from sage.matroids.chow_ring_ideal import AugmentedChowRingIdeal_fy

            sage: ch = AugmentedChowRingIdeal_fy(M=matroids.Wheel(3), R=QQ)
            sage: ch
            Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 
            6 elements with 16 bases of Feitchner-Yuzvinsky presentation
        """
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring_ideal import AugmentedChowRingIdeal_fy

            sage: I = AugmentedChowRingIdeal_fy(M=matroids.Wheel(3), R=QQ)
            sage: TestSuite(I).run()
        """
        self._matroid = M
        self._flats = [X for i in range(1, self._matroid.rank())
                for X in self._matroid.flats(i)]
        E = list(self._matroid.groundset())
        self._flats_generator = dict()
        try:
            names_groundset = ['A{}'.format(''.join(str(x))) for x in E]
            names_flats = ['B{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self._flats]
            poly_ring = PolynomialRing(R, names_groundset + names_flats)
        except ValueError:
            poly_ring = PolynomialRing(R, 'A', len(E) + len(self._flats))
        for i,x in enumerate(E):
            self._flats_generator[x] = poly_ring.gens()[i]
        for i,F in enumerate(self._flats):
            self._flats_generator[F] = poly_ring.gens()[len(E) + i]
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))
        
    
    def _gens_constructor(self, poly_ring):
        r"""
        Return the generators of augmented Chow ring ideal of
        Feitchner-Yuzvinsky presentation.

        EXAMPLES::

            sage: ch = AugmentedChowRingIdeal_fy(M=matroids.Wheel(3), R=QQ)
            sage: ch._gens_constructor() #WHERE IS OUTPUT?

        """
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for F in self._flats:
            for x in F:
                flats_containing[x].append(F)

        Q = list()
        for F in self._flats:
            for G in self._flats:
                    if not (F < G or G < F):
                        Q.append(self._flats_generator[F] * self._flats_generator[G])
        L = list()
        for x in E:
            term = poly_ring.zero()
            for F in self._flats:
                if F not in flats_containing[x]:
                    term += self._flats_generator[F]
            L.append(self._flats_generator[x] - term)
        return Q + L

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: ch = AugmentedChowRingIdeal_fy(M=matroids.Wheel(3), R=QQ)
            sage: ch
            Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 
            6 elements with 16 bases of Feitchner-Yuzvinsky presentation
        """
        return "Augmented Chow ring ideal of {} of Feitchner-Yuzvinsky presentation".format(self._matroid)
    
    def groebner_basis(self):
        r"""
        Returns the Groebner basis of the augmented Chow ring ideal.
        Return type - ``PolynomialSequence``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring_ideal import AugmentedChowRingIdeal_fy

            sage: ch = AugmentedChowRingIdeal_fy(M=matroids.catalog.Fano(), R=QQ)
            sage: ch.groebner_basis()
            Polynomial Sequence with 4116 Polynomials in 21 Variables
        """
        gb = []
        E = list(self._matroid.groundset())
        poly_ring = self.ring()
        for i in E:
            for F in self._flats:
                for G in self._flats:
                    term = poly_ring.zero()
                    term1 = poly_ring.zero()
                    for H in self._flats:
                        if i in Set(H):
                            term += self._flats_generator[H]
                        if H > F:
                            term1 += self._flats_generator[H]

                    gb.append(self._flats_generator[i] + term)
                    gb.append(term1**(self._matroid.rank(Set(F))) + 1)

                    if i in Set(F):
                        gb.append(self._flats_generator[i]*((term1)**self._matroid.rank(Set(F))))
            
                    elif not i in Set(F):
                        gb.append(self._flats_generator[i]*self._flats_generator[F])
                    
                    elif not (F < G or G < F):
                        gb.append(self._flats_generator[F]*self._flats_generator[G])
                    
                    elif G < F:
                        gb.append(self._flats_generator[G]*term1**(self._matroid.rank(Set(F))-self._matroid.rank(Set(G))))

        g_basis = PolynomialSequence(poly_ring, [gb])
        return g_basis

class AugmentedChowRingIdeal_atom_free(ChowRingIdeal):
    r"""
    The augmented Chow ring ideal in the atom-free
    presentation.

    INPUT:

    - ``M`` -- a matroid
    - ``R`` -- a ring

    OUTPUT: augmented Chow ring ideal of matroid `M` of atom-free
    presentation.

    EXAMPLES:

    Augmented Chow ring ideal of Wheel matroid of rank 3::

        sage: from sage.matroids.chow_ring_ideal import AugmentedChowRingIdeal_atom_free

        sage: ch = AugmentedChowRingIdeal_atom_free(M=matroids.Wheel(3), R=QQ)
        sage: ch
        Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 
        6 elements with 16 bases of atom-free presentation
    """ 
    def __init__(self, M, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring_ideal import AugmentedChowRingIdeal_atom_free

            sage: I = AugmentedChowRingIdeal_atom_free(M=matroids.Wheel(3), R=QQ)
            sage: TestSuite(I).run(skip="_test_category")
        """
        self._matroid = M
        self._flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(self._flats):
            for x in F:
                flats_containing[x].append(i)

        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self._flats]

        try:
            poly_ring = PolynomialRing(R, names) #self.ring
        except ValueError: # variables are not proper names
            poly_ring = PolynomialRing(R, 'A', len(self._flats))
        gens = poly_ring.gens()
        self._flats_generator = dict()
        self._flats_generator = dict(zip(self._flats, gens))
        MPolynomialIdeal.__init__(self, poly_ring, self._gens_constructor(poly_ring))
        

    def _gens_constructor(self, poly_ring):
        r"""
        Return the generators of augmented Chow ring ideal of
        atom-free presentation.

        EXAMPLES::

            sage: ch = AugmentedChowRingIdeal_atom_free(M=matroids.Wheel(3), R=QQ)
            sage: ch._gens_constructor() #WHERE IS OUTPUT?

        """
        E = list(self._matroid.groundset())
        Q = []
        flats_containing = {x: [] for x in E}
        for F in self._flats:
            for x in F:
                flats_containing[x].append(F)
        for F in self._flats:
            for G in self._flats:
                if not (G > F or F > G):
                        Q.append(self._flats_generator[F]*self._flats_generator[G])
                for x in E:
                    term = poly_ring.zero()
                    for H in flats_containing[x]:
                        term += self._flats_generator[H]
                    Q.append(term**2)

                    if F not in flats_containing[x]:
                        term = poly_ring.zero()
                        for H in flats_containing[x]:
                            term += self._flats_generator[H]
                        Q.append(self._flats_generator[F]*term)

        return Q
            
    def _repr_(self):
        r"""
        EXAMPLE::

            sage: ch = AugmentedChowRingIdeal_atom_free(M=matroids.Wheel(3), R=QQ)
            sage: ch
            Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 
            6 elements with 16 bases of atom-free presentation
        """
        return "Augmented Chow ring ideal of {} of atom-free presentation".format(self._matroid)
    
    def groebner_basis(self):
        """
        Returns the Groebner basis of the augmented Chow ring ideal.
        Return type - ``PolynomialSequence``.

        EXAMPLES::

            sage: from sage.matroids.chow_ring_ideal import AugmentedChowRingIdeal_atom_free
            sage: from sage.matroids.graphic_matroid import GraphicMatroid

            sage: M1 = GraphicMatroid(graphs.CycleGraph(3))
            sage: ch = AugmentedChowRingIdeal_atom_free(M=M1, R=QQ)
            sage: ch.groebner_basis()
            [A0^2, A0*A1, A0*A2, A0*A1, A1^2, A1*A2, A0*A2, A1*A2, A2^2]
        """
        gb = []
        flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        poly_ring = self.ring()
        if Set([]) in flats:
            flats.remove(Set([]))
        for F in flats:
            for G in flats:
                if not (F > G or G > F):
                    gb.append(self._flats_generator[F]*self._flats_generator[G])
                elif F < G:
                    term = poly_ring.zero()
                    for H in flats:
                        if H < F:
                            term += self._flats_generator[H]
                    gb.append(self._flats_generator[F]*(term**self._matroid.rank(Set(G)))*
                            (term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F)))))

        g_basis = PolynomialSequence(poly_ring, [gb])
        return g_basis

            

        

