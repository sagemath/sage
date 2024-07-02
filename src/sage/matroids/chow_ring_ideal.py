from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.matroid import Matroid
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

r"""
Chow ring ideals of matroids

AUTHORS:

- Travis Scrimshaw 
- Shriya M

These are the classes of Chow ring ideals for matroids. There are two classes 
created - Chow ring ideal and augmented Chow ring ideal. The augmented
Chow ring ideal has two different presentations implemented - the Feitchner-
Yuzvinsky presentation and atom-free presentation. Both classes have 
``grobner_basis()`` methods as well, as an explicit Groebner basis is known 
in each case.

REFERENCES

- :arxiv:`2309.14312`
- :arxiv:`2111.00393`
"""
#*****************************************************************************
#       Copyright (C) 2024 Travis Scrimshaw
#                     2024 Shriya M
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

class ChowRingIdeal(MPolynomialIdeal):
    r"""
    The class of Chow ring ideal, a multi-polynomial ideal. 
    Base class - ``MPolynomialIdeal``.

    INPUT:


    - `M` -- a matroid.
    - `R` -- a ring.

    OUTPUT: Chow ring ideal of matroid `M`.

    EXAMPLES::

    Chow ring ideal of uniform matroid of rank 3 on 6 elements::

        sage: ch = ChowRingIdeal(M=matroids.Uniform(3,6), R=QQ)
        sage: ch
        Chow ring ideal of U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
        {3: {{0, 1, 2, 3, 4, 5}}}
        sage: ch = ChowRingIdeal(M=matroids.catalog.Fano(), R=QQ)
        Chow ring ideal of Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
    """

    def __init__(self, M, R):
        self._matroid = M
        flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(flats):
            for x in F:
                flats_containing[x].append(i)

        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]

        try:
            self.poly_ring = PolynomialRing(R, names) #self.ring
        except ValueError: # variables are not proper names
            self.poly_ring = PolynomialRing(R, 'A', len(self.flats))
           
           
        gens = self.poly_ring.gens()
        self.flats_generator = dict(zip(flats, gens))

        
        Q = [gens[i] * gens[i+j+1] for i,F in enumerate(flats)
            for j,G in enumerate(flats[i+1:]) if not (F < G or G < F)]
        L = [sum(gens[i] for i in flats_containing[x])
            - sum(gens[i] for i in flats_containing[y])
            for j,x in enumerate(E) for y in E[j+1:]]
        

        MPolynomialIdeal.__init__(self, self.poly_ring, Q + L)

    def __repr__(self):
        return "Chow ring ideal of {}".format(self._matroid)

    def groebner_basis(self):
        r"""
        Returns the Groebner basis of the Chow ring ideal of consideration.
        Return type - ``PolynomialSequence``.

        EXAMPLES::

            sage: ch = ChowRingIdeal(M=BasisMatroid(groundset='abc', bases=['ab', 'ac']), R=QQ)
            sage: ch.groebner_basis()
            [Aa^2, Aa*Abc, Aa*Abc, Abc^2]
            sage: ch.groebner_basis().is_groebner()
            True
        
        Another example would be the Groebner basis of the Chow ring ideal of
        the Non-Fano matroid::

            sage: ch = ChowRingIdeal(M=matroids.catalog.NonFano(), R=QQ)
            sage: ch.groebner.basis()
            Polynomial Sequence with 592 Polynomials in 16 Variables
        """

        flats = list(self.flats_generator.keys())
        gb = list()
        for F in flats:
            for G in flats: 
                if not (F < G or G < F):
                    gb.append(self.flats_generator[F]*self.flats_generator[G])
                elif Set(F).is_empty():
                    term = self.poly_ring.zero()
                    for H in flats:
                        if H < F:
                            term += self.flats_generator[H]
                    gb.append(term**self._matroid.rank(Set(G)))
                elif F < G:
                    term = self.poly_ring.zero()
                    for H in flats:
                        if H < F:
                            term += self.flats_generator[H]
                        gb.append(term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F))))

        g_basis = PolynomialSequence(self.poly_ring, [gb])
        return g_basis


    def matroid(self):
        return self._matroid
    
    def flats_generator(self):
        return dict(self.flats_generator)



class AugmentedChowRingIdeal(MPolynomialIdeal):
    r"""
    The class of Chow ring ideal, a multi-polynomial ideal. 
    Base class - ``MPolynomialIdeal``.

    INPUT:


    - `M` -- a matroid.
    - `R` -- a ring.
    - ``atom_free`` -- a Boolean value, default value set to ``False``. When
    ``True``, it returns the atom-free presentation of the augmented Chow
    ring. If ``False``, it returns the Feitchner-Yuzvinsky presentation of the
    augmented Chow ring.

    OUTPUT: augmented Chow ring ideal of matroid `M`.

    EXAMPLES::

    Augmented Chow ring ideal of Wheel matroid of rank 3::

        sage: ch = AugumentedChowRingIdeal(M=matroids.Wheel(3), R=QQ)
        sage: ch
        Augmented Chow ring ideal of Wheel(3): Regular matroid of rank 3 on 
        6 elements with 16 bases of Feitchner-Yuzvinsky presentation
    """
    def __init__(self, M, R, atom_free=False):
        self.atom_free = atom_free
        self._matroid = M
        flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for F in flats:
            for x in F:
                flats_containing[x].append(F)
        self.flats_generator = dict()
        try:
            names_groundset = ['A{}'.format(''.join(str(x))) for x in E]
            names_flats = ['B{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]
            self.poly_ring = PolynomialRing(R, names_groundset + names_flats)
        except ValueError:
            self.poly_ring = PolynomialRing(R, 'A', len(E) + len(flats))
        gens = self.poly_ring.gens()
        for i,x in enumerate(E):
            self.flats_generator[x] = gens[i]
        for i,F in enumerate(flats):
            self.flats_generator[F] = gens[len(E) + i]


        le = len(E)
        Q = list()
        for i,F in enumerate(flats):
            for j,G in enumerate(flats):
                    if not (F < G or G < F):
                        Q.append(gens[le+i] * gens[le+j])
        
        for j,F in enumerate(flats):
            for k,x in enumerate(E):
                if F not in flats_containing[x]:
                    if self.atom_free:
                        term = self.poly_ring.zero()
                        for G in flats_containing[x]:
                            term += self.flats_generator[G]
                        Q.append(self.flats_generator[F]*term)
                    else:
                        Q.append(gens[k]*gens[le+j])
                   

        if self.atom_free:
            for i in E:
                term = self.poly_ring.zero()
                for F in flats_containing[i]:
                    term += self.flats_generator[F]
                Q.append(term**2)
            
            MPolynomialIdeal.__init__(self, self.poly_ring, Q)
        else:
            L = list()
            for i,x in enumerate(E):
                term = self.poly_ring.zero()
                for j,F in enumerate(flats):
                    if F not in flats_containing[x]:
                        term += gens[le+j]
                L.append(gens[i] - term)
            
            MPolynomialIdeal.__init__(self, self.poly_ring, Q + L)
    
    def __repr__(self):
        if self.atom_free:
            return "Augmented Chow ring ideal of {} of atom-free presentation".format(self._matroid)
        else:
            return "Augmented Chow ring ideal of {} of Feitchner-Yuzvinsky presentation".format(self._matroid)

    def matroid(self):
        return self._matroid
    
    def flat_generator(self):
        return dict(self.flats_generator)
    
    def groebner_basis(self):
        r"""
        Returns the Groebner basis of the augmented Chow ring ideal.
        Return type - ``PolynomialSequence``.

        EXAMPLES::

            sage: ch = AugmentedChowRingIdeal(M=matroids.catalog.Fano(), R=QQ)
            sage: ch.groebner_basis()
            Polynomial Sequence with 2744 Polynomials in 21 Variables
        
        Another example would be the Groebner basis of the augmented Chow ring
        ideal (atom-free presentation) of the Graphic matroid of cycle graph of
        3 vertices::

            sage: M1 = GraphicMatroid(graphs.CycleGraph(3))
            sage: ch = AugmentedChowRingIdeal(M=M1, R=QQ, atom_free=True)
            sage: ch.groebner_basis()
            [B0^2, B0*B1, B0*B2, B0*B1, B1^2, B1*B2, B0*B2, B1*B2, B2^2]
        """
        gb = []
        flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        if Set([]) in flats:
            flats.remove(Set([]))
        if self.atom_free:
            for F in flats:
                for G in flats:
                    if not (F > G or G > F):
                        gb.append(self.flats_generator[F]*self.flats_generator[G])
                    elif F < G:
                        term = self.poly_ring.zero()
                        for H in flats:
                            if H < F:
                                term += self.flats_generator[H]
                        gb.append(self.flats_generator[F]*(term**self._matroid.rank(Set(G)))*
                                (term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F)))))
        
        else:
            E = list(self._matroid.groundset())
            for i in E:
                for F in flats:
                    for G in flats:
                        term = self.poly_ring.zero()
                        for H in flats:
                            if i in Set(H):
                                term += self.flats_generator[H]
                        gb.append(self.flats_generator[i] + term)

                        if i in Set(F):
                            term = self.poly_ring.zero()
                            for H in flats:
                                if H < F:
                                    term += self.flats_generator[H]
                            gb.append(self.flats_generator[i]*(term**self._matroid.rank(Set(G)))*
                                (term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F)))))
                        
                        elif not i in Set(F):
                            gb.append(self.flats_generator[i]*self.flats_generator[F])
                        
                        elif not (F < G or G < F):
                            gb.append(self.flats_generator[F]*self.flats_generator[G])
                        
                        elif F < G:
                            term = self.poly_ring.zero()
                            for H in flats:
                                if H < F:
                                    term += self.flats_generator[H]
                            gb.append(term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F))))

        g_basis = PolynomialSequence(self.poly_ring, [gb])
        return g_basis



        




            
            

        

