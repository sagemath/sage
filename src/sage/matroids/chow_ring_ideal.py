#from sage.all import *
#from sage.rings.ideal import Ideal_generic
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.matroid import Matroid
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.morphism import _tensor_product_ring

class ChowRingIdeal(MPolynomialIdeal):
    def __init__(self, M, R=None, augmented=False):
        if R is None:
            raise ValueError("Ring not given")
        self.augmented = augmented
        self._matroid = M
        self.flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(self.flats):
            for x in F:
                flats_containing[x].append(i)
        self.names = dict()

        if self.augmented:
            try:
                P1 = PolynomialRing(R, 'A', len(E))
                names_flats = ['B{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self.flats]
                P2 = PolynomialRing(R, names_flats)
                self.poly_ring = _tensor_product_ring(P1, P2)

            except:
                P1 = PolynomialRing(R, 'A', len(E))
                P2 = PolynomialRing(R, 'B', len(self.flats))
                self.poly_ring = _tensor_product_ring(P1, P2)

            gens = self.poly_ring.gens()
            for i,x in enumerate(E):
                self.names[x] = gens[i]
            for i,F in enumerate(self.flats):
                self.names[F] = gens[len(E) + i]

            Q =[gens[len(E)+i] * gens[len(E)+i+j+1] for i,F in enumerate(self.flats)
                for j,G in enumerate(self.flats[i+1:]) if not (F < G or G < F)]
            Q.append([gens[i]*gens[len(E)+j] for i,x in enumerate(E) for j,F in enumerate(self.flats) if F not in flats_containing[x]])
            L = list()
            for i,x in enumerate(E):
                for j,F in enumerate(self.flats):
                    term = 0
                    if F not in flats_containing[x]:
                        term += gens[len(E)+j]
                L.append(gens[i] - term)
            print(L)
            self.gens = Q + L


        else: 
            try:
                names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self.flats]
                self.poly_ring = PolynomialRing(R, names)
                for F in self.flats:
                    for i in range(len(self.poly_ring.gens())):
                        self.names[F] = self.poly_ring.gens()[i]
            except ValueError: # variables are not proper names
                self.poly_ring = PolynomialRing(R, 'A', len(self.flats))
                for i in range(len(self.flats)):
                    self.names[self.flats[i]] = self.poly_ring.gens()[i]

            gens = self.poly_ring.gens()
            Q = [gens[i] * gens[i+j+1] for i,F in enumerate(self.flats)
                for j,G in enumerate(self.flats[i+1:]) if not (F < G or G < F)]
            L = [sum(gens[i] for i in flats_containing[x])
                - sum(gens[i] for i in flats_containing[y])
                for j,x in enumerate(E) for y in E[j+1:]]
            self.gens = Q + L
        

        MPolynomialIdeal.__init__(self, self.poly_ring, self.gens)

    def __repr__(self):
        if self.augmented:
            return "Augmented Chow ring ideal of {}".format(self._matroid)
        else:
            return "Chow ring ideal of {}".format(self._matroid)

    def groebner_basis(self):
        if self.augmented:
            print("Augmented gb")
        else:    
            gb = list()
            for F in self.flats:
                for G in self.flats:
                    if not (F < G or G < F):
                        gb.append(self.names[F]*self.names[G])
                    elif Set(F).is_empty():
                        term = 0
                        for H in self.flats:
                            if H < F:
                                term += self.names[H]
                        gb.append(term**self._matroid.rank(Set(G)))
                    elif F < G:
                        term = 0
                        for H in self.flats:
                            if H < F:
                                term += self.names[H]
                        gb.append(term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F))))

        return gb        



    
        
        




            
            

        

