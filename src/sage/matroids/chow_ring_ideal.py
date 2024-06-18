#from sage.all import *
#from sage.rings.ideal import Ideal_generic
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.matroid import Matroid
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.morphism import _tensor_product_ring

class ChowRingIdeal(MPolynomialIdeal):
    def __init__(self, M, R):
        self._matroid = M
        self.flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(self.flats):
            for x in F:
                flats_containing[x].append(i)
        self.names = dict()
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
        return "Chow ring ideal of {}".format(self._matroid)

    def groebner_basis(self):  
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


class AugmentedChowRingIdeal(MPolynomialIdeal):
    def __init__(self, M, R):
        self._matroid = M
        self.flats = [X for i in range(1, self._matroid.rank())
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(self.flats):
            for x in F:
                flats_containing[x].append(i)
        self.names = dict()
        #names_groundset = ['A{}'.format(''.join(str(x))) for x in E]
        #names_flats = ['B{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self.flats]
        self.poly_ring = PolynomialRing(R, 'A', len(E) + len(self.flats))
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
            term = 0
            for j,F in enumerate(self.flats):
                if F not in flats_containing[x]:
                    term += gens[len(E)+j]
        L.append(gens[i] - term)
        self.gens = Q + L
        MPolynomialIdeal.__init__(self, self.poly_ring, self.gens, coerce=False)
    
    def __repr__(self):
        return "Augmented Chow ring ideal of {}".format(self._matroid)

    
    def groebner_basis(self, atom_free=False):
        #add removal of empty flat line
        #list returned or iterator returned?
        gb = []
        if atom_free:
            for F in self.flats:
                for G in self.flats:
                    if not (F > G or G > F):
                        gb.append(self.names[F]*self.names[G])
                    elif F < G:
                        term = 0
                        for H in self.flats:
                            if H < F:
                                term += self.names[H]
                        gb.append(self.names[F]*(term**self._matroid.rank(Set(G)))*
                                (term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F)))))
        
        else:
            E = list(self._matroid.groundset())
            for i in E:
                for F in self.flats:
                    for G in self.flats:
                        term = 0
                        for H in self.flats:
                            if i in Set(H):
                                term += self.names[H]
                        gb.append(self.names[i] + term)

                        if i in Set(F):
                            term = 0
                            for H in self.flats:
                                if H < F:
                                    term += self.names[H]
                            gb.append(self.names[i]*(term**self._matroid.rank(Set(G)))*
                                (term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F)))))
                        
                        elif not i in Set(F):
                            gb.append(self.names[i]*self.names[F])
                        
                        elif not (F < G or G < F):
                            gb.append(self.names[F]*self.names[G])
                        
                        elif F < G:
                            term = 0
                            for H in self.flats:
                                if H < F:
                                    term += self.names[H]
                            gb.append(term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F))))

        return gb



        




            
            

        

