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
        self.flats = [X for i in range(1, self._matroid.rank()) #_flats. NOT NEEDED AS AN ATTRIBUTE. USE NAMES
                 for X in self._matroid.flats(i)]
        
        E = list(self._matroid.groundset())
        flats_containing = {x: [] for x in E}
        for i,F in enumerate(self.flats):
            for x in F:
                flats_containing[x].append(i)
        self.flat_generator = dict()
        try:
            names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self.flats]
            self.poly_ring = PolynomialRing(R, names) #self.ring
            for F in self.flats:
                for i in range(len(self.poly_ring.gens())):
                    self.flat_generator[F] = self.poly_ring.gens()[i] #change self.names to self.flat_generator
        except ValueError: # variables are not proper names
            self.poly_ring = PolynomialRing(R, 'A', len(self.flats))
            for i in range(len(self.flats)):
                self.flat_generator[self.flats[i]] = self.poly_ring.gens()[i]

        gens = self.poly_ring.gens()
        Q = [gens[i] * gens[i+j+1] for i,F in enumerate(self.flats)
            for j,G in enumerate(self.flats[i+1:]) if not (F < G or G < F)]
        L = [sum(gens[i] for i in flats_containing[x])
            - sum(gens[i] for i in flats_containing[y])
            for j,x in enumerate(E) for y in E[j+1:]]
        self.gens = Q + L
        

        MPolynomialIdeal.__init__(self, self.poly_ring, self.gens)

    def _repr_(self):
        return "Chow ring ideal of {}".format(self._matroid)

    def groebner_basis(self):  
        gb = list()
        for F in self.flats:
            for G in self.flats: #write from F not G
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


    def matroid(self):
        return Matroid(self._matroid)
    
    def flat_generator(self):
        return dict(self.flat_generator)

#get matroid method, called matroid, returning the matroid
#get names


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
        self.flats_generator = dict()
        #names_groundset = ['A{}'.format(''.join(str(x))) for x in E]
        #names_flats = ['B{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self.flats]
        self.poly_ring = PolynomialRing(R, 'A', len(E) + len(self.flats))
        gens = self.poly_ring.gens()
        for i,x in enumerate(E):
                self.flats_generator[x] = gens[i]
        for i,F in enumerate(self.flats):
            self.flats_generator[F] = gens[len(E) + i]

        print(gens)
    
        Q = list()
        for i,F in enumerate(self.flats):
            for j,G in enumerate(self.flats):
                    if not (F < G or G < F):
                        print(type(gens[len(E)+i] * gens[len(E)+i+j+1]))
                        Q.append(gens[len(E)+i] * gens[len(E)+i+j+1])
                    
        for j,F in enumerate(self.flats):
            for k,x in enumerate(E):
                if F not in flats_containing[x]:
                    print(type(gens[k]*gens[len(E)+j]))
                    Q.append(gens[k]*gens[len(E)+j])
        print("this is Q", Q)
                                            
                

        #debug Q and L. coerce must be true. Coerce allows you to add two different elements from two different rings
        #Q =[gens[len(E)+i] * gens[len(E)+i+j+1] for i,F in enumerate(self.flats)
                #for j,G in enumerate(self.flats[i+1:]) if not (F < G or G < F)]
        #Q.append([gens[i]*gens[len(E)+j] for i,x in enumerate(E) for j,F in enumerate(self.flats) if F not in flats_containing[x]])
        L = list()
        for i,x in enumerate(E):
            term = 0
            for j,F in enumerate(self.flats):
                if F not in flats_containing[x]:
                    term += gens[len(E)+j]
        L.append(gens[i] - term)
        self.gens = Q + L
        MPolynomialIdeal.__init__(self, self.poly_ring, self.gens, coerce=False)
    
    def _repr_short(self): #use single underscore
        return "Augmented Chow ring ideal of {}".format(self._matroid)

    
    def groebner_basis(self, atom_free=False):
        #list returned or iterator returned? - neither - polynomial_sequence_generic object
        gb = []
        flats = self.flats
        if Set([]) in flats:
            flats.remove(Set([]))
        if atom_free:
            for F in flats:
                for G in flats:
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
                for F in flats:
                    for G in flats:
                        term = 0
                        for H in flats:
                            if i in Set(H):
                                term += self.names[H]
                        gb.append(self.names[i] + term)

                        if i in Set(F):
                            term = 0
                            for H in flats:
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
                            for H in flats:
                                if H < F:
                                    term += self.names[H]
                            gb.append(term**(self._matroid.rank(Set(G))-self._matroid.rank(Set(F))))

        return gb



        




            
            

        

