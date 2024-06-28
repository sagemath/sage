from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.matroid import Matroid
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence


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

        names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in self.flats]

        try:
            self.poly_ring = PolynomialRing(R, names) #self.ring
        except ValueError: # variables are not proper names
            self.poly_ring = PolynomialRing(R, 'A', len(self.flats))
           
           
        gens = self.poly_ring.gens()
        self.flats_generator = dict(zip(self.flats, gens))

        
        Q = [gens[i] * gens[i+j+1] for i,F in enumerate(self.flats)
            for j,G in enumerate(self.flats[i+1:]) if not (F < G or G < F)]
        L = [sum(gens[i] for i in flats_containing[x])
            - sum(gens[i] for i in flats_containing[y])
            for j,x in enumerate(E) for y in E[j+1:]]
        

        MPolynomialIdeal.__init__(self, self.poly_ring, Q + L)

    def _repr_(self):
        return "Chow ring ideal of {}".format(self._matroid)

    def groebner_basis(self):  
        gb = list()
        for F in self.flats:
            for G in self.flats: 
                if not (F < G or G < F):
                    gb.append(self.flats_generator[F]*self.flats_generator[G])
                elif Set(F).is_empty():
                    term = self.poly_ring.zero()
                    for H in self.flats:
                        if H < F:
                            term += self.flats_generator[H]
                    gb.append(term**self._matroid.rank(Set(G)))
                elif F < G:
                    term = self.poly_ring.zero()
                    for H in self.flats:
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


        le = len(E)
        Q = list()
        for i,F in enumerate(self.flats):
            for j,G in enumerate(self.flats):
                    if not (F < G or G < F):
                        #print(type(gens[le+i] * gens[le+i+j+1]))
                        Q.append(gens[le+i] * gens[le+j])
                    
        for j,F in enumerate(self.flats):
            for k,x in enumerate(E):
                if F not in flats_containing[x]:
                    #print(type(gens[k]*gens[len(E)+j]))
                    Q.append(gens[k]*gens[le+j])
        #print("this is Q", Q)
                                            
                

        #debug Q and L. coerce must be true. Coerce allows you to add two different elements from two different rings
        #Q =[gens[len(E)+i] * gens[len(E)+i+j+1] for i,F in enumerate(self.flats)
                #for j,G in enumerate(self.flats[i+1:]) if not (F < G or G < F)]
        #Q.append([gens[i]*gens[len(E)+j] for i,x in enumerate(E) for j,F in enumerate(self.flats) if F not in flats_containing[x]])
        L = list()
        for i,x in enumerate(E):
            term = self.poly_ring.zero()
            for j,F in enumerate(self.flats):
                if F not in flats_containing[x]:
                    term += gens[le+j]
            L.append(gens[i] - term)
        
        MPolynomialIdeal.__init__(self, self.poly_ring, Q + L)
    
    def _repr_(self): #use single underscore
        return "Augmented Chow ring ideal of {}".format(self._matroid)

    def matroid(self):
        return self._matroid
    
    def flat_generator(self):
        return self.flats_generator
    
    


    
    def groebner_basis(self, atom_free=False):
        gb = []
        flats = self.flats
        if Set([]) in flats:
            flats.remove(Set([]))
        if atom_free:
            for F in flats:
                for G in flats:
                    if not (F > G or G > F):
                        gb.append(self.flats_generator[F]*self.flats_generator[G])
                    elif F < G:
                        term = self.poly_ring.zero()
                        for H in self.flats:
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



        




            
            

        

