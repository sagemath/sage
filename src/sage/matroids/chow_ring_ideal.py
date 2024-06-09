from sage.all import *
from sage.rings.ideal import Ideal_generic
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.matroids.matroid import Matroid
from sage.matroids.utilities import cmp_elements_key
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class ChowRingIdeal(MPolynomialIdeal):
    def __init__(self, M, R=None, augmented=False):
        if R is None:
            raise ValueError("Ring not given")
        self.augmented = augmented
        self.mat = M
        flats = [X for i in range(1, self.mat.rank())
                 for X in self.mat.flats(i)]
        try:
            self.names = ['A{}'.format(''.join(str(x) for x in sorted(F, key=cmp_elements_key))) for F in flats]
            self.poly_ring = PolynomialRing(R, self.names)
        except ValueError: # variables are not proper names
            self.poly_ring = PolynomialRing(R, 'A', len(flats))
            self.names = self.poly_ring.variable_names()
        self.gens = self.poly_ring.gens()
        Ideal_generic.__init__(self, self.poly_ring, self.gens, coerce=coerce)
        
        

    def groebner_basis(self):
        gb = []
        if self.augmented:
        
        else:



            
            

        

