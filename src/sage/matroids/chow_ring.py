from sage.matroids.chow_ring_ideal import *
from sage.rings.quotient_ring import QuotientRing_nc
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
import sage.misc.latex as latex

#chow rings class needs to be written properly
#gens of chow ring ideal must be debugged
#gb must be a poly sequence
#commented and documented
#use is_groebner()
#matroids - try all of them


class ChowRing(QuotientRing_nc):
    def __init__(self, R, M, augmented):
        self._matroid = M
        self._augmented = augmented
        if augmented:
            self._ideal = AugmentedChowRingIdeal(M, R)
        else:
            self._ideal = ChowRingIdeal(M, R) #check method to get ring
        QuotientRing_nc.__init__(self, R=self._ideal.poly_ring, I=self._ideal, names=self.poly_ring.variable_names(), category=GradedAlgebrasWithBasis(R))

    def _repr_(self):
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        return "%s/%s" % (latex.latex(self.poly_ring), latex.latex(self._ideal))  

