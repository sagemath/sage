from sage.matroids.chow_ring_ideal import *
from sage.rings.quotient_ring import QuotientRing_nc
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
import sage.misc.latex as latex

class ChowRing(QuotientRing_nc):
    def __init__(self, R, M):
        self._matroid = M
        self._ideal = ChowRingIdeal(M, R)
        self.poly_ring = self._ideal.poly_ring
        QuotientRing_nc.__init__(self, R=R, I=self._ideal, names=self.poly_ring.variable_names, category=GradedAlgebrasWithBasis)

    def _repr_short(self):
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        return "%s/%s" % (latex.latex(self.poly_ring), latex.latex(self._ideal))  

class AugmentedChowRing(QuotientRing_nc):
    def __init__(self, R, M):
        self._matroid = M
        self._ideal = AugmentedChowRingIdeal(M, R)
        self.poly_ring = self._ideal.poly_ring
        QuotientRing_nc.__init__(self, R, self._ideal, names=self.poly_ring.variable_names, category=GradedAlgebrasWithBasis)

    def _repr_short(self):
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        return "%s/%s" % (latex.latex(self.poly_ring), latex.latex(self._ideal))