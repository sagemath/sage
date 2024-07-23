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
        self._matroid = M
        self._augmented = augmented
        if augmented:
            if presentation=='fy':
                self._ideal = AugmentedChowRingIdeal_fy(M, R)
            elif presentation=='atom-free':
                self._ideal = AugmentedChowRingIdeal_fy(M, R)
        else:
            self._ideal = ChowRingIdeal_nonaug(M, R) #check method to get ring
        QuotientRing_nc.__init__(self, R=self._ideal.poly_ring, I=self._ideal, names=self._ideal.poly_ring.variable_names(), category=GradedAlgebrasWithBasis(R))

    def _repr_(self):
        return "Chow ring of {}".format(self._matroid)    

    def _latex_(self):
        import sage.misc.latex as latex
        return "%s/%s" % (latex.latex(self.poly_ring), latex.latex(self._ideal))
