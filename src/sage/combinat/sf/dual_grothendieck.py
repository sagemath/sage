"""
Dual Grothendieck Symmetric Functions
"""

from sage.categories.hopf_algebras import HopfAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partition, Partitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.skew_tableau import SemistandardSkewTableaux
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class DualGrothendiecks(CombinatorialFreeModule):
    """
    Dual (symmetric) Grothendieck functions.

    EXAMPLES::

        sage: G = DualGrothendiecks(QQ)
        sage: G._s.print_options(sorting_reverse=True)
        sage: G[2,1].schur_expand()
        s[2, 1] + beta*s[2]
        sage: G[1,1,1].schur_expand()
        s[1, 1, 1] + 2*beta*s[1, 1] + beta^2*s[1]
        sage: G[3,3,3].schur_expand()
        s[3, 3, 3] + 2*beta*s[3, 3, 2] + 3*beta^2*s[3, 3, 1] + 4*beta^3*s[3, 3]
         + beta^2*s[3, 2, 2] + 2*beta^3*s[3, 2, 1] + 3*beta^4*s[3, 2]
         + beta^4*s[3, 1, 1] + 2*beta^5*s[3, 1] + beta^6*s[3]
        sage: G[2,1] * G[2,1]
        dG[4, 2] + dG[4, 1, 1] - beta*dG[4, 1] + dG[3, 3] + 2*dG[3, 2, 1]
         - 3*beta*dG[3, 2] + dG[3, 1, 1, 1] - 3*beta*dG[3, 1, 1]
         + 2*beta^2*dG[3, 1] + dG[2, 2, 2] + dG[2, 2, 1, 1] - 3*beta*dG[2, 2, 1]
         + beta^2*dG[2, 2] - beta*dG[2, 1, 1, 1] + 2*beta^2*dG[2, 1, 1] - beta^3*dG[2, 1]
    """
    def __init__(self, base_ring):
        poly = PolynomialRing(base_ring, 'beta')
        beta = poly.gen()
        cat = HopfAlgebras(poly).Filtered().WithBasis()
        self._s = SymmetricFunctions(poly).s()
        super().__init__(poly, Partitions(), prefix="dG", bracket=False, category=cat)

        def expand(la):
            if not la:
                return self._s.zero()
            m = sum(la)
            ell = len(la)
            ret = {}
            for n in range(m, la[0]-1, -1):
                for mu in Partitions(n, outer=la, inner=[la[0]]):
                    coeff = beta ** (m - n)
                    coeff *= sum(1 for T in SemistandardSkewTableaux([la, mu], max_entry=ell-1)
                                 if all(all(val is None for val in row) or row[-1] <= ind
                                        for ind, row in enumerate(T)))
                    ret[mu] = coeff
            return self._s._from_dict(ret)

        self._phi = self.module_morphism(on_basis=expand, codomain=self._s, unitriangular="upper")
        self._phi_inv = ~self._phi

    def product(self, left, right):
        return self._phi_inv(self._phi(left) * self._phi(right))

    def __getitem__(self, la):
        """
        EXAMPLES::

            sage: G[2]
            dG[2]
            sage: G[[]]
            dG[]
            sage: G[2,1]
            dG[2, 1]
        """
        if la in ZZ:
            la = Partition([la])
        else:
            la = Partition(la)
        return self.monomial(la)

    class Element(CombinatorialFreeModule.Element):
        def schur_expand(self):
            return self.parent()._phi(self)
