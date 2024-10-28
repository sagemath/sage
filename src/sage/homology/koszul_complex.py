"""
Koszul Complexes
"""
########################################################################
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.combination import rank
from sage.arith.misc import binomial
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex_class

import itertools


class KoszulComplex(ChainComplex_class, UniqueRepresentation):
    r"""
    A Koszul complex.

    Let `R` be a ring and consider `x_1, x_2, \ldots, x_n \in R`. The
    *Koszul complex* `K_*(x_1, \ldots, x_n)` is given by defining a
    chain complex structure on the exterior algebra `\bigwedge^n R` with
    the basis `e_{i_1} \wedge \cdots \wedge e_{i_a}`. The differential is
    given by

    .. MATH::

        \partial(e_{i_1} \wedge \cdots \wedge e_{i_a}) =
        \sum_{r=1}^a (-1)^{r-1} x_{i_r} e_{i_1} \wedge \cdots \wedge
        \hat{e}_{i_r} \wedge \cdots \wedge e_{i_a},

    where `\hat{e}_{i_r}` denotes the omitted factor.

    Alternatively we can describe the Koszul complex by considering the
    basic complex `K_{x_i}`

    .. MATH::

        0 \rightarrow R \xrightarrow{x_i} R \rightarrow 0.

    Then the Koszul complex is given by
    `K_*(x_1, \ldots, x_n) = \bigotimes_i K_{x_i}`.

    INPUT:

    - ``R`` -- the base ring
    - ``elements`` -- tuple of elements of `R`

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: K = KoszulComplex(R, [x,y])
        sage: ascii_art(K)
                                [-y]
                    [x y]       [ x]
         0 <-- C_0 <------ C_1 <----- C_2 <-- 0
        sage: K = KoszulComplex(R, [x,y,z])
        sage: ascii_art(K)
                                  [-y -z  0]       [ z]
                                  [ x  0 -z]       [-y]
                    [x y z]       [ 0  x  y]       [ x]
         0 <-- C_0 <-------- C_1 <----------- C_2 <----- C_3 <-- 0
        sage: K = KoszulComplex(R, [x+y*z,x+y-z])
        sage: ascii_art(K)
                                                [-x - y + z]
                    [  y*z + x x + y - z]       [   y*z + x]
         0 <-- C_0 <---------------------- C_1 <------------- C_2 <-- 0

    REFERENCES:

    - :wikipedia:`Koszul_complex`
    """
    @staticmethod
    def __classcall_private__(cls, R=None, elements=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: R.<x,y,z> = QQ[]
            sage: K1 = KoszulComplex(R, [x,y,z])
            sage: K2 = KoszulComplex(R, (x,y,z))
            sage: K3 = KoszulComplex((x,y,z))
            sage: K1 is K2 and K2 is K3
            True

        Check some corner cases::

            sage: K1 = KoszulComplex(ZZ)
            sage: K2 = KoszulComplex(())
            sage: K3 = KoszulComplex(ZZ, [])
            sage: K1 is K2 and K2 is K3
            True
            sage: K1 is KoszulComplex()
            True
        """
        if elements is None:
            if R is None:
                R = ()
            elements = R
            if not elements:
                R = ZZ  # default to ZZ as the base ring if no elements are given
            elif isinstance(R, Parent):
                elements = ()
            else:
                R = elements[0].parent()
        elif R is None:  # elements is not None
            R = elements[0].parent()
        return super().__classcall__(cls, R, tuple(elements))

    def __init__(self, R, elements):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: K = KoszulComplex(R, [x,y])
            sage: TestSuite(K).run()
        """
        # Generate the differentials
        self._elements = elements
        n = len(elements)
        I = list(range(n))
        diff = {}
        zero = R.zero()
        for i in I:
            M = matrix(R, binomial(n, i), binomial(n, i+1), zero)
            j = 0
            for comb in itertools.combinations(I, i+1):
                for k, val in enumerate(comb):
                    r = rank(comb[:k] + comb[k+1:], n, False)
                    M[r, j] = (-1)**k * elements[val]
                j += 1
            M.set_immutable()
            diff[i+1] = M
        diff[0] = matrix(R, 0, 1, zero)
        diff[0].set_immutable()
        diff[n+1] = matrix(R, 1, 0, zero)
        diff[n+1].set_immutable()
        ChainComplex_class.__init__(self, ZZ, ZZ(-1), R, diff)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: KoszulComplex(R, [x,y,z])
            Koszul complex defined by (x, y, z) over
             Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: KoszulComplex(ZZ, [])
            Trivial Koszul complex over Integer Ring
        """
        if not self._elements:
            return "Trivial Koszul complex over {}".format(self.base_ring())
        return "Koszul complex defined by {} over {}".format(self._elements, self.base_ring())
