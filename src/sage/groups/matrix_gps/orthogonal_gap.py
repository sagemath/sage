r"""
Orthogonal Linear Groups with GAP
"""

from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.groups.matrix_gps.named_group_gap import NamedMatrixGroup_gap
from sage.groups.matrix_gps.orthogonal import OrthogonalMatrixGroup_generic
from sage.misc.cachefunc import cached_method


class OrthogonalMatrixGroup_gap(OrthogonalMatrixGroup_generic, NamedMatrixGroup_gap, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The general or special orthogonal group in GAP.

    TESTS:

    Check that :trac:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
        sage: G = GO(3,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """
    @cached_method
    def invariant_bilinear_form(self):
        """
        Return the symmetric bilinear form preserved by the orthogonal
        group.

        OUTPUT:

        A matrix `M` such that, for every group element g, the
        identity `g m g^T = m` holds. In characteristic different from
        two, this uniquely determines the orthogonal group.

        EXAMPLES::

            sage: G = GO(4, GF(7), -1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]

            sage: G = GO(4, GF(7), +1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 6 0]
            [0 0 0 2]

            sage: G = SO(4, GF(7), -1)
            sage: G.invariant_bilinear_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]

        TESTS::

            sage: G.invariant_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]
        """
        m = self.gap().InvariantBilinearForm()['matrix'].matrix()
        m.set_immutable()
        return m

    invariant_form = invariant_bilinear_form # alias (analogues to symplectic and unitary cases)

    @cached_method
    def invariant_quadratic_form(self):
        r"""
        Return the quadratic form preserved by the orthogonal group.

        OUTPUT:

        The matrix `Q` defining "orthogonal" as follows. The matrix
        determines a quadratic form `q` on the natural vector space
        `V`, on which `G` acts, by `q(v) = v Q v^t`. A matrix `M` is
        an element of the orthogonal group if `q(v) = q(v M)` for all
        `v \in V`.

        EXAMPLES::

            sage: G = GO(4, GF(7), -1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = GO(4, GF(7), +1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]

            sage: G = GO(4, QQ)
            sage: G.invariant_quadratic_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: G = SO(4, GF(7), -1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 1 0]
            [0 0 0 1]

        TESTS::

            sage: GO(4, GF(7), -1).invariant_form()
            [0 1 0 0]
            [1 0 0 0]
            [0 0 2 0]
            [0 0 0 2]
        """
        m = self.gap().InvariantQuadraticForm()['matrix'].matrix()
        m.set_immutable()
        return m
