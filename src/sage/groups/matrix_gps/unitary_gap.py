r"""
Unitary Groups `GU(n,q)` and `SU(n,q)` with GAP
"""

from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.groups.matrix_gps.named_group_gap import NamedMatrixGroup_gap
from sage.groups.matrix_gps.unitary import UnitaryMatrixGroup_generic
from sage.misc.cachefunc import cached_method


class UnitaryMatrixGroup_gap(UnitaryMatrixGroup_generic, NamedMatrixGroup_gap, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The general or special unitary group in GAP.

    TESTS:

    Check that :trac:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
        sage: G = GU(3,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """

    @cached_method
    def invariant_form(self):
        """
        Return the hermitian form preserved by the unitary group.

        OUTPUT:

        A square matrix describing the bilinear form

        EXAMPLES::

            sage: G32=GU(3,2)
            sage: G32.invariant_form()
            [0 0 1]
            [0 1 0]
            [1 0 0]
        """
        d = self.degree()
        R = self.base_ring()
        # note that self.gap().InvariantSesquilinearForm()['matrix'].matrix().base_ring() != R for example for self = GU(3.2)
        # therefore we have to coerce into the right matrix space
        from sage.matrix.constructor import matrix
        m = matrix(R, d, d, self.gap().InvariantSesquilinearForm()['matrix'].matrix())
        m.set_immutable()
        return m
