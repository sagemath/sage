"""
Symplectic Linear Groups with GAP
"""

from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.groups.matrix_gps.named_group_gap import NamedMatrixGroup_gap
from sage.groups.matrix_gps.symplectic import SymplecticMatrixGroup_generic
from sage.misc.cachefunc import cached_method


class SymplecticMatrixGroup_gap(SymplecticMatrixGroup_generic, NamedMatrixGroup_gap, FinitelyGeneratedMatrixGroup_gap):
    r"""
    Symplectic group in GAP.

    EXAMPLES::

        sage: Sp(2,4)
        Symplectic Group of degree 2 over Finite Field in a of size 2^2

        sage: latex(Sp(4,5))
        \text{Sp}_{4}(\Bold{F}_{5})

    TESTS:

    Check that :trac:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
        sage: G = Sp(4,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """

    @cached_method
    def invariant_form(self):
        """
        Return the quadratic form preserved by the symplectic group.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: Sp(4, GF(3)).invariant_form()
            [0 0 0 1]
            [0 0 1 0]
            [0 2 0 0]
            [2 0 0 0]
        """
        m = self.gap().InvariantBilinearForm()['matrix'].matrix()
        m.set_immutable()
        return m
