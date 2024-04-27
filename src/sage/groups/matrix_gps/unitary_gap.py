# sage.doctest: needs sage.rings.finite_rings
r"""
Unitary Groups `GU(n,q)` and `SU(n,q)` with GAP
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#                     2018 Sebastian Oehms
#                     2018 Travis Scrimshaw
#                     2023 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.groups.matrix_gps.named_group_gap import NamedMatrixGroup_gap
from sage.groups.matrix_gps.unitary import UnitaryMatrixGroup_generic
from sage.misc.cachefunc import cached_method


class UnitaryMatrixGroup_gap(UnitaryMatrixGroup_generic, NamedMatrixGroup_gap, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The general or special unitary group in GAP.

    TESTS:

    Check that :issue:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
        sage: G = GU(3,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """

    @cached_method
    def invariant_form(self):
        """
        Return the hermitian form preserved by the unitary group.

        OUTPUT: a square matrix describing the bilinear form

        EXAMPLES::

            sage: G32 = GU(3,2)
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
