# sage_setup: distribution = sagemath-gap
"""
Symplectic Linear Groups with GAP
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#                     2013 Volker Braun <vbraun.name@gmail.com>
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
from sage.groups.matrix_gps.symplectic import SymplecticMatrixGroup_generic
from sage.misc.cachefunc import cached_method


class SymplecticMatrixGroup_gap(SymplecticMatrixGroup_generic, NamedMatrixGroup_gap, FinitelyGeneratedMatrixGroup_gap):
    r"""
    Symplectic group in GAP.

    EXAMPLES::

        sage: Sp(2,4)                                                                   # needs sage.rings.finite_rings
        Symplectic Group of degree 2 over Finite Field in a of size 2^2

        sage: latex(Sp(4,5))
        \text{Sp}_{4}(\Bold{F}_{5})

    TESTS:

    Check that :issue:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
        sage: G = Sp(4,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """

    @cached_method
    def invariant_form(self):
        """
        Return the quadratic form preserved by the symplectic group.

        OUTPUT: a matrix

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
