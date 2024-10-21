# sage_setup: distribution = sagemath-gap
# sage.doctest: needs sage.libs.gap sage.modules
"""
Heisenberg Group

AUTHORS:

- Hilder Vitor Lima Pereira (2017-08): initial version
"""

#*****************************************************************************
#    Copyright (C) 2017 Hilder Vitor Lima Pereira <hilder.vitor at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.latex import latex
from sage.matrix.matrix_space import MatrixSpace
from sage.categories.groups import Groups
from sage.categories.rings import Rings
from sage.rings.integer_ring import ZZ
from copy import copy

class HeisenbergGroup(UniqueRepresentation, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The Heisenberg group of degree `n`.

    Let `R` be a ring, and let `n` be a positive integer.
    The Heisenberg group of degree `n` over `R` is a multiplicative
    group whose elements are matrices with the following form:

    .. MATH::

        \begin{pmatrix}
        1 & x^T & z \\
        0 & I_n & y \\
        0 &  0  & 1
        \end{pmatrix},

    where `x` and `y` are column vectors in `R^n`, `z` is a scalar in `R`,
    and `I_n` is the identity matrix of size `n`.

    INPUT:

    - ``n`` -- the degree of the Heisenberg group

    - ``R`` -- (default: `\ZZ`) the ring `R` or a positive integer as
      a shorthand for the ring `\ZZ/R\ZZ`

    EXAMPLES::

        sage: H = groups.matrix.Heisenberg(); H
        Heisenberg group of degree 1 over Integer Ring
        sage: H.gens()
        (
        [1 1 0]  [1 0 0]  [1 0 1]
        [0 1 0]  [0 1 1]  [0 1 0]
        [0 0 1], [0 0 1], [0 0 1]
        )
        sage: X, Y, Z = H.gens()
        sage: Z * X * Y**-1
        [ 1  1  0]
        [ 0  1 -1]
        [ 0  0  1]
        sage: X * Y * X**-1 * Y**-1 == Z
        True

        sage: H = groups.matrix.Heisenberg(R=5); H
        Heisenberg group of degree 1 over Ring of integers modulo 5
        sage: H = groups.matrix.Heisenberg(n=3, R=13); H
        Heisenberg group of degree 3 over Ring of integers modulo 13

    REFERENCES:

    - :wikipedia:`Heisenberg_group`
    """
    @staticmethod
    def __classcall_private__(cls, n=1, R=0):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: H1 = groups.matrix.Heisenberg(n=2, R=5)
            sage: H2 = groups.matrix.Heisenberg(n=2, R=ZZ.quo(5))
            sage: H1 is H2
            True

            sage: H1 = groups.matrix.Heisenberg(n=2)
            sage: H2 = groups.matrix.Heisenberg(n=2, R=ZZ)
            sage: H1 is H2
            True
        """
        if n not in ZZ or n <= 0:
            raise TypeError("degree of Heisenberg group must be a positive integer")
        if R in ZZ:
            if R == 0:
                R = ZZ
            elif R > 1:
                R = ZZ.quo(R)
            else:
                raise ValueError("R must be a positive integer")
        elif R is not ZZ and R not in Rings().Finite():
            raise NotImplementedError("R must be a finite ring or ZZ")
        return super().__classcall__(cls, n, R)

    def __init__(self, n=1, R=0):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: H = groups.matrix.Heisenberg(n=2, R=5)
            sage: TestSuite(H).run()  # long time
            sage: H.category()
            Category of finitely generated finite enumerated groups
            sage: H = groups.matrix.Heisenberg(n=2, R=4)
            sage: TestSuite(H).run()  # long time
            sage: H = groups.matrix.Heisenberg(n=3)
            sage: TestSuite(H).run(max_runs=30, skip='_test_elements')  # long time
            sage: H = groups.matrix.Heisenberg(n=2, R=GF(4))
            sage: TestSuite(H).run()  # long time

        TESTS::

            sage: groups.matrix.Heisenberg(n=2, R=ZZ).category()
            Category of finitely generated infinite enumerated groups
        """
        def elementary_matrix(i, j, val, MS):
            elm = copy(MS.one())
            elm[i,j] = val
            elm.set_immutable()
            return elm

        self._n = n
        self._ring = R
        # We need the generators of the ring as a commutative additive group
        if self._ring is ZZ:
            ring_gens = [self._ring.one()]
        else:
            if self._ring.cardinality() == self._ring.characteristic():
                ring_gens = [self._ring.one()]
            else:
                # This is overkill, but is the only way to ensure
                #   we get all of the elements
                ring_gens = list(self._ring)

        dim = ZZ(n + 2)
        MS = MatrixSpace(self._ring, dim)
        gens_x = [elementary_matrix(0, j, gen, MS)
                  for j in range(1, dim-1) for gen in ring_gens]
        gens_y = [elementary_matrix(i, dim-1, gen, MS)
                  for i in range(1, dim-1) for gen in ring_gens]
        gen_z = [elementary_matrix(0, dim-1, gen, MS) for gen in ring_gens]
        gens = gens_x + gens_y + gen_z

        from sage.libs.gap.libgap import libgap
        gap_gens = [libgap(single_gen) for single_gen in gens]
        gap_group = libgap.Group(gap_gens)

        cat = Groups().FinitelyGenerated()
        if self._ring in Rings().Finite():
            cat = cat.Finite()
        else:
            cat = cat.Infinite()

        FinitelyGeneratedMatrixGroup_gap.__init__(self, ZZ(dim), self._ring,
                                                  gap_group, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: groups.matrix.Heisenberg()
            Heisenberg group of degree 1 over Integer Ring
        """
        return "Heisenberg group of degree {} over {}".format(self._n, self._ring)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: H = groups.matrix.Heisenberg()
            sage: latex(H)
            H_{1}({\Bold{Z}})
        """
        return "H_{{{}}}({{{}}})".format(self._n, latex(self._ring))

    def order(self):
        """
        Return the order of ``self``.

        EXAMPLES::

            sage: H = groups.matrix.Heisenberg()
            sage: H.order()
            +Infinity
            sage: H = groups.matrix.Heisenberg(n=4)
            sage: H.order()
            +Infinity
            sage: H = groups.matrix.Heisenberg(R=3)
            sage: H.order()
            27
            sage: H = groups.matrix.Heisenberg(n=2, R=3)
            sage: H.order()
            243
            sage: H = groups.matrix.Heisenberg(n=2, R=GF(4))
            sage: H.order()
            1024
        """
        if self._ring is ZZ:
            from sage.rings.infinity import Infinity
            return Infinity
        else:
            return ZZ(self._ring.cardinality() ** (2*self._n + 1))

    cardinality = order

    def center(self):
        """
        Return the center of ``self``.

        This is the subgroup generated by the `z`, the matrix with a `1`
        in the upper right corner and along the diagonal.

        EXAMPLES::

            sage: H = groups.matrix.Heisenberg(2)
            sage: H.center()
            Subgroup with 1 generators (
            [1 0 0 1]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            ) of Heisenberg group of degree 2 over Integer Ring

            sage: H = groups.matrix.Heisenberg(3, 4)
            sage: H.center()
            Subgroup with 1 generators (
            [1 0 0 0 1]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            ) of Heisenberg group of degree 3 over Ring of integers modulo 4
        """
        return self.subgroup([self.gens()[-1]])
