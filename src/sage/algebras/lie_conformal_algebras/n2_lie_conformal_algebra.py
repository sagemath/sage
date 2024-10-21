# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules sage.rings.number_field
r"""
N=2 Super Lie Conformal Algebra

The `N=2` super Lie conformal algebra is an extension of the Virasoro
Lie conformal algebra (with generators `L,C`) by an even generator `J`
which is primary of conformal weight `1` and two odd generators
`G_1,G_2` which are primary of conformal weight `3/2`. The remaining
`\lambda`-brackets are given by:

.. MATH::

    [J_\lambda J] &= \frac{\lambda}{3} C, \\
    [J_\lambda G_1] &= G_1, \\
    [J_\lambda G_2] &= -G_2, \\
    [{G_1}_\lambda G_1] &= [{G_2}_\lambda G_2 ] = 0, \\
    [{G_1}_\lambda G_2] &= L + \frac{1}{2} TJ + \lambda J + \frac{\lambda^2}{6}C.

AUTHORS:

- Reimundo Heluani (2020-06-03): Initial implementation.
"""
# *****************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .graded_lie_conformal_algebra import GradedLieConformalAlgebra


class N2LieConformalAlgebra(GradedLieConformalAlgebra):
    """
    The N=2 super Lie conformal algebra.

    INPUT:

    - ``R`` -- a commutative ring; the base ring of this super
      Lie conformal algebra

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: F.<x> = NumberField(x^2 - 2)
        sage: R = lie_conformal_algebras.N2(F); R
        The N=2 super Lie conformal algebra over Number Field in x with defining polynomial x^2 - 2
        sage: R.inject_variables()
        Defining L, J, G1, G2, C
        sage: G1.bracket(G2)
        {0: L + 1/2*TJ, 1: J, 2: 1/3*C}
        sage: G2.bracket(G1)
        {0: L - 1/2*TJ, 1: -J, 2: 1/3*C}
        sage: G1.degree()
        3/2
        sage: J.degree()
        1

    The topological twist is a Virasoro vector with central
    charge 0::

        sage: L2 = L - 1/2*J.T()
        sage: L2.bracket(L2) == {0: L2.T(), 1: 2*L2}
        True

    The sum of the fermions is a generator of the Neveu-Schwarz
    Lie conformal algebra::

        sage: G = (G1 + G2)
        sage: G.bracket(G)
        {0: 2*L, 2: 2/3*C}
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: V = lie_conformal_algebras.N2(QQ)
            sage: TestSuite(V).run()
        """
        n2dict = {('L', 'L'): {0: {('L', 1): 1},
                               1: {('L', 0): 2},
                               3: {('C', 0): R(2).inverse_of_unit()}},
                  ('L', 'G1'): {0: {('G1', 1): 1},
                                1: {('G1', 0): 3 * R(2).inverse_of_unit()}},
                  ('L', 'G2'): {0: {('G2', 1): 1},
                                1: {('G2', 0): 3 * R(2).inverse_of_unit()}},
                  ('G1', 'G2'): {0: {('L', 0): 1, ('J', 1): R(2).inverse_of_unit()},
                                 1: {('J', 0): 1},
                                 2: {('C', 0): R(3).inverse_of_unit()}},
                  ('L', 'J'): {0: {('J', 1): 1}, 1: {('J', 0): 1}},
                  ('J', 'J'): {1: {('C', 0): R(3).inverse_of_unit()}},
                  ('J', 'G1'): {0: {('G1', 0): 1}},
                  ('J', 'G2'): {0: {('G2', 0): -1}}}
        from sage.rings.rational_field import QQ
        weights = (2, 1, QQ(3) / 2, QQ(3) / 2)
        parity = (0, 0, 1, 1)
        GradedLieConformalAlgebra.__init__(self, R, n2dict,
                                           names=('L', 'J', 'G1', 'G2'),
                                           central_elements=('C',),
                                           weights=weights, parity=parity)

    def _repr_(self):
        """
        The name of this Lie conformal algebra.

        EXAMPLES::

            sage: R = lie_conformal_algebras.N2(QQbar); R
            The N=2 super Lie conformal algebra over Algebraic Field
        """
        return f"The N=2 super Lie conformal algebra over {self.base_ring()}"
