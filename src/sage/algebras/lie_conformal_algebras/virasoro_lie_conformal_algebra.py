# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Virasoro Lie Conformal Algebra

The Virasoro Lie conformal algebra is generated by `L` and a central
element `C`. The `\lambda`-brackets are given by:

.. MATH::

    [L_\lambda L] = T L + 2 \lambda L + \frac{\lambda^3}{12} C.

It is an H-graded Lie conformal algebra with `L` of degree `2`.

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .graded_lie_conformal_algebra import GradedLieConformalAlgebra


class VirasoroLieConformalAlgebra(GradedLieConformalAlgebra):
    """
    The Virasoro Lie Conformal algebra over `R`.

    INPUT:

    - ``R`` -- a commutative ring; behaviour is undefined if `R` is
      not a Field of characteristic zero.

    EXAMPLES::

        sage: Vir = lie_conformal_algebras.Virasoro(QQ)
        sage: Vir.category()
        Category of H-graded finitely generated Lie conformal algebras with basis over Rational Field
        sage: Vir.inject_variables()
        Defining L, C
        sage: L.bracket(L)
        {0: TL, 1: 2*L, 3: 1/2*C}

    TESTS::

        sage: Vir.gens()
        (L, C)
    """
    def __init__(self, R):
        """
        Initialize self.

        TESTS::

            sage: V = lie_conformal_algebras.Virasoro(QQ)
            sage: TestSuite(V).run()
        """
        virdict = {('L', 'L'): {0: {('L', 1): 1},
                                1: {('L', 0): 2},
                                3: {('C', 0): R(2).inverse_of_unit()}}}
        GradedLieConformalAlgebra.__init__(self, R, virdict,
            names=('L',), central_elements=('C',), weights=(2,))

    def _repr_(self):
        """
        The name of this Lie conformal algebra.

        EXAMPLES::

            sage: lie_conformal_algebras.Virasoro(QQbar)
            The Virasoro Lie conformal algebra over Algebraic Field
        """
        return "The Virasoro Lie conformal algebra over {}".format(
                                                            self.base_ring())
