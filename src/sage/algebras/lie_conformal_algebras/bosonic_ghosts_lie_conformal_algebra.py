# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Bosonic Ghosts Lie Conformal Algebra

The *Bosonic-ghosts* or `\beta-\gamma`-system Lie conformal algebra
with `2n` generators is the H-graded Lie conformal algebra generated
by `\beta_i, \gamma_i, i = 1,\ldots,n` and a central element `K`, with
non-vanishing `\lambda`-brackets:

.. MATH::

    [{\beta_i}_\lambda \gamma_j] = \delta_{ij} K.

The generators `\beta_i` have degree `1` while the generators `\gamma_i`
have degree `0`.

AUTHORS:

- Reimundo Heluani (2020-06-15): Initial implementation.
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

from sage.matrix.special import identity_matrix
from sage.structure.indexed_generators import standardize_names_index_set
from .graded_lie_conformal_algebra import GradedLieConformalAlgebra


class BosonicGhostsLieConformalAlgebra(GradedLieConformalAlgebra):
    r"""
    The Bosonic ghosts or `\beta-\gamma`-system Lie conformal
    algebra.

    INPUT:

    - ``R`` -- a commutative ring
    - ``ngens`` -- an even positive Integer (default: ``2``); the
      number of non-central generators of this Lie conformal
      algebra.
    - ``names`` -- list of ``str``; alternative names for the
      generators
    - ``index_set`` -- an enumerated set; an indexing set for the
      generators

    OUTPUT:

    The Bosonic Ghosts Lie conformal algebra with generators
    `\beta_i,\gamma_i, i=1,\ldots,n` and `K`, where `2n` is
    ``ngens``.

    EXAMPLES::

        sage: R = lie_conformal_algebras.BosonicGhosts(QQ); R
        The Bosonic ghosts Lie conformal algebra with generators (beta, gamma, K) over Rational Field
        sage: R.inject_variables(); beta.bracket(gamma)
        Defining beta, gamma, K
        {0: K}
        sage: beta.degree()
        1
        sage: gamma.degree()
        0

        sage: R = lie_conformal_algebras.BosonicGhosts(QQbar, ngens = 4, names = 'abcd'); R
        The Bosonic ghosts Lie conformal algebra with generators (a, b, c, d, K) over Algebraic Field
        sage: R.structure_coefficients()
        Finite family {('a', 'c'): ((0, K),),  ('b', 'd'): ((0, K),),  ('c', 'a'): ((0, -K),),  ('d', 'b'): ((0, -K),)}

    TESTS::

        sage: lie_conformal_algebras.BosonicGhosts(AA).category()
        Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Real Field
    """

    def __init__(self, R, ngens=2, names=None, index_set=None):
        """
        Initialize ``self``.

        TESTS::

            sage: V = lie_conformal_algebras.BosonicGhosts(QQ)
            sage: TestSuite(V).run()
        """
        from sage.rings.integer_ring import ZZ
        try:
            assert (ngens in ZZ and ngens > 0 and not ngens % 2)
        except AssertionError:
            raise ValueError("ngens should be an even positive integer, " +
                             "got {}".format(ngens))
        latex_names = None
        half = ngens // 2
        if (names is None) and (index_set is None):
            from sage.misc.defaults import variable_names as varnames
            from sage.misc.defaults import latex_variable_names as laxnames
            names = varnames(half, 'beta') + varnames(half, 'gamma')
            latex_names = tuple(laxnames(half, r'\beta') +
                                laxnames(half, r'\gamma')) + ('K',)

        names, index_set = standardize_names_index_set(names=names,
                                                      index_set=index_set,
                                                      ngens=ngens)
        A = identity_matrix(R, half)
        from sage.matrix.special import block_matrix
        gram_matrix = block_matrix([[R.zero(), A], [-A, R.zero()]])
        ghostsdict = {(i, j): {0: {('K', 0): gram_matrix[index_set.rank(i),
                                                         index_set.rank(j)]}}
                      for i in index_set for j in index_set}
        weights = (1,) * half + (0,) * half
        super().__init__(R,
                         ghostsdict, names=names,
                         latex_names=latex_names,
                         index_set=index_set,
                         weights=weights,
                         central_elements=('K',))

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: lie_conformal_algebras.BosonicGhosts(QQbar)
            The Bosonic ghosts Lie conformal algebra with generators (beta, gamma, K) over Algebraic Field
        """
        return "The Bosonic ghosts Lie conformal algebra with generators {} "\
               "over {}".format(self.gens(),self.base_ring())
