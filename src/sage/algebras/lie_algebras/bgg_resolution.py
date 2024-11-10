r"""
BGG Resolutions

AUTHORS:

- Travis Scrimshaw (2024-01-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2024 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_function
from sage.structure.unique_representation import UniqueRepresentation
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.homology.chain_complex import ChainComplex_class


class BGGResolution(UniqueRepresentation, ChainComplex_class):
    r"""
    The BGG resolution of a simple module.

    We realize the BGG resolution as a chain complex, where the `(-1)`-th
    factor corresponds to  the finite dimensional simple module `L_{\lambda}`
    and the `i`-th factor (`i \geq 0`) corresponds to

    .. MATH::

        M_i := \bigoplus_{\substack{w \in W \\ \ell(w) = i}} M_{w\lambda}.

    Since the morphisms can be defined in terms of the image of the
    highest weight vectors, we only encode this information as a
    (finite) chain complex. We do not include the final natural projection
    map `p: M_{\lambda} \to L_{\lambda}` since the highest weight vector of
    weight `\lambda` only occurs in `M_{\lambda}` and `L_{\lambda}`.

    INPUT:

    - ``L`` -- a simple module

    EXAMPLES::

        sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
        sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
        sage: L = g.simple_module(La[1] + 4*La[2])
        sage: res = L.bgg_resolution()
        sage: ascii_art(res)
                                [ 1 -1]       [1]
                    [1 1]       [-1  1]       [1]
         0 <-- C_0 <------ C_1 <-------- C_2 <---- C_3 <-- 0

        sage: g = LieAlgebra(QQ, cartan_type=['D', 4])
        sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
        sage: L = g.simple_module(La[1] + La[2] + 3*La[3])
        sage: res = L.bgg_resolution()
        sage: w0 = WeylGroup(g.cartan_type(), prefix='s').long_element()
        sage: all(res.differential(i) * res.differential(i+1) == 0
        ....:     for i in range(w0.length()))
        True
    """
    def __init__(self, L):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[2])
            sage: res = L.bgg_resolution()
            sage: TestSuite(res).run()
        """
        from sage.combinat.root_system.weyl_group import WeylGroup
        ct = L.lie_algebra().cartan_type()
        self._cartan_type = ct
        self._simple = L
        self._W = WeylGroup(ct, prefix='s')
        # finish initialization
        R = self._simple.base_ring()
        differentials, mod_order = build_differentials(self._W)
        differentials = {deg: mat.change_ring(R) for deg, mat in differentials.items()}
        for deg in differentials:
            differentials[deg].set_immutable()
        self._module_order = mod_order
        super().__init__(ZZ, -ZZ.one(), R, differentials)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1])
            sage: L.bgg_resolution()
            BGG resolution of Simple module with highest weight Lambda[1]
             of Lie algebra of ['B', 2] in the Chevalley basis
        """
        return "BGG resolution of " + repr(self._simple)

    def simple_module(self):
        r"""
        Return the simple module `L_{\lambda}` defining ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['C', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[2])
            sage: res = L.bgg_resolution()
            sage: res.simple_module() is L
            True
        """
        return self._simple

    def module_order(self, i):
        r"""
        Return a tuple of Weyl group elements of length ``i``
        determining the ordering of the direct sum defining the
        differential matrix.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1])
            sage: res = L.bgg_resolution()
            sage: [res.module_order(i) for i in range(7)]
            [[1],
             [s2, s1],
             [s2*s1, s1*s2],
             [s1*s2*s1, s2*s1*s2],
             [s1*s2*s1*s2, s2*s1*s2*s1],
             [s1*s2*s1*s2*s1, s2*s1*s2*s1*s2],
             [s2*s1*s2*s1*s2*s1]]
        """
        if i not in self._module_order:
            return []
        return self._module_order[i]


@cached_function
def build_differentials(W):
    r"""
    Construct the differentials for the BGG resolution corresponding
    to the Weyl group `W`.

    ALGORITHM:

    We use the fact that (parabolic) Bruhat order is built locally
    from squares, all values defining the differential are `+1` or `-1`,
    and the product over the two different paths must sum to `0`.
    This is outlined in Ch. 6 of [Humphreys08]_.

    This only depends on the Coxeter group `W`. There is no stabilizer
    for any dominant integral weight `\lambda` undert the dot action
    (i.e., the stabilizer of `\lambda + \rho` is empty).

    EXAMPLES::

        sage: from sage.algebras.lie_algebras.bgg_resolution import build_differentials
        sage: W = WeylGroup(['B', 2], prefix='s')
        sage: D, O = build_differentials(W)
        sage: D
        {0: [],
         1: [-1  1],
         2: [1 1]
            [1 1],
         3: [ 1 -1]
            [-1  1],
         4: [1]
            [1],
         5: []}
        sage: O
        {0: [1],
         1: [s2, s1],
         2: [s2*s1, s1*s2],
         3: [s2*s1*s2, s1*s2*s1],
         4: [s2*s1*s2*s1]}
    """
    from itertools import combinations
    w0 = W.long_element()
    maxlen = w0.length()
    module_order = {i: [] for i in range(maxlen+1)}
    for w in W:
        module_order[w.length()].append(w)

    one = ZZ.one()
    # Set the initial step
    prev = {w: (j, frozenset([0])) for j, w in enumerate(module_order[maxlen-1])}
    prev_mat = matrix(ZZ, [[one]]*len(module_order[maxlen-1]), immutable=True)
    differentials = {maxlen: prev_mat}
    for i in range(maxlen-2, -1, -1):
        mat = matrix.zero(ZZ, len(module_order[i]), len(prev))
        cur = {}
        for j, w in enumerate(module_order[i]):
            # build the data to find the squares
            covers = frozenset([prev[v][0] for v in w.bruhat_upper_covers()])
            cur[w] = (j, covers)
            # set the indices in each square
            for v, vp in combinations(w.bruhat_upper_covers(), 2):
                # get the (unique) value at the top of the square
                dv = prev[v]
                dvp = prev[vp]
                vind = dv[0]
                vpind = dvp[0]
                for uind in dv[1] & dvp[1]:
                    # set the entries corresponding to the square
                    if not mat[j, vind]:
                        if not mat[j, vpind]:
                            mat[j, vpind] = one
                        mat[j, vind] = -mat[j, vpind] * prev_mat[vpind, uind] * prev_mat[vind, uind]
                    elif not mat[j, vpind]:
                        mat[j, vpind] = -mat[j, vind] * prev_mat[vpind, uind] * prev_mat[vind, uind]
                    else:
                        assert mat[j, vpind] * prev_mat[vpind, uind] + mat[j, vind] * prev_mat[vind, uind] == 0
        differentials[i+1] = mat
        prev = cur
        prev_mat = mat
    differentials[0] = matrix.zero(ZZ, 0, 1)
    differentials[maxlen+1] = matrix.zero(ZZ, 1, 0)
    return differentials, module_order
