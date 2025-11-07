# sage.doctest: needs sage.combinat sage.modules
r"""
Star-crystal structure on `B(\infty)`

AUTHORS:

- Ben Salisbury: initial version
- Travis Scrimshaw: initial version
"""

#*****************************************************************************
#       Copyright (C) 2016 Ben Salisbury <ben DOT salisbury AT cmich DOT edu>
#                          Travis Scrimshaw <tscrimsh AT umn DOT edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.combinat.crystals.elementary_crystals import ElementaryCrystal


class StarCrystal(UniqueRepresentation, Parent):
    r"""
    The star-crystal or `*`-crystal version of a highest weight crystal.

    The `*`-crystal structure on `B(\infty)` is the structure induced by
    the algebra antiautomorphism `* \colon U_q(\mathfrak{g}) \longrightarrow
    U_q(\mathfrak{g})` that stabilizes the negative half `U_q^-(\mathfrak{g})`.
    It is defined by

    .. MATH::

        E_i^* = E_i , \ \ \
        F_i^* = F_i , \ \ \
        q^* = q, \ \ \
        (q^h)^* = q^{-h},

    where `E_i` and `F_i` are the Chevalley generators of `U_q(\mathfrak{g})`
    and `h` is an element of the Cartan subalgebra.

    The induced operation on the crystal `B(\infty)` is called the
    *Kashiwara involution*.  Its implementation here is based on the
    recursive algorithm from Theorem 2.2.1 of [Ka1993]_, which states
    that for any `i \in I` there is a unique strict crystal embedding

    .. MATH::

        \Psi_i\colon B(\infty) \longrightarrow B_i \otimes B(\infty)

    such that

    - `u_{\infty} \mapsto b_i(0) \otimes u_{\infty}`, where `u_{\infty}`
      is the highest weight vector in `B(\infty)`;

    - if `\Psi_i(b) = f_i^mb_i(0) \otimes b_0`, then
      `\Psi_i(f_i^*b) = f_i^{m+1}b_i(0) \otimes b_0`
      and `\varepsilon_i(b^*) = m`;

    - the image of `\Psi_i` is `\{f_i^mb_i(0)\otimes b :
      \varepsilon_i(b^*) = 0, \ m\ge 0\}`.

    Here, `B_i` is the `i`-th elementary crystal.  See
    :class:`~sage.combinat.crystals.elementary_crystals.ElementaryCrystal`
    for more information.

    INPUT:

    - ``Binf`` -- a crystal from
      :class:`~sage.combinat.crystals.catalog_infinity_crystals`

    EXAMPLES::

        sage: B = crystals.infinity.Tableaux(['A',2])
        sage: Bstar = crystals.infinity.Star(B)
        sage: mg = Bstar.highest_weight_vector()
        sage: mg
        [[1, 1], [2]]
        sage: mg.f_string([1,2,1,2,2])
        [[1, 1, 1, 1, 1, 2, 2], [2, 3, 3, 3]]
    """

    def __init__(self, Binf):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['A',2])
            sage: Bstar = crystals.infinity.Star(B)
            sage: TestSuite(Bstar).run(max_runs=40)
            sage: TestSuite(Bstar).run(max_runs=1000) # long time
        """
        self._Binf = Binf
        self._cartan_type = Binf.cartan_type()
        Parent.__init__(self, category=HighestWeightCrystals().Infinite())
        self.module_generators = (self(self._Binf.module_generators[0]),)
        t0 = Binf.highest_weight_vector()
        B = {i: ElementaryCrystal(Binf.cartan_type(),i) for i in self.index_set()}
        self._tens = {i: B[i].tensor(Binf) for i in self.index_set()}
        gens = {i: self._tens[i](B[i](0), t0) for i in self.index_set()}
        self._embedding = {i: Binf.crystal_morphism({t0: gens[i]}) for i in self.index_set()}
        self._pullback = {i: self._tens[i].crystal_morphism({gens[i]: t0}) for i in self.index_set()}

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(3)
            sage: Ystar = crystals.infinity.Star(Y)
            sage: Ystar
            Star-crystal version of Crystal of generalized Young walls of type ['A', 3, 1]
        """
        return "Star-crystal version of %s" % self._Binf

    class Element(ElementWrapper):

        def e(self, i):
            r"""
            Return the action of `e_i^*` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: RC = crystals.infinity.RiggedConfigurations(['E',6,1])
                sage: RCstar = crystals.infinity.Star(RC)
                sage: nuJ = RCstar.module_generators[0].f_string([0,4,6,1,2])
                sage: ascii_art(nuJ.e(1))
                -1[ ]-1   (/)   0[ ]1   (/)   -1[ ]-1   (/)   -2[ ]-1

                sage: M = crystals.infinity.NakajimaMonomials(['B',2,1])
                sage: Mstar = crystals.infinity.Star(M)
                sage: m = Mstar.module_generators[0].f_string([0,1,2,2,1,0])
                sage: m.e(1)
                Y(0,0)^-1 Y(0,2)^-1 Y(1,1) Y(1,2)^-1 Y(2,1)^2
            """
            P = self.parent()
            image = P._embedding[i](self.value)
            if image[0].e(i)._m > 0:
                return None
            return P(P._pullback[i]( P._tens[i](image[0].e(i),image[1]) ))

        def f(self, i):
            r"""
            Return the action of `f_i^*` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: T = crystals.infinity.Tableaux("G2")
                sage: Tstar = crystals.infinity.Star(T)
                sage: t = Tstar.module_generators[0].f_string([1,2,1,1,2])
                sage: t
                [[1, 1, 1, 2, 0], [2, 3]]

                sage: M = crystals.infinity.NakajimaMonomials(['B',2,1])
                sage: Mstar = crystals.infinity.Star(M)
                sage: m = Mstar.module_generators[0].f_string([0,1,2,2,1,0])
                sage: m
                Y(0,0)^-1 Y(0,2)^-1 Y(1,0)^-1 Y(1,2)^-1 Y(2,0)^2 Y(2,1)^2
            """
            P = self.parent()
            image = P._embedding[i](self.value)
            return P(P._pullback[i]( P._tens[i](image[0].f(i),image[1]) ))

        def weight(self):
            r"""
            Return the weight of ``self``.

            EXAMPLES::

                sage: RC = crystals.infinity.RiggedConfigurations(['E',6,1])
                sage: RCstar = crystals.infinity.Star(RC)
                sage: nuJ = RCstar.module_generators[0].f_string([0,4,6,1,2])
                sage: nuJ.weight()
                -Lambda[0] - 2*Lambda[1] + 2*Lambda[3] - Lambda[4]
                 + 2*Lambda[5] - 2*Lambda[6] - delta
            """
            return self.value.weight()

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i^*` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: Y = crystals.infinity.GeneralizedYoungWalls(3)
                sage: Ystar = crystals.infinity.Star(Y)
                sage: y = Ystar.module_generators[0].f_string([0,1,3,2,1,0])
                sage: [y.epsilon(i) for i in y.index_set()]
                [1, 0, 1, 0]

                sage: RC = crystals.infinity.RiggedConfigurations(['E',6,1])
                sage: RCstar = crystals.infinity.Star(RC)
                sage: nuJ = RCstar.module_generators[0].f_string([0,4,6,1,2])
                sage: [nuJ.epsilon(i) for i in nuJ.index_set()]
                [0, 1, 1, 0, 0, 0, 1]
            """
            ep = -1
            while self is not None:
                ep += 1
                self = self.e(i)
            return ep

        def phi(self, i):
            r"""
            Return `\varphi_i^*` of ``self``.

            For `b \in B(\infty)`,

            .. MATH::

                \varphi_i^*(b) = \varepsilon_i^*(b) + \langle h_i,
                \mathrm{wt}(b) \rangle,

            where `h_i` is a simple coroot.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: T = crystals.infinity.Tableaux("A2")
                sage: Tstar = crystals.infinity.Star(T)
                sage: t = Tstar.module_generators[0].f_string([1,2,1,1,2])
                sage: [t.phi(i) for i in t.index_set()]
                [-3, 1]

                sage: M = crystals.infinity.NakajimaMonomials(['B',2,1])
                sage: Mstar = crystals.infinity.Star(M)
                sage: m = Mstar.module_generators[0].f_string([0,1,2,2,1,0])
                sage: [m.phi(i) for i in m.index_set()]
                [-1, -1, 4]
            """
            P = self.parent().weight_lattice_realization()
            ac = P.simple_coroot(i)
            return P(self.weight()).scalar(ac) + self.epsilon(i)

        def jump(self, i):
            r"""
            Return the `i`-jump of ``self``.

            For `b \in B(\infty)`,

            .. MATH::

                \operatorname{jump}_i(b) = \varepsilon_i(b) + \varepsilon_i^*(b)
                + \langle h_i, \mathrm{wt}(b) \rangle,

            where `h_i` is a simple coroot.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: RC = crystals.infinity.RiggedConfigurations("D4")
                sage: RCstar = crystals.infinity.Star(RC)
                sage: nu0star = RCstar.module_generators[0]
                sage: nustar = nu0star.f_string([2,1,3,4,2])
                sage: [nustar.jump(i) for i in RC.index_set()]
                [0, 1, 0, 0]
                sage: nustar = nu0star.f_string([2,1,3,4,2,2,1,3,2]) # long time
                sage: [nustar.jump(i) for i in RC.index_set()] # long time
                [1, 0, 1, 2]
            """
            P = self.parent().weight_lattice_realization()
            ac = P.simple_coroot(i)
            return P(self.value.weight()).scalar(ac) + self.epsilon(i) + self.value.epsilon(i)
