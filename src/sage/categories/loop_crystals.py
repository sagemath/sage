r"""
Loop Crystals
"""
# ****************************************************************************
#  Copyright (C) 2015   Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.tensor import TensorProductsCategory
from sage.categories.map import Map
from sage.graphs.dot2tex_utils import have_dot2tex
from sage.functions.other import ceil

class LoopCrystals(Category_singleton):
    r"""
    The category of `U_q'(\mathfrak{g})`-crystals, where `\mathfrak{g}`
    is of affine type.

    The category is called loop crystals as we can also consider them
    as crystals corresponding to the loop algebra `\mathfrak{g}_0[t]`,
    where `\mathfrak{g}_0` is the corresponding classical type.

    EXAMPLES::

        sage: from sage.categories.loop_crystals import LoopCrystals
        sage: C = LoopCrystals()
        sage: C
        Category of loop crystals
        sage: C.super_categories()
        [Category of crystals]
        sage: C.example()
        Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(1,1)

    TESTS::

        sage: TestSuite(C).run()
        sage: B = FiniteCrystals().example()
        sage: TestSuite(B).run()
    """
    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.loop_crystals import LoopCrystals
            sage: LoopCrystals().super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    def example(self, n = 3):
        """
        Return an example of Kirillov-Reshetikhin crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: from sage.categories.loop_crystals import LoopCrystals
            sage: B = LoopCrystals().example(); B
            Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(1,1)
        """
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal
        return KirillovReshetikhinCrystal(['A', n, 1], 1, 1)

    class ParentMethods:
        def weight_lattice_realization(self):
            """
            Return the weight lattice realization used to express weights
            of elements in ``self``.

            The default is to use the non-extended affine weight lattice.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C.weight_lattice_realization()
                Ambient space of the Root system of type ['A', 5]
                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.weight_lattice_realization()
                Weight lattice of the Root system of type ['A', 2, 1]
            """
            F = self.cartan_type().root_system()
            return F.weight_lattice(extended=False)

        def digraph(self, subset=None, index_set=None):
            r"""
            Return the :class:`DiGraph` associated to ``self``.

            INPUT:

            - ``subset`` -- (optional) a subset of vertices for
              which the digraph should be constructed

            - ``index_set`` -- (optional) the index set to draw arrows

            .. SEEALSO::

                :meth:`sage.categories.crystals.Crystals.ParentMethods.digraph`

            EXAMPLES::

                sage: C = crystals.KirillovReshetikhin(['D',4,1], 2, 1)
                sage: G = C.digraph()
                sage: G.latex_options()  # optional - dot2tex
                LaTeX options for Digraph on 29 vertices:
                {...'edge_options': <function ... at ...>...}
                sage: view(G, tightpage=True)  # optional - dot2tex graphviz, not tested (opens external window)
            """
            G = Crystals().parent_class.digraph(self, subset, index_set)
            if have_dot2tex():
                def eopt(u_v_label):
                    return {"backward": u_v_label[2] == 0}
                G.set_latex_options(edge_options=eopt)
            return G

# TODO: Should we make "regular" an axiom?
class RegularLoopCrystals(Category_singleton):
    r"""
    The category of regular `U_q'(\mathfrak{g})`-crystals, where
    `\mathfrak{g}` is of affine type.
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.loop_crystals import RegularLoopCrystals
            sage: RegularLoopCrystals().super_categories()
            [Category of regular crystals,
             Category of loop crystals]
        """
        return [RegularCrystals(), LoopCrystals()]

    class ElementMethods:
        def classical_weight(self):
            """
            Return the classical weight of ``self``.

            EXAMPLES::

                sage: R = RootSystem(['A',2,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
                sage: hw = LS.classically_highest_weight_vectors()
                sage: [(v.weight(), v.classical_weight()) for v in hw]
                [(-2*Lambda[0] + 2*Lambda[1], (2, 0, 0)),
                 (-Lambda[0] + Lambda[2], (1, 1, 0))]
            """
            CT = self.cartan_type().classical()
            I0 = CT.index_set()
            La = CT.root_system().ambient_space().fundamental_weights()
            return sum(La[i] * (self.phi(i) - self.epsilon(i)) for i in I0)

class KirillovReshetikhinCrystals(Category_singleton):
    """
    Category of Kirillov-Reshetikhin crystals.
    """
    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.loop_crystals import KirillovReshetikhinCrystals
            sage: KirillovReshetikhinCrystals().super_categories()
            [Category of finite regular loop crystals]
        """
        return [RegularLoopCrystals().Finite()]

    class ParentMethods:
        @abstract_method
        def r(self):
            r"""
            Return the value `r` in ``self`` written as `B^{r,s}`.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',3,1], 2,4)
                sage: K.r()
                2
            """

        @abstract_method
        def s(self):
            r"""
            Return the value `s` in ``self`` written as `B^{r,s}`.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',3,1], 2,4)
                sage: K.s()
                4
            """

        @abstract_method(optional=True)
        def classical_decomposition(self):
            """
            Return the classical decomposition of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',3,1], 2,2)
                sage: K.classical_decomposition()
                The crystal of tableaux of type ['A', 3] and shape(s) [[2, 2]]
            """

        @cached_method
        def classically_highest_weight_vectors(self):
            """
            Return the classically highest weight elements of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['E',6,1],1,1)
                sage: K.classically_highest_weight_vectors()
                ([(1,)],)
            """
            I0 = self.cartan_type().classical().index_set()
            return tuple([x for x in self if x.is_highest_weight(I0)])

        # TODO: This is duplicated in tensor product category
        def cardinality(self):
            """
            Return the cardinality of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['E',6,1], 1,1)
                sage: K.cardinality()
                27
                sage: K = crystals.KirillovReshetikhin(['C',6,1], 4,3)
                sage: K.cardinality()
                4736732
            """
            CWLR = self.cartan_type().classical().root_system().ambient_space()
            return sum(CWLR.weyl_dimension(mg.classical_weight())
                       for mg in self.classically_highest_weight_vectors())

        @cached_method
        def maximal_vector(self):
            r"""
            Return the unique element of classical weight `s \Lambda_r`
            in ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['C',2,1],1,2)
                sage: K.maximal_vector()
                [[1, 1]]
                sage: K = crystals.KirillovReshetikhin(['E',6,1],1,1)
                sage: K.maximal_vector()
                [(1,)]

                sage: K = crystals.KirillovReshetikhin(['D',4,1],2,1)
                sage: K.maximal_vector()
                [[1], [2]]

            TESTS:

            Check that :trac:`23028` is fixed::

                sage: ct = CartanType(['A',8,2]).dual()
                sage: K = crystals.KirillovReshetikhin(ct, 4, 1)
                sage: K.maximal_vector()
                [[1], [2], [3], [4]]
                sage: K = crystals.KirillovReshetikhin(ct, 4, 2)
                sage: K.maximal_vector()
                [[1, 1], [2, 2], [3, 3], [4, 4]]
            """
            R = self.weight_lattice_realization()
            Lambda = R.fundamental_weights()
            r = self.r()
            s = self.s()
            if self.cartan_type().dual().type() == 'BC':
                if self.cartan_type().rank() - 1 == r:
                    weight = 2*s*Lambda[r] - s*Lambda[0]
                else:
                    weight = s*Lambda[r] - s*Lambda[0]
            else:
                weight = s*Lambda[r] - s*Lambda[0] * Lambda[r].level() / Lambda[0].level()

            # First check the module generators as it is likely to be in here
            for b in self.module_generators:
                if b.weight() == weight:
                    return b

            # Otherwise check all of the elements
            for b in self:
                if b not in self.module_generators and b.weight() == weight:
                    return b

            assert False, "BUG: invalid Kirillov-Reshetikhin crystal"

        def module_generator(self):
            r"""
            Return the unique module generator of classical weight
            `s \Lambda_r` of the Kirillov-Reshetikhin crystal `B^{r,s}`.

            EXAMPLES::

                sage: La = RootSystem(['G',2,1]).weight_space().fundamental_weights()
                sage: K = crystals.ProjectedLevelZeroLSPaths(La[1])
                sage: K.module_generator()
                (-Lambda[0] + Lambda[1],)
            """
            return self.maximal_vector()

        # TODO: Should this be in one of the super categories?
        def affinization(self):
            """
            Return the corresponding affinization crystal of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.affinization()
                Affinization of Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1)

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1, model='KR')
                sage: K.affinization()
                Affinization of Kirillov-Reshetikhin tableaux of type ['A', 2, 1] and shape (1, 1)
            """
            from sage.combinat.crystals.affinization import AffinizationOfCrystal
            return AffinizationOfCrystal(self)

        @cached_method
        def R_matrix(self, K):
            r"""
            Return the combinatorial `R`-matrix of ``self`` to ``K``.

            The *combinatorial* `R`-*matrix* is the affine crystal
            isomorphism `R : L \otimes K \to K \otimes L` which maps
            `u_{L} \otimes u_K` to `u_K \otimes u_{L}`, where `u_K`
            is the unique element in `K = B^{r,s}` of weight
            `s\Lambda_r - s c \Lambda_0` (see :meth:`maximal_vector`).

            INPUT:

            - ``self`` -- a crystal `L`
            - ``K`` -- a Kirillov-Reshetikhin crystal of the same type as `L`

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                sage: L = crystals.KirillovReshetikhin(['A',2,1],1,2)
                sage: f = K.R_matrix(L)
                sage: [[b,f(b)] for b in crystals.TensorProduct(K,L)]
                [[[[[1]], [[1, 1]]], [[[1, 1]], [[1]]]],
                 [[[[1]], [[1, 2]]], [[[1, 1]], [[2]]]],
                 [[[[1]], [[2, 2]]], [[[1, 2]], [[2]]]],
                 [[[[1]], [[1, 3]]], [[[1, 1]], [[3]]]],
                 [[[[1]], [[2, 3]]], [[[1, 2]], [[3]]]],
                 [[[[1]], [[3, 3]]], [[[1, 3]], [[3]]]],
                 [[[[2]], [[1, 1]]], [[[1, 2]], [[1]]]],
                 [[[[2]], [[1, 2]]], [[[2, 2]], [[1]]]],
                 [[[[2]], [[2, 2]]], [[[2, 2]], [[2]]]],
                 [[[[2]], [[1, 3]]], [[[2, 3]], [[1]]]],
                 [[[[2]], [[2, 3]]], [[[2, 2]], [[3]]]],
                 [[[[2]], [[3, 3]]], [[[2, 3]], [[3]]]],
                 [[[[3]], [[1, 1]]], [[[1, 3]], [[1]]]],
                 [[[[3]], [[1, 2]]], [[[1, 3]], [[2]]]],
                 [[[[3]], [[2, 2]]], [[[2, 3]], [[2]]]],
                 [[[[3]], [[1, 3]]], [[[3, 3]], [[1]]]],
                 [[[[3]], [[2, 3]]], [[[3, 3]], [[2]]]],
                 [[[[3]], [[3, 3]]], [[[3, 3]], [[3]]]]]

                sage: K = crystals.KirillovReshetikhin(['D',4,1],1,1)
                sage: L = crystals.KirillovReshetikhin(['D',4,1],2,1)
                sage: f = K.R_matrix(L)
                sage: T = crystals.TensorProduct(K,L)
                sage: b = T( K(rows=[[1]]), L(rows=[]) )
                sage: f(b)
                [[[2], [-2]], [[1]]]

            Alternatively, one can compute the combinatorial `R`-matrix
            using the isomorphism method of digraphs::

                sage: K1 = crystals.KirillovReshetikhin(['A',2,1],1,1)
                sage: K2 = crystals.KirillovReshetikhin(['A',2,1],2,1)
                sage: T1 = crystals.TensorProduct(K1,K2)
                sage: T2 = crystals.TensorProduct(K2,K1)
                sage: T1.digraph().is_isomorphic(T2.digraph(), edge_labels=True, certificate=True) #todo: not implemented (see #10904 and #10549)
                (True, {[[[1]], [[2], [3]]]: [[[1], [3]], [[2]]], [[[3]], [[2], [3]]]: [[[2], [3]], [[3]]],
                [[[3]], [[1], [3]]]: [[[1], [3]], [[3]]], [[[1]], [[1], [3]]]: [[[1], [3]], [[1]]], [[[1]],
                [[1], [2]]]: [[[1], [2]], [[1]]], [[[2]], [[1], [2]]]: [[[1], [2]], [[2]]], [[[3]],
                [[1], [2]]]: [[[2], [3]], [[1]]], [[[2]], [[1], [3]]]: [[[1], [2]], [[3]]], [[[2]], [[2], [3]]]: [[[2], [3]], [[2]]]})
            """
            from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
            T1 = TensorProductOfCrystals(self, K)
            T2 = TensorProductOfCrystals(K, self)
            gen1 = T1(self.maximal_vector(), K.maximal_vector())
            gen2 = T2(K.maximal_vector(), self.maximal_vector())
            return T1.crystal_morphism({gen1: gen2}, check=False)

        @cached_method
        def local_energy_function(self, B):
            r"""
            Return the local energy function of ``self`` and ``B``.

            See
            :class:`~sage.categories.loop_crystals.LocalEnergyFunction`
            for a definition.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',6,2], 2,1)
                sage: Kp = crystals.KirillovReshetikhin(['A',6,2], 1,1)
                sage: H = K.local_energy_function(Kp); H
                Local energy function of
                 Kirillov-Reshetikhin crystal of type ['BC', 3, 2] with (r,s)=(2,1)
                tensor
                 Kirillov-Reshetikhin crystal of type ['BC', 3, 2] with (r,s)=(1,1)
            """
            return LocalEnergyFunction(self, B)

        @cached_method
        def b_sharp(self):
            r"""
            Return the element `b^{\sharp}` of ``self``.

            Let `B` be a KR crystal. The element `b^{\sharp}` is the unique
            element such that `\varphi(b^{\sharp}) = \ell \Lambda_0` with
            `\ell = \min \{ \langle c, \varphi(b) \rangle \mid b \in B \}`.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',6,2], 2,1)
                sage: K.b_sharp()
                []
                sage: K.b_sharp().Phi()
                Lambda[0]

                sage: K = crystals.KirillovReshetikhin(['C',3,1], 1,3)
                sage: K.b_sharp()
                [[-1]]
                sage: K.b_sharp().Phi()
                2*Lambda[0]

                sage: K = crystals.KirillovReshetikhin(['D',6,2], 2,2)
                sage: K.b_sharp() # long time
                []
                sage: K.b_sharp().Phi() # long time
                2*Lambda[0]
            """
            ell = float('inf')
            bsharp = None
            for b in self:
                phi = b.Phi()
                if phi.support() == [0] and phi[0] < ell:
                    bsharp = b
                    ell = phi[0]
            return bsharp

        def is_perfect(self, ell=None):
            r"""
            Check if ``self`` is a perfect crystal of level ``ell``.

            A crystal `\mathcal{B}` is perfect of level `\ell` if:

            #. `\mathcal{B}` is isomorphic to the crystal graph of a
               finite-dimensional `U_q'(\mathfrak{g})`-module.
            #. `\mathcal{B} \otimes \mathcal{B}` is connected.
            #. There exists a `\lambda\in X`, such that
               `\mathrm{wt}(\mathcal{B}) \subset \lambda + \sum_{i\in I}
               \ZZ_{\le 0} \alpha_i` and there is a unique element in
               `\mathcal{B}` of classical weight `\lambda`.
            #. For all `b \in \mathcal{B}`,
               `\mathrm{level}(\varepsilon (b)) \geq \ell`.
            #. For all `\Lambda` dominant weights of level `\ell`, there
               exist unique elements `b_{\Lambda}, b^{\Lambda} \in
               \mathcal{B}`, such that `\varepsilon(b_{\Lambda}) =
               \Lambda = \varphi(b^{\Lambda})`.

            Points (1)-(3) are known to hold. This method checks
            points (4) and (5).

            If ``self`` is the Kirillov-Reshetikhin crystal `B^{r,s}`,
            then it was proven for non-exceptional types in [FOS2010]_
            that it is perfect if and only if `s/c_r` is an integer
            (where `c_r` is a constant related to the type of the crystal).

            It is conjectured this is true for all affine types.

            INPUT:

            - ``ell`` -- (default: `s / c_r`) integer; the level

            REFERENCES:

            [FOS2010]_

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.is_perfect()
                True

                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 1)
                sage: K.is_perfect()
                False

                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
                sage: K.is_perfect()
                True

                sage: K = crystals.KirillovReshetikhin(['E',6,1], 1,3)
                sage: K.is_perfect()
                True

            TESTS:

            Check that this works correctly for `B^{n,s}`
            of type `A_{2n}^{(2)\dagger}` (:trac:`24364`)::

                sage: K = crystals.KirillovReshetikhin(CartanType(['A',6,2]).dual(), 3,1)
                sage: K.is_perfect()
                True
                sage: K.is_perfect(1)
                True

            .. TODO::

                Implement a version for tensor products of KR crystals.
            """
            from sage.rings.integer_ring import ZZ
            if ell is None:
                if (self.cartan_type().dual().type() == 'BC'
                    and self.cartan_type().rank() - 1 == self.r()):
                    return True
                ell = self.s() / self.cartan_type().c()[self.r()]
                if ell not in ZZ:
                    return False

            if ell not in ZZ:
                raise ValueError("perfectness not defined for non-integral levels")

            # [FOS2010]_ check
            if self.cartan_type().classical().type() not in ['E','F','G']:
                if (self.cartan_type().dual().type() == 'BC'
                    and self.cartan_type().rank() - 1 == self.r()):
                    return ell == self.s()
                return ell == self.s() / self.cartan_type().c()[self.r()]

            # Check by definition
            # TODO: This is duplicated from ProjectedLevelZeroLSPaths, combine the two methods.
            # TODO: Similarly, don't duplicate in the tensor product category, maybe
            #   move this to the derived affine category?
            MPhi = []
            for b in self:
                p = b.Phi().level()
                assert p == b.Epsilon().level()
                if p < ell:
                    return False
                if p == ell:
                    MPhi += [b]
            weights = []
            I = self.index_set()
            rank = len(I)
            La = self.weight_lattice_realization().basis()
            from sage.combinat.integer_vector import IntegerVectors
            for n in range(1, ell+1):
                for c in IntegerVectors(n, rank):
                    w = sum(c[i]*La[i] for i in I)
                    if w.level() == ell:
                        weights.append(w)
            return sorted(b.Phi() for b in MPhi) == sorted(weights)

        def level(self):
            r"""
            Return the level of ``self`` when ``self`` is a perfect crystal.

            .. SEEALSO::

                :meth:`~sage.categories.loop_crystals.KirillovReshetikhinCrystals.ParentMethods.is_perfect`

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.level()
                1
                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
                sage: K.level()
                1
                sage: K = crystals.KirillovReshetikhin(['D',4,1], 1, 3)
                sage: K.level()
                3

                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 1)
                sage: K.level()
                Traceback (most recent call last):
                ...
                ValueError: this crystal is not perfect

            TESTS:

            Check that this works correctly for `B^{n,s}`
            of type `A_{2n}^{(2)\dagger}` (:trac:`24364`)::

                sage: ct = CartanType(['A',6,2]).dual()
                sage: K1 = crystals.KirillovReshetikhin(ct, 3,1)
                sage: K1.level()
                1
                sage: K4 = crystals.KirillovReshetikhin(ct, 3,4)
                sage: K4.level()
                4
            """
            if not self.is_perfect():
                raise ValueError("this crystal is not perfect")
            if (self.cartan_type().dual().type() == 'BC'
                and self.cartan_type().rank() - 1 == self.r()):
                return self.s()
            return self.s() / self.cartan_type().c()[self.r()]

        def q_dimension(self, q=None, prec=None, use_product=False):
            """
            Return the `q`-dimension of ``self``.

            The `q`-dimension of a KR crystal is defined as the `q`-dimension of
            the underlying classical crystal.

            EXAMPLES::

                sage: KRC = crystals.KirillovReshetikhin(['A',2,1], 2,2)
                sage: KRC.q_dimension()
                q^4 + q^3 + 2*q^2 + q + 1
                sage: KRC = crystals.KirillovReshetikhin(['D',4,1], 2,1)
                sage: KRC.q_dimension()
                q^10 + q^9 + 3*q^8 + 3*q^7 + 4*q^6 + 4*q^5 + 4*q^4 + 3*q^3 + 3*q^2 + q + 2
            """
            return self.classical_decomposition().q_dimension(q, prec, use_product)

    class ElementMethods:
        def lusztig_involution(self):
            r"""
            Return the result of the classical Lusztig involution on ``self``.

            EXAMPLES::

                sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KR')
                sage: mg = KRT.module_generators[1]
                sage: mg.lusztig_involution()
                [[-2, -2, 1], [-1, -1, 2]]
                sage: elt = mg.f_string([2,1,3,2]); elt
                [[3, -2, 1], [4, -1, 2]]
                sage: elt.lusztig_involution()
                [[-4, -2, 1], [-3, -1, 2]]
            """
            Cl = self.parent().cartan_type().classical()
            I = Cl.index_set()
            aut = Cl.opposition_automorphism()
            hw = self.to_highest_weight(I)[1]
            hw.reverse()
            return self.to_lowest_weight(I)[0].e_string(aut[i] for i in hw)

        @cached_method
        def energy_function(self):
            r"""
            Return the energy function of ``self``.

            Let `B` be a KR crystal. Let `b^{\sharp}` denote the unique
            element such that `\varphi(b^{\sharp}) = \ell \Lambda_0` with
            `\ell = \min \{ \langle c, \varphi(b) \mid b \in B \}`. Let
            `u_B` denote the maximal element of `B`. The *energy* of
            `b \in B` is given by

            .. MATH::

                D(b) = H(b \otimes b^{\sharp}) - H(u_B \otimes b^{\sharp}),

            where `H` is the :meth:`local energy function
            <sage.categories.loop_crystals.KirillovReshetikhinCrystals.ParentMethods.local_energy_function>`.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1)
                sage: for x in K.classically_highest_weight_vectors():
                ....:    x, x.energy_function()
                ([], 1)
                ([[1], [2]], 0)

                sage: K = crystals.KirillovReshetikhin(['D',4,3], 1,2)
                sage: for x in K.classically_highest_weight_vectors():
                ....:    x, x.energy_function()
                ([], 2)
                ([[1]], 1)
                ([[1, 1]], 0)
            """
            B = self.parent()
            bsharp = B.b_sharp()
            T = B.tensor(B)
            H = B.local_energy_function(B)
            return H(T(self, bsharp)) - H(T(B.maximal_vector(), bsharp))

    class TensorProducts(TensorProductsCategory):
        """
        The category of tensor products of Kirillov-Reshetikhin crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: from sage.categories.loop_crystals import KirillovReshetikhinCrystals
                sage: KirillovReshetikhinCrystals().TensorProducts().extra_super_categories()
                [Category of finite regular loop crystals]
            """
            return [RegularLoopCrystals().Finite()]

        class ParentMethods:
            @cached_method
            def maximal_vector(self):
                """
                Return the maximal vector of ``self``.

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: T.maximal_vector()
                    [[[1]], [[1]], [[1]]]
                """
                return self(*[K.maximal_vector() for K in self.crystals])

            @cached_method
            def classically_highest_weight_vectors(self):
                r"""
                Return the classically highest weight elements of ``self``.

                This works by using a backtracking algorithm since if
                `b_2 \otimes b_1` is classically highest weight then `b_1`
                is classically highest weight.

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: T.classically_highest_weight_vectors()
                    ([[[1]], [[1]], [[1]]],
                     [[[2]], [[1]], [[1]]],
                     [[[1]], [[2]], [[1]]],
                     [[[3]], [[2]], [[1]]])
                """
                n = len(self.crystals)
                I0 = self.cartan_type().classical().index_set()
                it = [ iter(self.crystals[-1].classically_highest_weight_vectors()) ]
                path = []
                ret = []
                while it:
                    try:
                        x = next(it[-1])
                    except StopIteration:
                        it.pop()
                        if path:
                            path.pop(0)
                        continue

                    b = self.element_class(self, [x] + path)
                    if not b.is_highest_weight(index_set=I0):
                        continue
                    path.insert(0, x)
                    if len(path) == n:
                        ret.append(b)
                        path.pop(0)
                    else:
                        it.append( iter(self.crystals[-len(path)-1]) )
                return tuple(ret)

            # TODO: This is duplicated in KR crystals category
            def cardinality(self):
                """
                Return the cardinality of ``self``.

                EXAMPLES::

                    sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
                    sage: RC.cardinality()
                    100
                    sage: len(RC.list())
                    100

                    sage: RC = RiggedConfigurations(['E', 7, 1], [[1,1]])
                    sage: RC.cardinality()
                    134
                    sage: len(RC.list())
                    134

                    sage: RC = RiggedConfigurations(['B', 3, 1], [[2,2],[1,2]])
                    sage: RC.cardinality()
                    5130
                """
                CWLR = self.cartan_type().classical().root_system().ambient_space()
                return sum(CWLR.weyl_dimension(mg.classical_weight())
                           for mg in self.classically_highest_weight_vectors())

            def one_dimensional_configuration_sum(self, q=None, group_components=True):
                r"""
                Compute the one-dimensional configuration sum of ``self``.

                INPUT:

                - ``q`` -- (default: ``None``) a variable or ``None``;
                  if ``None``, a variable `q` is set in the code
                - ``group_components`` -- (default: ``True``) boolean; if
                  ``True``, then the terms are grouped by classical component

                The one-dimensional configuration sum is the sum of the
                weights of all elements in the crystal weighted by the
                energy function.

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: T.one_dimensional_configuration_sum()
                    B[-2*Lambda[1] + 2*Lambda[2]] + (q+1)*B[-Lambda[1]]
                     + (q+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
                     + B[-2*Lambda[2]] + (q+1)*B[Lambda[2]]
                    sage: R.<t> = ZZ[]
                    sage: T.one_dimensional_configuration_sum(t, False)
                    B[-2*Lambda[1] + 2*Lambda[2]] + (t+1)*B[-Lambda[1]]
                     + (t+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
                     + B[-2*Lambda[2]] + (t+1)*B[Lambda[2]]

                    sage: R = RootSystem(['A',2,1])
                    sage: La = R.weight_space().basis()
                    sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
                    sage: LS.one_dimensional_configuration_sum() == T.one_dimensional_configuration_sum() # long time
                    True

                TESTS::

                    sage: K1 = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: K2 = crystals.KirillovReshetikhin(['A',2,1],2,1)
                    sage: T = crystals.TensorProduct(K1,K2)
                    sage: T.one_dimensional_configuration_sum() == T.one_dimensional_configuration_sum(group_components=False)
                    True

                    sage: RC = RiggedConfigurations(['A',3,1],[[1,1],[1,2]])
                    sage: B = crystals.KirillovReshetikhin(['A',3,1],1,1)
                    sage: B1 = crystals.KirillovReshetikhin(['A',3,1],1,2)
                    sage: T = crystals.TensorProduct(B,B1)
                    sage: RC.fermionic_formula() == T.one_dimensional_configuration_sum()
                    True
                """
                if q is None:
                    from sage.rings.rational_field import QQ
                    q = QQ['q'].gens()[0]
                P0 = self.weight_lattice_realization().classical()
                B = P0.algebra(q.parent())
                if group_components:
                    G = self.digraph(index_set=self.cartan_type().classical().index_set())
                    C = G.connected_components()
                    return B.sum(q**(c[0].energy_function())*B.sum(B(P0(b.weight())) for b in c)
                                 for c in C)
                return B.sum(q**(b.energy_function())*B(P0(b.weight())) for b in self)

        class ElementMethods:
            def energy_function(self, algorithm=None):
                r"""
                Return the energy function of ``self``.

                ALGORITHM:

                .. RUBRIC:: definition

                Let `T` be a tensor product of Kirillov-Reshetikhin
                crystals. Let `R_i` and `H_i` be the combinatorial
                `R`-matrix and local energy functions, respectively, acting
                on the `i` and `i+1` factors. Let `D_B` be the energy
                function of a single Kirillov-Reshetikhin crystal. The
                *energy function* is given by

                .. MATH::

                    D = \sum_{j > i} H_i R_{i+1} R_{i+2} \cdots R_{j-1}
                    + \sum_j D_B R_1 R_2 \cdots R_{j-1},

                where `D_B` acts on the rightmost factor.

                .. RUBRIC:: grading

                If  ``self`` is an element of `T`, a tensor product of
                perfect crystals of the same level, then use the affine
                grading to determine the energy. Specifically, let `g`
                denote the affine grading of ``self`` and `d` the affine
                grading of the maximal vector in `T`. Then the energy
                of ``self`` is given by `d - g`.

                For more details, see Theorem 7.5 in [ST2011]_.

                INPUT:

                - ``algorithm`` -- (default: ``None``) use one of the
                  following algorithms to determine the energy function:

                  * ``'definition'`` - use the definition of the energy
                    function;
                  * ``'grading'`` - use the affine grading;

                  if not specified, then this uses ``'grading'`` if all
                  factors are perfect of the same level and otherwise
                  this uses ``'definition'``

                OUTPUT: an integer

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = T.classically_highest_weight_vectors()
                    sage: for b in hw:
                    ....:     print("{} {}".format(b, b.energy_function()))
                    [[[1]], [[1]], [[1]]] 0
                    [[[2]], [[1]], [[1]]] 1
                    [[[1]], [[2]], [[1]]] 2
                    [[[3]], [[2]], [[1]]] 3

                    sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: hw = T.classically_highest_weight_vectors()
                    sage: for b in hw:
                    ....:     print("{} {}".format(b, b.energy_function()))
                    [[], []] 4
                    [[[1, 1]], []] 3
                    [[], [[1, 1]]] 1
                    [[[1, 1]], [[1, 1]]] 0
                    [[[1, 2]], [[1, 1]]] 1
                    [[[2, 2]], [[1, 1]]] 2
                    [[[-1, -1]], [[1, 1]]] 2
                    [[[1, -1]], [[1, 1]]] 2
                    [[[2, -1]], [[1, 1]]] 2

                    sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 1)
                    sage: T = crystals.TensorProduct(K)
                    sage: t = T.module_generators[0]
                    sage: t.energy_function('grading')
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: all crystals in the tensor product
                     need to be perfect of the same level

                TESTS::

                    sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
                    sage: K2 = crystals.KirillovReshetikhin(['C',2,1], 2, 2)
                    sage: T = tensor([K, K2])
                    sage: hw = T.classically_highest_weight_vectors()
                    sage: all(b.energy_function() == b.energy_function(algorithm='definition')
                    ....:     for b in hw)
                    True
                """
                C = self.parent().crystals[0]
                ell = ceil(C.s()/C.cartan_type().c()[C.r()])
                is_perfect = all(ell == K.s()/K.cartan_type().c()[K.r()]
                                 for K in self.parent().crystals)
                if algorithm is None:
                    if is_perfect:
                        algorithm = 'grading'
                    else:
                        algorithm = 'definition'

                if algorithm == 'grading':
                    if not is_perfect:
                        raise NotImplementedError("all crystals in the tensor product need to be perfect of the same level")
                    d = self.parent().maximal_vector().affine_grading()
                    return d - self.affine_grading()

                if algorithm == 'definition':
                    # Setup
                    from sage.rings.integer_ring import ZZ
                    energy = ZZ.zero()
                    R_mats = [[K.R_matrix(Kp) for Kp in self.parent().crystals[i+1:]]
                              for i,K in enumerate(self.parent().crystals)]
                    H_funcs = [[K.local_energy_function(Kp) for Kp in self.parent().crystals[i+1:]]
                               for i,K in enumerate(self.parent().crystals)]

                    for i,b in enumerate(self):
                        for j,R in enumerate(R_mats[i]):
                            H = H_funcs[i][j]
                            bp = self[i+j+1]
                            T = R.domain()
                            t = T(b, bp)
                            energy += H(t)
                            b = R(t)[1]
                        energy += b.energy_function()  # D contribution
                    return energy
                else:
                    raise ValueError("invalid algorithm")

            def affine_grading(self):
                r"""
                Return the affine grading of ``self``.

                The affine grading is calculated by finding a path
                from ``self`` to a ground state path (using the helper method
                :meth:`e_string_to_ground_state`) and counting the number
                of affine Kashiwara operators `e_0` applied on the way.

                OUTPUT: an integer

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]
                    sage: t.affine_grading()
                    1

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = T.classically_highest_weight_vectors()
                    sage: for b in hw:
                    ....:     print("{} {}".format(b, b.affine_grading()))
                    [[[1]], [[1]], [[1]]] 3
                    [[[2]], [[1]], [[1]]] 2
                    [[[1]], [[2]], [[1]]] 1
                    [[[3]], [[2]], [[1]]] 0

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = T.classically_highest_weight_vectors()
                    sage: for b in hw:
                    ....:     print("{} {}".format(b, b.affine_grading()))
                    [[[1]], [[1]], [[1]]] 2
                    [[[2]], [[1]], [[1]]] 1
                    [[[-1]], [[1]], [[1]]] 1
                    [[[1]], [[2]], [[1]]] 1
                    [[[-2]], [[2]], [[1]]] 0
                    [[[1]], [[-1]], [[1]]] 0
                """
                return self.e_string_to_ground_state().count(0)

            @cached_method
            def e_string_to_ground_state(self):
                r"""
                Return a string of integers in the index set
                `(i_1, \ldots, i_k)` such that `e_{i_k} \cdots e_{i_1}`
                of ``self`` is the ground state.

                This method calculates a path from ``self`` to a ground
                state path using Demazure arrows as defined in Lemma 7.3
                in [ST2011]_.

                OUTPUT: a tuple of integers `(i_1, \ldots, i_k)`

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]
                    sage: t.e_string_to_ground_state()
                    (0, 2)

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]; t
                    [[[1]], [[1]]]
                    sage: t.e_string_to_ground_state()
                    (0,)
                    sage: x = t.e(0)
                    sage: x.e_string_to_ground_state()
                    ()
                    sage: y = t.f_string([1,2,1,1,0]); y
                    [[[2]], [[1]]]
                    sage: y.e_string_to_ground_state()
                    ()

                TESTS:

                Check that :trac:`22882` is fixed::

                    sage: K = crystals.KirillovReshetikhin(CartanType(['A',6,2]).dual(), 1,1)
                    sage: T = tensor([K,K,K])
                    sage: hw = [x for x in T if x.is_highest_weight([1,2,3])]
                    sage: gs = T(K(0), K(0), K(0))
                    sage: all(elt.e_string(elt.e_string_to_ground_state()) == gs
                    ....:     for elt in hw)
                    True
                    sage: all(elt.energy_function() == elt.energy_function('definition')
                    ....:     for elt in hw)
                    True
                """
                ell = max(ceil(K.s()/K.cartan_type().c()[K.r()])
                          for K in self.parent().crystals)
                if self.cartan_type().dual().type() == 'BC':
                    I = self.cartan_type().index_set()
                    for i in I[:-1]:
                        if self.epsilon(i) > 0:
                            return (i,) + (self.e(i)).e_string_to_ground_state()
                    if self.epsilon(I[-1]) > ell:
                        return (I[-1],) + (self.e(I[-1])).e_string_to_ground_state()
                    return ()

                I = self.cartan_type().classical().index_set()
                for i in I:
                    if self.epsilon(i) > 0:
                        return (i,) + self.e(i).e_string_to_ground_state()
                if self.epsilon(0) > ell:
                    return (0,) + self.e(0).e_string_to_ground_state()
                return ()


#####################################################################
## Local energy function

class LocalEnergyFunction(Map):
    r"""
    The local energy function.

    Let `B` and `B'` be Kirillov-Reshetikhin crystals with maximal
    vectors `u_B` and `u_{B'}` respectively. The *local energy function*
    `H : B \otimes B' \to \ZZ` is the function which satisfies

    .. MATH::

        H(e_0(b \otimes b')) = H(b \otimes b') + \begin{cases}
        1 & \text{if } i = 0 \text{ and LL}, \\
        -1 & \text{if } i = 0 \text{ and RR}, \\
        0 & \text{otherwise,}
        \end{cases}

    where LL (resp. RR) denote `e_0` acts on the left (resp. right)
    on both `b \otimes b'` and `R(b \otimes b')`, and
    normalized by `H(u_B \otimes u_{B'}) = 0`.

    INPUT:

    - ``B`` -- a Kirillov-Reshetikhin crystal
    - ``Bp`` -- a Kirillov-Reshetikhin crystal
    - ``normalization`` -- (default: 0) the normalization value

    EXAMPLES::

        sage: K = crystals.KirillovReshetikhin(['C',2,1], 1,2)
        sage: K2 = crystals.KirillovReshetikhin(['C',2,1], 2,1)
        sage: H = K.local_energy_function(K2)
        sage: T = tensor([K, K2])
        sage: hw = T.classically_highest_weight_vectors()
        sage: for b in hw:
        ....:     b, H(b)
        ([[], [[1], [2]]], 1)
        ([[[1, 1]], [[1], [2]]], 0)
        ([[[2, -2]], [[1], [2]]], 1)
        ([[[1, -2]], [[1], [2]]], 1)

    REFERENCES:

    [KKMMNN1992]_
    """
    def __init__(self, B, Bp, normalization=0):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A',7,2], 1,2)
            sage: K2 = crystals.KirillovReshetikhin(['A',7,2], 2,1)
            sage: H = K.local_energy_function(K2)
            sage: TestSuite(H).run(skip=['_test_category', '_test_pickling'])

        TESTS:

        Check that :trac:`23014` is fixed::

            sage: La = RootSystem(['G',2,1]).weight_space().fundamental_weights()
            sage: K = crystals.ProjectedLevelZeroLSPaths(La[1])
            sage: H = K.local_energy_function(K)
            sage: hw = H.domain().classically_highest_weight_vectors()
            sage: [H(x) for x in hw]
            [0, 1, 2, 1]
        """
        from sage.rings.integer_ring import ZZ
        self._B = B
        self._Bp = Bp
        self._R_matrix = self._B.R_matrix(self._Bp)
        T = B.tensor(Bp)
        self._known_values = {T(*[K.maximal_vector() for K in T.crystals]):
                              ZZ(normalization)}
        self._I0 = T.cartan_type().classical().index_set()
        from sage.categories.homset import Hom
        Map.__init__(self, Hom(T, ZZ))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['A', 6, 2], 2, 1)
            sage: Kp = crystals.KirillovReshetikhin(['A', 6, 2], 1, 1)
            sage: H = K.local_energy_function(Kp); H
            Local energy function of
             Kirillov-Reshetikhin crystal of type ['BC', 3, 2] with (r,s)=(2,1)
            tensor
             Kirillov-Reshetikhin crystal of type ['BC', 3, 2] with (r,s)=(1,1)
        """
        return "Local energy function of {} tensor {}".format(self._B, self._Bp)

    def _call_(self, x):
        """
        Return the local energy of ``x``.

        EXAMPLES::

            sage: K = crystals.KirillovReshetikhin(['B',4,1], 1,2)
            sage: K2 = crystals.KirillovReshetikhin(['B',4,1], 2,1)
            sage: H = K.local_energy_function(K2)
            sage: T = tensor([K, K2])
            sage: hw = [x for x in T if x.is_highest_weight([1,2])]
            sage: H(hw[0])
            1
        """
        # Setup variables
        visited = {x: 0}
        check0 = [x]

        # Helper function
        def to_classical_hw(cur):
            for i in self._I0:
                b = cur.e(i)
                if b is not None and b not in visited:
                    visited[b] = visited[cur] # No change
                    return b
            return None # is classically HW or all have been visited

        cur = x
        # Get the affine node (it might not be 0 if the type
        #   has been relabeled)
        i0 = x.parent().cartan_type().special_node()
        while cur not in self._known_values:
            # We first go towards the classically highest weight since
            #   the maximal vector is classically highest weight
            b = to_classical_hw(cur)

            # If classically HW, then try 0 arrows
            while b is None:
                b = check0.pop()
                c = b.e(i0)
                # If there is no 0 arrow or we have already seen c, move along
                if c is None or c in visited:
                    b = None
                    continue

                bp = self._R_matrix(b)
                cp = bp.e(i0)
                if b[1] == c[1] and bp[1] == cp[1]: # LL case
                    visited[c] = visited[b] + 1
                elif b[0] == c[0] and bp[0] == cp[0]: # RR case
                    visited[c] = visited[b] - 1
                else:
                    visited[c] = visited[b] # Otherwise no change
                b = c

            cur = b
            check0.append(b)

        baseline = self._known_values[cur] - visited[cur]
        for y in visited:
            self._known_values[y] = baseline + visited[y]

        return self._known_values[x]

