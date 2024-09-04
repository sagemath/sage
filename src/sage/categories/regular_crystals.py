# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.combinat sage.graphs
r"""
Regular Crystals
"""
#*****************************************************************************
#  Copyright (C) 2013    Anne Schilling <anne at math.ucdavis.edu>
#                        Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.tensor import TensorProductsCategory


class RegularCrystals(Category_singleton):
    r"""
    The category of regular crystals.

    A crystal is called *regular* if every vertex `b` satisfies

    .. MATH::

        \varepsilon_i(b) = \max\{ k \mid e_i^k(b) \neq 0 \} \quad \text{and}
        \quad \varphi_i(b) = \max\{ k \mid f_i^k(b) \neq 0 \}.

    .. NOTE::

        Regular crystals are sometimes referred to as *normal*. When only one
        of the conditions (on either `\varphi_i` or `\varepsilon_i`) holds,
        these crystals are sometimes called *seminormal* or *semiregular*.

    EXAMPLES::

        sage: C = RegularCrystals()
        sage: C
        Category of regular crystals
        sage: C.super_categories()
        [Category of crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = RegularCrystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: RegularCrystals().super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    def example(self, n=3):
        """
        Return an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = RegularCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
        return Crystals().example(n)

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of regular crystals defines no new
        structure: it only relates `\varepsilon_a` and `\varphi_a` to
        `e_a` and `f_a` respectively.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`CategoryWithAxiom`?

        EXAMPLES::

            sage: RegularCrystals().additional_structure()
        """
        return None

    class MorphismMethods:
        def is_isomorphism(self):
            """
            Check if ``self`` is a crystal isomorphism, which is true
            if and only if this is a strict embedding with the same number
            of connected components.

            EXAMPLES::

                sage: A21 = RootSystem(['A',2,1])
                sage: La = A21.weight_space(extended=True).fundamental_weights()
                sage: B = crystals.LSPaths(La[0])
                sage: La = A21.weight_lattice(extended=True).fundamental_weights()
                sage: C = crystals.GeneralizedYoungWalls(2, La[0])
                sage: H = Hom(B, C)
                sage: from sage.categories.highest_weight_crystals import HighestWeightCrystalMorphism
                sage: class Psi(HighestWeightCrystalMorphism):
                ....:     def is_strict(self):
                ....:         return True
                sage: psi = Psi(H, C.module_generators); psi
                ['A', 2, 1] Crystal morphism:
                  From: The crystal of LS paths of type ['A', 2, 1] and weight Lambda[0]
                  To:   Highest weight crystal of generalized Young walls
                        of Cartan type ['A', 2, 1] and highest weight Lambda[0]
                  Defn: (Lambda[0],) |--> []
                sage: psi.is_isomorphism()
                True
            """
            return (self.is_strict()
                    and self.domain().number_of_connected_components() ==
                        self.codomain().number_of_connected_components())

    class ParentMethods:

        # TODO: this could be a method in Crystals.Algebras.ElementMethods, so that
        # one could do:
        #
        #   sage: C = crystals.Tableaux(['A',2], shape=[2,1])
        #   sage: M = C.algebra(QQ)
        #   sage: m = M.an_element()
        #   sage: m.demazure_operator([1,4,2])
        def demazure_operator(self, element, reduced_word):
            r"""
            Return the application of Demazure operators `D_i` for `i` from
            ``reduced_word`` on ``element``.

            INPUT:

            - ``element`` -- an element of a free module indexed by the
              underlying crystal
            - ``reduced_word`` -- a reduced word of the Weyl group of the
              same type as the underlying crystal

            OUTPUT: an element of the free module indexed by the underlying crystal

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = CombinatorialFreeModule(QQ, T)                                # needs sage.modules
                sage: t = T.highest_weight_vector()
                sage: b = 2*C(t)                                                        # needs sage.modules
                sage: T.demazure_operator(b,[1,2,1])                                    # needs sage.modules
                2*B[[[1, 1], [2]]] + 2*B[[[1, 2], [2]]] + 2*B[[[1, 3], [2]]]
                + 2*B[[[1, 1], [3]]] + 2*B[[[1, 2], [3]]] + 2*B[[[1, 3], [3]]]
                + 2*B[[[2, 2], [3]]] + 2*B[[[2, 3], [3]]]

            The Demazure operator is idempotent::

                sage: # needs sage.modules
                sage: T = crystals.Tableaux("A1", shape=[4])
                sage: C = CombinatorialFreeModule(QQ, T)
                sage: b = C(T.module_generators[0]); b
                B[[[1, 1, 1, 1]]]
                sage: e = T.demazure_operator(b,[1]); e
                B[[[1, 1, 1, 1]]] + B[[[1, 1, 1, 2]]] + B[[[1, 1, 2, 2]]]
                + B[[[1, 2, 2, 2]]] + B[[[2, 2, 2, 2]]]
                sage: e == T.demazure_operator(e,[1])
                True
                sage: all(T.demazure_operator(T.demazure_operator(C(t),[1]),[1])
                ....:      == T.demazure_operator(C(t),[1]) for t in T)
                True
            """
            M = element.parent()
            for i in reversed(reduced_word):
                element = M.linear_combination((c.demazure_operator_simple(i), coeff)
                                               for c, coeff in element)
            return element

        def demazure_subcrystal(self, element, reduced_word, only_support=True):
            r"""
            Return the subcrystal corresponding to the application of
            Demazure operators `D_i` for `i` from ``reduced_word`` on
            ``element``.

            INPUT:

            - ``element`` -- an element of a free module indexed by the
              underlying crystal
            - ``reduced_word`` -- a reduced word of the Weyl group of the
              same type as the underlying crystal
            - ``only_support`` -- boolean (default: ``True``); only include
              arrows corresponding to the support of ``reduced_word``

            OUTPUT: the Demazure subcrystal

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T.highest_weight_vector()
                sage: S = T.demazure_subcrystal(t, [1,2])
                sage: list(S)
                [[[1, 1], [2]], [[1, 2], [2]], [[1, 1], [3]],
                 [[1, 2], [3]], [[2, 2], [3]]]
                sage: S = T.demazure_subcrystal(t, [2,1])
                sage: list(S)
                [[[1, 1], [2]], [[1, 2], [2]], [[1, 1], [3]],
                 [[1, 3], [2]], [[1, 3], [3]]]

            We construct an example where we don't only want the arrows
            indicated by the support of the reduced word::

                sage: K = crystals.KirillovReshetikhin(['A',1,1], 1, 2)
                sage: mg = K.module_generator()
                sage: S = K.demazure_subcrystal(mg, [1])
                sage: S.digraph().edges(sort=True)
                [([[1, 1]], [[1, 2]], 1), ([[1, 2]], [[2, 2]], 1)]
                sage: S = K.demazure_subcrystal(mg, [1], only_support=False)
                sage: S.digraph().edges(sort=True)
                [([[1, 1]], [[1, 2]], 1),
                 ([[1, 2]], [[1, 1]], 0),
                 ([[1, 2]], [[2, 2]], 1),
                 ([[2, 2]], [[1, 2]], 0)]
            """
            from sage.combinat.free_module import CombinatorialFreeModule
            from sage.rings.rational_field import QQ
            C = CombinatorialFreeModule(QQ, self)
            D = self.demazure_operator(C(element), reduced_word)
            if only_support:
                index_set = tuple(frozenset(reduced_word))
            else:
                index_set = self.cartan_type().index_set()
            return self.subcrystal(contained=D.support(), generators=[element],
                                   index_set=index_set)

        def _test_stembridge_local_axioms(self, index_set=None, verbose=False, complete=False, **options):
            r"""
            This implements tests for the Stembridge local characterization
            on the finite crystal ``self``.

            The current implementation only uses the rules for simply-laced
            types.  Crystals of other types should still pass the test, but
            expansion of this test to non-simply laced type would be desirable.

            One can specify an index set smaller than the full index set of
            the crystal, using the option ``index_set``.

            Running with ``verbose=True`` will print each node for which a
            local axiom test applies.

            Running with ``complete=True`` will continue to run the test past
            the first failure of the local axioms.  This is probably only
            useful in conjunction with the verbose option, to see all places
            where the local axioms fail.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',3], shape=[2,1])
                sage: T._test_stembridge_local_axioms()
                True
                sage: T._test_stembridge_local_axioms(verbose=True)
                True
                sage: T._test_stembridge_local_axioms(index_set=[1,3])
                True

                sage: B = Crystals().example(choice='naive')
                sage: B._test_stembridge_local_axioms()
                Traceback (most recent call last):
                ...
                AssertionError: None
            """
            tester = self._tester(**options)
            goodness = True
            i = 0
            for x in self:
                goodness = x._test_stembridge_local_axioms(index_set, verbose)
                if not goodness and not complete:
                    tester.fail()
                i += 1
                if i > tester._max_runs:
                    return
            tester.assertTrue(goodness)
            return goodness

        def dual_equivalence_graph(self, X=None, index_set=None, directed=True):
            r"""
            Return the dual equivalence graph indexed by ``index_set``
            on the subset ``X`` of ``self``.

            Let `b \in B` be an element of weight `0`, so `\varepsilon_j(b)
            = \varphi_j(b)` for all `j \in I`, where `I` is the indexing
            set. We say `b'` is an `i`-elementary dual equivalence
            transformation of `b` (where `i \in I`) if

            * `\varepsilon_i(b) = 1` and `\varepsilon_{i-1}(b) = 0`, and
            * `b' = f_{i-1} f_i e_{i-1} e_i b`.

            We can do the inverse procedure by interchanging `i` and `i-1`
            above.

            .. NOTE::

                If the index set is not an ordered interval, we let
                `i - 1` mean the index appearing before `i` in `I`.

            This definition comes from [As2008]_ Section 4 (where our
            `\varphi_j(b)` and `\varepsilon_j(b)` are denoted by
            `\epsilon(b, j)` and `-\delta(b, j)`, respectively).

            The dual equivalence graph of `B` is defined to be the
            colored graph whose vertices are the elements of `B` of
            weight `0`, and whose edges of color `i` (for `i \in I`)
            connect pairs `\{ b, b' \}` such that `b'` is an
            `i`-elementary dual equivalence transformation of `b`.

            .. NOTE::

                This dual equivalence graph is a generalization of
                `\mathcal{G}\left(\mathcal{X}\right)` in [As2008]_
                Section 4 except we do not require
                `\varepsilon_i(b) = 0, 1` for all `i`.

            This definition can be generalized by choosing a subset `X`
            of the set of all vertices of `B` of weight `0`, and
            restricting the dual equivalence graph to the vertex set
            `X`.

            INPUT:

            - ``X`` -- (optional) the vertex set `X` (default:
              the whole set of vertices of ``self`` of weight `0`)
            - ``index_set`` -- (optional) the index set `I`
              (default: the whole index set of ``self``); this has
              to be a subset of the index set of ``self`` (as a list
              or tuple)
            - ``directed`` -- boolean (default: ``True``); whether to have the
              dual equivalence graph be directed, where the head of
              an edge `b - b'` is `b` and the tail is
              `b' = f_{i-1} f_i e_{i-1} e_i b`)

            .. SEEALSO::

                :meth:`sage.combinat.partition.Partition.dual_equivalence_graph`

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',3], shape=[2,2])
                sage: G = T.dual_equivalence_graph()
                sage: G.edges(sort=True)
                [([[1, 3], [2, 4]], [[1, 2], [3, 4]], 2),
                 ([[1, 2], [3, 4]], [[1, 3], [2, 4]], 3)]
                sage: T = crystals.Tableaux(['A',4], shape=[3,2])
                sage: G = T.dual_equivalence_graph()
                sage: G.edges(sort=True)
                [([[1, 3, 5], [2, 4]], [[1, 3, 4], [2, 5]], 4),
                 ([[1, 3, 5], [2, 4]], [[1, 2, 5], [3, 4]], 2),
                 ([[1, 3, 4], [2, 5]], [[1, 2, 4], [3, 5]], 2),
                 ([[1, 2, 5], [3, 4]], [[1, 3, 5], [2, 4]], 3),
                 ([[1, 2, 4], [3, 5]], [[1, 2, 3], [4, 5]], 3),
                 ([[1, 2, 3], [4, 5]], [[1, 2, 4], [3, 5]], 4)]

                sage: T = crystals.Tableaux(['A',4], shape=[3,1])
                sage: G = T.dual_equivalence_graph(index_set=[1,2,3])
                sage: G.vertices(sort=True)
                [[[1, 3, 4], [2]], [[1, 2, 4], [3]], [[1, 2, 3], [4]]]
                sage: G.edges(sort=True)
                [([[1, 3, 4], [2]], [[1, 2, 4], [3]], 2),
                 ([[1, 2, 4], [3]], [[1, 2, 3], [4]], 3)]

            TESTS::

                sage: T = crystals.Tableaux(['A',4], shape=[3,1])
                sage: G = T.dual_equivalence_graph(index_set=[2,3])
                sage: G.edges(sort=True)
                [([[1, 2, 4], [3]], [[1, 2, 3], [4]], 3),
                 ([[2, 4, 5], [3]], [[2, 3, 5], [4]], 3)]
                sage: G.vertices(sort=True)
                [[[1, 3, 4], [2]],
                 [[1, 2, 4], [3]],
                 [[2, 4, 5], [3]],
                 [[1, 2, 3], [4]],
                 [[2, 3, 5], [4]],
                 [[1, 1, 1], [5]],
                 [[1, 1, 5], [5]],
                 [[1, 5, 5], [5]],
                 [[2, 3, 4], [5]]]
            """
            if index_set is None:
                index_set = self.index_set()

            def wt_zero(x):
                for i in index_set:
                    if x.epsilon(i) != x.phi(i):
                        return False
                return True

            if X is None:
                X = [x for x in self if wt_zero(x)]
                checker = lambda x: True
            elif any(not wt_zero(x) for x in X):
                raise ValueError("the elements are not all weight 0")
            else:
                checker = lambda x: x in X

            edges = []
            for x in X:
                for k, i in enumerate(index_set[1:]):
                    im = index_set[k]
                    if x.epsilon(i) == 1 and x.epsilon(im) == 0:
                        y = x.e(i).e(im).f(i).f(im)
                        if checker(y):
                            edges.append([x, y, i])
            from sage.graphs.digraph import DiGraph
            G = DiGraph([X, edges], format='vertices_and_edges', immutable=True)
            from sage.graphs.dot2tex_utils import have_dot2tex
            if have_dot2tex():
                G.set_latex_options(format='dot2tex', edge_labels=True,
                                    color_by_label=self.cartan_type()._index_set_coloring)
            return G

    class ElementMethods:

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).epsilon(1)
                0
                sage: C(2).epsilon(1)
                1
            """
            assert i in self.index_set()
            x = self.e(i)
            eps = 0
            while x is not None:
                x = x.e(i)
                eps = eps + 1
            return eps

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).phi(1)
                1
                sage: C(2).phi(1)
                0
            """
            assert i in self.index_set()
            x = self.f(i)
            phi = 0
            while x is not None:
                x = x.f(i)
                phi += 1
            return phi

        def weight(self):
            """
            Return the weight of this crystal element.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).weight()
                (1, 0, 0, 0, 0, 0)
            """
            return self.Phi() - self.Epsilon()

        def demazure_operator_simple(self, i, ring=None):
            r"""
            Return the Demazure operator `D_i` applied to ``self``.

            INPUT:

            - ``i`` -- an element of the index set of the underlying crystal
            - ``ring`` -- (default: ``QQ``) a ring

            OUTPUT:

            An element of the ``ring``-free module indexed by the underlying
            crystal.

            Let `r = \langle \mathrm{wt}(b), \alpha^{\vee}_i \rangle`, then
            `D_i(b)` is defined as follows:

            - If `r \geq 0`, this returns the sum of the elements obtained
              from ``self`` by application of `f_i^k` for `0 \leq k \leq r`.
            - If `r < 0`, this returns the opposite of the sum of the
              elements obtained by application of `e_i^k` for `0 < k < -r`.

            REFERENCES:

            - [Li1995]_

            - [Ka1993]_

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,2],[2]])
                sage: t.demazure_operator_simple(2)
                B[[[1, 2], [2]]] + B[[[1, 3], [2]]] + B[[[1, 3], [3]]]
                sage: t.demazure_operator_simple(2).parent()
                Algebra of The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
                        over Integer Ring

                sage: t.demazure_operator_simple(1)
                0

                sage: K = crystals.KirillovReshetikhin(['A',2,1],2,1)
                sage: t = K(rows=[[3],[2]])
                sage: t.demazure_operator_simple(0)
                B[[[1, 2]]] + B[[[2, 3]]]

            TESTS::

                sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                sage: x = K.an_element(); x
                [[1]]
                sage: x.demazure_operator_simple(0)
                0
                sage: x.demazure_operator_simple(0, ring = QQ).parent()
                Algebra of Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1)
                        over Rational Field
            """
            from sage.rings.integer_ring import ZZ
            if ring is None:
                ring = ZZ
            C = self.parent().algebra(ring)
            r = self.phi(i) - self.epsilon(i)
            if r >= 0:
                l = [self]
                element = self
                for k in range(r):
                    element = element.f(i)
                    l.append(element)
                return C.sum_of_monomials(l)
            else:
                l = []
                element = self
                for k in range(-r-1):
                    element = element.e(i)
                    l.append(element)
                return - C.sum_of_monomials(l)

        def stembridgeDelta_depth(self,i,j):
            r"""
            Return the difference in the `j`-depth of ``self`` and `e_i`
            of ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-depth of a crystal node `x` is `-\varepsilon_i(x)`.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,2],[2]])
                sage: t.stembridgeDelta_depth(1,2)
                0
                sage: s = T(rows=[[2,3],[3]])
                sage: s.stembridgeDelta_depth(1,2)
                -1
            """
            if self.e(i) is None:
                return 0
            return -self.e(i).epsilon(j) + self.epsilon(j)

        def stembridgeDelta_rise(self,i,j):
            r"""
            Return the difference in the `j`-rise of ``self`` and `e_i` of
            ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-rise of a crystal node `x` is `\varphi_i(x)`.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,2],[2]])
                sage: t.stembridgeDelta_rise(1,2)
                -1
                sage: s = T(rows=[[2,3],[3]])
                sage: s.stembridgeDelta_rise(1,2)
                0
            """
            if self.e(i) is None:
                return 0
            return self.e(i).phi(j) - self.phi(j)

        def stembridgeDel_depth(self,i,j):
            r"""
            Return the difference in the `j`-depth of ``self`` and `f_i` of
            ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-depth of a crystal node `x` is `\varepsilon_i(x)`.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,1],[2]])
                sage: t.stembridgeDel_depth(1,2)
                0
                sage: s = T(rows=[[1,3],[3]])
                sage: s.stembridgeDel_depth(1,2)
                -1
            """
            if self.f(i) is None:
                return 0
            return -self.epsilon(j) + self.f(i).epsilon(j)

        def stembridgeDel_rise(self,i,j):
            r"""
            Return the difference in the `j`-rise of ``self`` and `f_i` of
            ``self``, where `i` and `j` are in the index set of the
            underlying crystal. This function is useful for checking the
            Stembridge local axioms for crystal bases.

            The `i`-rise of a crystal node `x` is `\varphi_i(x)`.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,1],[2]])
                sage: t.stembridgeDel_rise(1,2)
                -1
                sage: s = T(rows=[[1,3],[3]])
                sage: s.stembridgeDel_rise(1,2)
                0
            """
            if self.f(i) is None:
                return 0
            return self.phi(j)-self.f(i).phi(j)

        def stembridgeTriple(self,i,j):
            r"""
            Let `A` be the Cartan matrix of the crystal, `x` a crystal element,
            and let `i` and `j` be in the index set of the crystal.
            Further, set
            ``b=stembridgeDelta_depth(x,i,j)``, and
            ``c=stembridgeDelta_rise(x,i,j))``.
            If ``x.e(i)`` is non-empty, this function returns the triple
            `( A_{ij}, b, c )`; otherwise it returns ``None``.
            By the Stembridge local characterization of crystal bases,
            one should have `A_{ij}=b+c`.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,1],[2]])
                sage: t.stembridgeTriple(1,2)
                sage: s = T(rows=[[1,2],[2]])
                sage: s.stembridgeTriple(1,2)
                (-1, 0, -1)

                sage: T = crystals.Tableaux(['B',2], shape=[2,1])
                sage: t = T(rows=[[1,2],[2]])
                sage: t.stembridgeTriple(1,2)
                (-2, 0, -2)
                sage: s = T(rows=[[-1,-1],[0]])
                sage: s.stembridgeTriple(1,2)
                (-2, -2, 0)
                sage: u = T(rows=[[0,2],[1]])
                sage: u.stembridgeTriple(1,2)
                (-2, -1, -1)
            """
            if self.e(i) is None:
                return None
            b = self.stembridgeDelta_depth(i,j)
            c = self.stembridgeDelta_rise(i,j)
            dd = self.cartan_type().dynkin_diagram()
            a = dd[j,i]
            return (a, b, c)

        def _test_stembridge_local_axioms(self, index_set=None, verbose=False, **options):
            r"""
            This implements tests for the Stembridge local characterization
            on the element of a crystal ``self``.

            The current implementation only uses the axioms for simply-laced
            types.  Crystals of other types should still pass the test, but
            in non-simply-laced types, passing is not a guarantee that the
            crystal arises from a representation.

            One can specify an index set smaller than the full index set of
            the crystal, using the option ``index_set``.

            Running with ``verbose=True`` will print warnings when a test fails.

            REFERENCES:

            - [Ste2003]_

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: t = T(rows=[[1,1],[2]])
                sage: t._test_stembridge_local_axioms()
                True
                sage: t._test_stembridge_local_axioms(index_set=[1,3])
                True
                sage: t._test_stembridge_local_axioms(verbose=True)
                True
            """
            tester = self._tester(**options)
            goodness = True
            if index_set is None:
                index_set = self.index_set()

            from sage.combinat.subset import Subsets

            for (i,j) in Subsets(index_set, 2):
                if self.e(i) is not None and self.e(j) is not None:
                    triple = self.stembridgeTriple(i,j)
                    #Test axioms P3 and P4.
                    if not triple[0] == triple[1]+triple[2] or triple[1] > 0 or triple[2] > 0:
                        if verbose:
                            print('Warning: Failed axiom P3 or P4 at vector ', self, 'i,j=', i, j, 'Stembridge triple:', self.stembridgeTriple(i, j))
                            goodness = False
                        else:
                            tester.fail()
                    if self.stembridgeDelta_depth(i,j) == 0:
                        #check E_i E_j(x)= E_j E_i(x)
                        if self.e(i).e(j) != self.e(j).e(i) or self.e(i).e(j).stembridgeDel_rise(j, i) != 0:
                            if verbose:
                                print('Warning: Failed axiom P5 at: vector ', self, 'i,j=', i, j, 'Stembridge triple:', self.stembridgeTriple(i, j))
                                goodness = False
                            else:
                                tester.fail()
                    if self.stembridgeDelta_depth(i,j) == -1 and self.stembridgeDelta_depth(j,i) == -1:
                        #check E_i E_j^2 E_i (x)= E_j E_i^2 E_j (x)
                        y1 = self.e(j).e(i).e(i).e(j)
                        y2 = self.e(j).e(i).e(i).e(j)
                        a = y1.stembridgeDel_rise(j, i)
                        b = y2.stembridgeDel_rise(i, j)
                        if y1 != y2 or a != -1 or b != -1:
                            if verbose:
                                print('Warning: Failed axiom P6 at: vector ', self, 'i,j=', i, j, 'Stembridge triple:', self.stembridgeTriple(i, j))
                                goodness = False
                            else:
                                tester.fail()
            tester.assertTrue(goodness)
            return goodness

        def dual_equivalence_class(self, index_set=None):
            r"""
            Return the dual equivalence class indexed by ``index_set``
            of ``self``.

            The dual equivalence class of an element `b \in B`
            is the set of all elements of `B` reachable from
            `b` via sequences of `i`-elementary dual equivalence
            relations (i.e., `i`-elementary dual equivalence
            transformations and their inverses) for `i` in the index
            set of `B`.

            For this to be well-defined, the element `b` has to be
            of weight `0` with respect to `I`; that is, we need to have
            `\varepsilon_j(b) = \varphi_j(b)` for all `j \in I`.

            See [As2008]_. See also :meth:`dual_equivalence_graph` for
            a definition of `i`-elementary dual equivalence
            transformations.

            INPUT:

            - ``index_set`` -- (optional) the index set `I`
              (default: the whole index set of the crystal); this has
              to be a subset of the index set of the crystal (as a list
              or tuple)

            OUTPUT:

            The dual equivalence class of ``self`` indexed by the
            subset ``index_set``. This class is returned as an
            undirected edge-colored multigraph. The color of an edge
            is the index `i` of the dual equivalence relation it
            encodes.

            .. SEEALSO::

                - :meth:`~sage.categories.regular_crystals.RegularCrystals.ParentMethods.dual_equivalence_graph`
                - :meth:`sage.combinat.partition.Partition.dual_equivalence_graph`

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',3], shape=[2,2])
                sage: G = T(2,1,4,3).dual_equivalence_class()
                sage: G.edges(sort=True)
                [([[1, 3], [2, 4]], [[1, 2], [3, 4]], 2),
                 ([[1, 3], [2, 4]], [[1, 2], [3, 4]], 3)]
                sage: T = crystals.Tableaux(['A',4], shape=[3,2])
                sage: G = T(2,1,4,3,5).dual_equivalence_class()
                sage: G.edges(sort=True,sort_vertices=True)
                [([[1, 3, 4], [2, 5]], [[1, 3, 5], [2, 4]], 4),
                 ([[1, 3, 4], [2, 5]], [[1, 2, 4], [3, 5]], 2),
                 ([[1, 2, 5], [3, 4]], [[1, 3, 5], [2, 4]], 2),
                 ([[1, 2, 5], [3, 4]], [[1, 3, 5], [2, 4]], 3),
                 ([[1, 2, 3], [4, 5]], [[1, 2, 4], [3, 5]], 3),
                 ([[1, 2, 3], [4, 5]], [[1, 2, 4], [3, 5]], 4)]
            """
            if index_set is None:
                index_set = self.index_set()

            for i in index_set:
                if self.epsilon(i) != self.phi(i):
                    raise ValueError("the element is not weight 0")

            visited = set()
            todo = {self}
            edges = []
            while todo:
                x = todo.pop()
                visited.add(x)
                for k, i in enumerate(index_set[1:]):
                    im = index_set[k]
                    if x.epsilon(i) == 1 and x.epsilon(im) == 0:
                        y = x.e(i).e(im).f(i).f(im)
                        if [y, x, i] not in edges:
                            edges.append([x, y, i])
                        if y not in visited:
                            todo.add(y)
                    if x.epsilon(i) == 0 and x.epsilon(im) == 1:
                        y = x.e(im).e(i).f(im).f(i)
                        if [y, x, i] not in edges:
                            edges.append([x, y, i])
                        if y not in visited:
                            todo.add(y)
            from sage.graphs.graph import Graph
            G = Graph([visited, edges], format='vertices_and_edges',
                      immutable=True, multiedges=True)
            from sage.graphs.dot2tex_utils import have_dot2tex
            if have_dot2tex():
                G.set_latex_options(format='dot2tex', edge_labels=True,
                                    color_by_label=self.cartan_type()._index_set_coloring)
            return G

    class TensorProducts(TensorProductsCategory):
        """
        The category of regular crystals constructed by tensor
        product of regular crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: RegularCrystals().TensorProducts().extra_super_categories()
                [Category of regular crystals]
            """
            return [self.base_category()]
