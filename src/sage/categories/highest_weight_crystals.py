# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.combinat sage.graphs
r"""
Highest Weight Crystals
"""
# ****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import (Crystals, CrystalHomset,
                                      CrystalMorphismByGenerators)
from sage.categories.tensor import TensorProductsCategory


class HighestWeightCrystals(Category_singleton):
    """
    The category of highest weight crystals.

    A crystal is highest weight if it is acyclic; in particular, every
    connected component has a unique highest weight element, and that
    element generate the component.

    EXAMPLES::

        sage: C = HighestWeightCrystals()
        sage: C
        Category of highest weight crystals
        sage: C.super_categories()
        [Category of crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = HighestWeightCrystals().example()
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

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: HighestWeightCrystals().super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    def example(self):
        """
        Return an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = HighestWeightCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
        return Crystals().example()

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of highest weight crystals defines no
        additional structure: it only guarantees the existence of a
        unique highest weight element in each component.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`CategoryWithAxiom`?

        EXAMPLES::

            sage: HighestWeightCrystals().additional_structure()
        """
        return None

    class ParentMethods:

        @cached_method
        def highest_weight_vectors(self):
            r"""
            Return the highest weight vectors of ``self``.

            This default implementation selects among the module
            generators those that are highest weight, and caches the result.
            A crystal element `b` is highest weight if `e_i(b)=0` for all `i` in the
            index set.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C.highest_weight_vectors()
                (1,)

            ::

                sage: C = crystals.Letters(['A',2])
                sage: T = crystals.TensorProduct(C, C, C, generators=[[C(2),C(1),C(1)],
                ....:                                                 [C(1),C(2),C(1)]])
                sage: T.highest_weight_vectors()
                ([2, 1, 1], [1, 2, 1])
            """
            return tuple(g for g in self.module_generators if g.is_highest_weight())

        def highest_weight_vector(self):
            r"""
            Return the highest weight vector if there is a single one;
            otherwise, raises an error.

            Caveat: this assumes that :meth:`.highest_weight_vectors`
            returns a list or tuple.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C.highest_weight_vector()
                1
            """
            hw = self.highest_weight_vectors()
            if len(hw) == 1:
                return hw[0]
            else:
                raise RuntimeError("The crystal does not have exactly one highest weight vector")

        # TODO: Not every highest weight crystal is a lowest weight crystal
        @cached_method
        def lowest_weight_vectors(self):
            r"""
            Return the lowest weight vectors of ``self``.

            This default implementation selects among all elements of the crystal
            those that are lowest weight, and cache the result.
            A crystal element `b` is lowest weight if `f_i(b)=0` for all `i` in the
            index set.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C.lowest_weight_vectors()
                (6,)

            ::

                sage: C = crystals.Letters(['A',2])
                sage: T = crystals.TensorProduct(C, C, C,generators=[[C(2),C(1),C(1)],
                ....:                                                [C(1),C(2),C(1)]])
                sage: T.lowest_weight_vectors()
                ([3, 2, 3], [3, 3, 2])
            """
            return tuple(g for g in self if g.is_lowest_weight())

        def __iter__(self, index_set=None, max_depth=float("inf")):
            """
            Return the iterator of ``self``.

            INPUT:

            - ``index_set`` -- (default: ``None``) the index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (default: infinity) the maximum depth to build

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2,1],[0,1,0])
                sage: sorted([p for p in C.__iter__(max_depth=3)], key=str)
                [(-Lambda[0] + 2*Lambda[2] - delta,),
                 (-Lambda[0] + Lambda[1] + 1/2*Lambda[2] - delta, Lambda[0] - 1/2*Lambda[2]),
                 (1/2*Lambda[0] + Lambda[1] - Lambda[2] - 1/2*delta, -1/2*Lambda[0] + Lambda[2] - 1/2*delta),
                 (2*Lambda[0] - Lambda[2],),
                 (Lambda[0] - Lambda[1] + Lambda[2],),
                 (Lambda[1],)]
                sage: [p for p in C.__iter__(index_set=[0, 1], max_depth=3)]
                [(Lambda[1],), (Lambda[0] - Lambda[1] + Lambda[2],), (-Lambda[0] + 2*Lambda[2] - delta,)]
            """
            if index_set is None:
                index_set = self.index_set()
            from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
            return RecursivelyEnumeratedSet(self.module_generators,
                           lambda x: [x.f(i) for i in index_set],
                           structure='graded',
                           max_depth=max_depth).breadth_first_search_iterator()

        @cached_method
        def q_dimension(self, q=None, prec=None, use_product=False):
            r"""
            Return the `q`-dimension of ``self``.

            Let `B(\lambda)` denote a highest weight crystal. Recall that
            the degree of the `\mu`-weight space of `B(\lambda)` (under
            the principal gradation) is equal to
            `\langle \rho^{\vee}, \lambda - \mu \rangle` where
            `\langle \rho^{\vee}, \alpha_i \rangle = 1` for all `i \in I`
            (in particular, take `\rho^{\vee} = \sum_{i \in I} h_i`).

            The `q`-dimension of a highest weight crystal `B(\lambda)` is
            defined as

            .. MATH::

                \dim_q B(\lambda) := \sum_{j \geq 0} \dim(B_j) q^j,

            where `B_j` denotes the degree `j` portion of `B(\lambda)`. This
            can be expressed as the product

            .. MATH::

                \dim_q B(\lambda) = \prod_{\alpha^{\vee} \in \Delta_+^{\vee}}
                \left( \frac{1 - q^{\langle \lambda + \rho, \alpha^{\vee}
                \rangle}}{1 - q^{\langle \rho, \alpha^{\vee} \rangle}}
                \right)^{\mathrm{mult}\, \alpha},

            where `\Delta_+^{\vee}` denotes the set of positive coroots.
            Taking the limit as `q \to 1` gives the dimension of `B(\lambda)`.
            For more information, see [Ka1990]_ Section 10.10.

            INPUT:

            - ``q`` -- the (generic) parameter `q`

            - ``prec`` -- (default: ``None``) the precision of the power
              series ring to use if the crystal is not known to be finite
              (i.e. the number of terms returned).
              If ``None``, then the result is returned as a lazy power series.

            - ``use_product`` -- boolean (default: ``False``); if we have a
              finite crystal and ``True``, use the product formula

            EXAMPLES::

                sage: C = crystals.Tableaux(['A',2], shape=[2,1])
                sage: qdim = C.q_dimension(); qdim
                q^4 + 2*q^3 + 2*q^2 + 2*q + 1
                sage: qdim(1)
                8
                sage: len(C) == qdim(1)
                True
                sage: C.q_dimension(use_product=True) == qdim
                True
                sage: C.q_dimension(prec=20)
                q^4 + 2*q^3 + 2*q^2 + 2*q + 1
                sage: C.q_dimension(prec=2)
                2*q + 1

                sage: R.<t> = QQ[]
                sage: C.q_dimension(q=t^2)
                t^8 + 2*t^6 + 2*t^4 + 2*t^2 + 1

                sage: C = crystals.Tableaux(['A',2], shape=[5,2])
                sage: C.q_dimension()
                q^10 + 2*q^9 + 4*q^8 + 5*q^7 + 6*q^6 + 6*q^5
                 + 6*q^4 + 5*q^3 + 4*q^2 + 2*q + 1

                sage: C = crystals.Tableaux(['B',2], shape=[2,1])
                sage: qdim = C.q_dimension(); qdim
                q^10 + 2*q^9 + 3*q^8 + 4*q^7 + 5*q^6 + 5*q^5
                 + 5*q^4 + 4*q^3 + 3*q^2 + 2*q + 1
                sage: qdim == C.q_dimension(use_product=True)
                True

                sage: C = crystals.Tableaux(['D',4], shape=[2,1])
                sage: C.q_dimension()
                q^16 + 2*q^15 + 4*q^14 + 7*q^13 + 10*q^12 + 13*q^11
                 + 16*q^10 + 18*q^9 + 18*q^8 + 18*q^7 + 16*q^6 + 13*q^5
                 + 10*q^4 + 7*q^3 + 4*q^2 + 2*q + 1

            We check with a finite tensor product::

                sage: TP = crystals.TensorProduct(C, C)
                sage: TP.cardinality()
                25600
                sage: qdim = TP.q_dimension(use_product=True); qdim # long time
                q^32 + 2*q^31 + 8*q^30 + 15*q^29 + 34*q^28 + 63*q^27 + 110*q^26
                 + 175*q^25 + 276*q^24 + 389*q^23 + 550*q^22 + 725*q^21
                 + 930*q^20 + 1131*q^19 + 1362*q^18 + 1548*q^17 + 1736*q^16
                 + 1858*q^15 + 1947*q^14 + 1944*q^13 + 1918*q^12 + 1777*q^11
                 + 1628*q^10 + 1407*q^9 + 1186*q^8 + 928*q^7 + 720*q^6
                 + 498*q^5 + 342*q^4 + 201*q^3 + 117*q^2 + 48*q + 26
                sage: qdim(1) # long time
                25600
                sage: TP.q_dimension() == qdim # long time
                True

            The `q`-dimensions of infinite crystals are returned
            as formal power series::

                sage: C = crystals.LSPaths(['A',2,1], [1,0,0])
                sage: C.q_dimension(prec=5)
                1 + q + 2*q^2 + 2*q^3 + 4*q^4 + O(q^5)
                sage: C.q_dimension(prec=10)
                1 + q + 2*q^2 + 2*q^3 + 4*q^4 + 5*q^5 + 7*q^6
                 + 9*q^7 + 13*q^8 + 16*q^9 + O(q^10)
                sage: qdim = C.q_dimension(); qdim
                1 + q + 2*q^2 + 2*q^3 + 4*q^4 + 5*q^5 + 7*q^6 + O(q^7)
                sage: qdim[:16]
                [1, 1, 2, 2, 4, 5, 7, 9, 13, 16, 22, 27, 36, 44, 57, 70]
            """
            from sage.rings.integer_ring import ZZ
            WLR = self.weight_lattice_realization()
            I = self.index_set()
            mg = self.highest_weight_vectors()
            max_deg = float('inf') if prec is None else prec - 1

            def iter_by_deg(gens):
                next = set(gens)
                deg = -1
                while next and deg < max_deg:
                    deg += 1
                    yield len(next)
                    todo = next
                    next = set()
                    while todo:
                        x = todo.pop()
                        for i in I:
                            y = x.f(i)
                            if y is not None:
                                next.add(y)
                # def iter_by_deg

            from sage.categories.finite_crystals import FiniteCrystals
            if self in FiniteCrystals():
                if q is None:
                    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                    q = PolynomialRing(ZZ, 'q').gen(0)

                if use_product:
                    # Since we are in the classical case, all roots occur with multiplicity 1
                    pos_coroots = [x.associated_coroot() for x in WLR.positive_roots()]
                    rho = WLR.rho()
                    P = q.parent()
                    ret = P.zero()
                    for v in self.highest_weight_vectors():
                        hw = v.weight()
                        ret += P.prod((1 - q**(rho+hw).scalar(ac)) / (1 - q**rho.scalar(ac))
                                      for ac in pos_coroots)
                    # We do a cast since the result would otherwise live in the fraction field
                    return P(ret)

            elif prec is None:
                # If we're here, we may not be a finite crystal.
                # In fact, we're probably infinite.
                from sage.rings.lazy_series_ring import LazyPowerSeriesRing
                if q is None:
                    P = LazyPowerSeriesRing(ZZ, names='q')
                else:
                    P = q.parent()
                if not isinstance(P, LazyPowerSeriesRing):
                    raise TypeError("the parent of q must be a lazy power series ring")
                ret = P(iter_by_deg(mg))
                return ret

            from sage.rings.power_series_ring import PowerSeriesRing, PowerSeriesRing_generic
            if q is None:
                q = PowerSeriesRing(ZZ, 'q', default_prec=prec).gen(0)
            P = q.parent()
            ret = P.sum(c * q**deg for deg, c in enumerate(iter_by_deg(mg)))
            if ret.degree() == max_deg and isinstance(P, PowerSeriesRing_generic):
                ret = P(ret, prec)
            return ret

        # TODO: This is not correct if a factor has multiple heads (i.e., we
        #   should have a category for uniqueness of highest/lowest weights)
        connected_components_generators = highest_weight_vectors

        def _Hom_(self, Y, category=None, **options):
            r"""
            Return the homset from ``self`` to ``Y`` in the
            category ``category``.

            INPUT:

            - ``Y`` -- a crystal
            - ``category`` -- a subcategory of :class:`HighestWeightCrysals`()
              or ``None``

            The sole purpose of this method is to construct the homset as a
            :class:`~sage.categories.highest_weight_crystals.HighestWeightCrystalHomset`.
            If ``category`` is specified and is not a subcategory of
            :class:`HighestWeightCrystals`, a :exc:`TypeError` is raised
            instead

            This method is not meant to be called directly. Please use
            :func:`sage.categories.homset.Hom` instead.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: H = B._Hom_(B)
                sage: H
                Set of Crystal Morphisms from The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
                 to The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
                sage: type(H)
                <class 'sage.categories.highest_weight_crystals.HighestWeightCrystalHomset_with_category'>

            TESTS:

            Check that we fallback first to trying a crystal homset
            (:issue:`19458`)::

                sage: Binf = crystals.infinity.Tableaux(['A',2])
                sage: Bi = crystals.elementary.Elementary(Binf.cartan_type(), 1)
                sage: tens = Bi.tensor(Binf)
                sage: Hom(Binf, tens)
                Set of Crystal Morphisms from ...
            """
            if category is None:
                category = self.category()
            elif not category.is_subcategory(Crystals()):
                raise TypeError("{} is not a subcategory of Crystals()".format(category))
            if Y not in Crystals():
                raise TypeError("{} is not a crystal".format(Y))
            return HighestWeightCrystalHomset(self, Y, category=category, **options)

        def digraph(self, subset=None, index_set=None, depth=None):
            """
            Return the DiGraph associated to ``self``.

            INPUT:

            - ``subset`` -- (optional) a subset of vertices for
              which the digraph should be constructed

            - ``index_set`` -- (optional) the index set to draw arrows

            - ``depth`` -- the depth to draw; optional only for finite crystals

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2], shape=[2,1])
                sage: T.digraph()
                Digraph on 8 vertices
                sage: S = T.subcrystal(max_depth=2)
                sage: len(S)
                5
                sage: G = T.digraph(subset=list(S))
                sage: G.is_isomorphic(T.digraph(depth=2), edge_labels=True)
                True

            TESTS:

            The following example demonstrates the speed improvement.
            The speedup in non-affine types is small however::

                sage: depth = 5
                sage: C = crystals.AlcovePaths(['A',2,1], [1,1,0])
                sage: general_digraph = Crystals().parent_class.digraph
                sage: S = C.subcrystal(max_depth=depth, direction='lower')
                sage: %timeit C.digraph(depth=depth) # not tested
                10 loops, best of 3: 48.9 ms per loop
                sage: %timeit general_digraph(C, subset=S) # not tested
                10 loops, best of 3: 96.5 ms per loop
                sage: G1 = C.digraph(depth=depth)
                sage: G2 = general_digraph(C, subset=S)
                sage: G1.is_isomorphic(G2, edge_labels=True)
                True
            """
            if subset is not None:
                return Crystals().parent_class.digraph(self, subset, index_set)

            if self not in Crystals().Finite() and depth is None:
                raise NotImplementedError("crystals not known to be finite must"
                                          " specify either the subset or depth")

            from sage.graphs.digraph import DiGraph
            if index_set is None:
                index_set = self.index_set()

            rank = 0
            d = {g: {} for g in self.module_generators}
            visited = set(d.keys())

            while depth is None or rank < depth:
                recently_visited = set()
                for x in visited:
                    d.setdefault(x, {}) # does nothing if there's a default
                    for i in index_set:
                        xfi = x.f(i)
                        if xfi is not None:
                            d[x][xfi] = i
                            recently_visited.add(xfi)
                if not recently_visited: # No new nodes, nothing more to do
                    break
                rank += 1
                visited = recently_visited

            G = DiGraph(d)
            from sage.graphs.dot2tex_utils import have_dot2tex
            if have_dot2tex():
                G.set_latex_options(format='dot2tex',
                                    edge_labels=True,
                                    color_by_label=self.cartan_type()._index_set_coloring)
            return G

    class ElementMethods:
        def string_parameters(self, word=None):
            r"""
            Return the string parameters of ``self`` corresponding to the
            reduced word ``word``.

            Given a reduced expression `w = s_{i_1} \cdots s_{i_k}`,
            the string parameters of `b \in B` corresponding to `w`
            are `(a_1, \ldots, a_k)` such that

            .. MATH::

                \begin{aligned}
                e_{i_m}^{a_m} \cdots e_{i_1}^{a_1} b & \neq 0 \\
                e_{i_m}^{a_m+1} \cdots e_{i_1}^{a_1} b & = 0
                \end{aligned}

            for all `1 \leq m \leq k`.

            For connected components isomorphic to `B(\lambda)` or
            `B(\infty)`, if `w = w_0` is the longest element of the
            Weyl group, then the path determined by the string
            parametrization terminates at the highest weight vector.

            INPUT:

            - ``word`` -- a word in the alphabet of the index set; if not
              specified and we are in finite type, then this will be some
              reduced expression for the long element determined by the
              Weyl group

            EXAMPLES::

                sage: B = crystals.infinity.NakajimaMonomials(['A',3])
                sage: mg = B.highest_weight_vector()
                sage: w0 = [1,2,1,3,2,1]
                sage: mg.string_parameters(w0)
                [0, 0, 0, 0, 0, 0]
                sage: mg.f_string([1]).string_parameters(w0)
                [1, 0, 0, 0, 0, 0]
                sage: mg.f_string([1,1,1]).string_parameters(w0)
                [3, 0, 0, 0, 0, 0]
                sage: mg.f_string([1,1,1,2,2]).string_parameters(w0)
                [1, 2, 2, 0, 0, 0]
                sage: mg.f_string([1,1,1,2,2]) == mg.f_string([1,1,2,2,1])
                True
                sage: x = mg.f_string([1,1,1,2,2,1,3,3,2,1,1,1])
                sage: x.string_parameters(w0)
                [4, 1, 1, 2, 2, 2]
                sage: x.string_parameters([3,2,1,3,2,3])
                [2, 3, 7, 0, 0, 0]
                sage: x == mg.f_string([1]*7 + [2]*3 + [3]*2)
                True

            ::

                sage: B = crystals.infinity.Tableaux("A5")
                sage: b = B(rows=[[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,6,6,6,6,6,6],
                ....:             [2,2,2,2,2,2,2,2,2,4,5,5,5,6],
                ....:             [3,3,3,3,3,3,3,5],
                ....:             [4,4,4,6,6,6],
                ....:             [5,6]])
                sage: b.string_parameters([1,2,1,3,2,1,4,3,2,1,5,4,3,2,1])
                [0, 1, 1, 1, 1, 0, 4, 4, 3, 0, 11, 10, 7, 7, 6]

                sage: B = crystals.infinity.Tableaux("G2")
                sage: b = B(rows=[[1,1,1,1,1,3,3,0,-3,-3,-2,-2,-1,-1,-1,-1],[2,3,3,3]])
                sage: b.string_parameters([2,1,2,1,2,1])
                [5, 13, 11, 15, 4, 4]
                sage: b.string_parameters([1,2,1,2,1,2])
                [7, 12, 15, 8, 10, 0]

            ::

                sage: C = crystals.Tableaux(['C',2], shape=[2,1])
                sage: mg = C.highest_weight_vector()
                sage: lw = C.lowest_weight_vectors()[0]
                sage: lw.string_parameters([1,2,1,2])
                [1, 2, 3, 1]
                sage: lw.string_parameters([2,1,2,1])
                [1, 3, 2, 1]
                sage: lw.e_string([2,1,1,1,2,2,1]) == mg
                True
                sage: lw.e_string([1,2,2,1,1,1,2]) == mg
                True

            TESTS::

                sage: B = crystals.infinity.NakajimaMonomials(['B',3])
                sage: mg = B.highest_weight_vector()
                sage: mg.string_parameters()
                [0, 0, 0, 0, 0, 0, 0, 0, 0]
                sage: w0 = WeylGroup(['B',3]).long_element().reduced_word()
                sage: def f_word(params):
                ....:     return reversed([index for i, index in enumerate(w0)
                ....:                      for _ in range(params[i])])
                sage: all(mg.f_string( f_word(x.value.string_parameters(w0)) ) == x.value
                ....:     for x in B.subcrystal(max_depth=4))
                True

                sage: B = crystals.infinity.NakajimaMonomials(['A',2,1])
                sage: mg = B.highest_weight_vector()
                sage: mg.string_parameters()
                Traceback (most recent call last):
                ...
                ValueError: the word must be specified because the
                 Weyl group is not finite
            """
            if word is None:
                if not self.cartan_type().is_finite():
                    raise ValueError("the word must be specified because"
                                     " the Weyl group is not finite")
                from sage.combinat.root_system.weyl_group import WeylGroup
                word = WeylGroup(self.cartan_type()).long_element().reduced_word()
            x = self
            params = []
            for i in word:
                count = 0
                y = x.e(i)
                while y is not None:
                    x = y
                    y = x.e(i)
                    count += 1
                params.append(count)
            return params

    class TensorProducts(TensorProductsCategory):
        """
        The category of highest weight crystals constructed by tensor
        product of highest weight crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: HighestWeightCrystals().TensorProducts().extra_super_categories()
                [Category of highest weight crystals]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            Implement operations on tensor products of crystals.
            """
            @cached_method
            def highest_weight_vectors(self):
                r"""
                Return the highest weight vectors of ``self``.

                This works by using a backtracing algorithm since if
                `b_2 \otimes b_1` is highest weight then `b_1` is
                highest weight.

                EXAMPLES::

                    sage: C = crystals.Tableaux(['D',4], shape=[2,2])
                    sage: D = crystals.Tableaux(['D',4], shape=[1])
                    sage: T = crystals.TensorProduct(D, C)
                    sage: T.highest_weight_vectors()
                    ([[[1]], [[1, 1], [2, 2]]],
                     [[[3]], [[1, 1], [2, 2]]],
                     [[[-2]], [[1, 1], [2, 2]]])
                    sage: L = filter(lambda x: x.is_highest_weight(), T)
                    sage: tuple(L) == T.highest_weight_vectors()
                    True

                TESTS:

                We check this works with Kashiwara's convention for
                tensor products::

                    sage: C = crystals.Tableaux(['B',3], shape=[2,2])
                    sage: D = crystals.Tableaux(['B',3], shape=[1])
                    sage: T = crystals.TensorProduct(D, C)
                    sage: T.options(convention='Kashiwara')
                    sage: T.highest_weight_vectors()
                    ([[[1, 1], [2, 2]], [[1]]],
                     [[[1, 1], [2, 2]], [[3]]],
                     [[[1, 1], [2, 2]], [[-2]]])
                    sage: T.options._reset()
                    sage: T.highest_weight_vectors()
                    ([[[1]], [[1, 1], [2, 2]]],
                     [[[3]], [[1, 1], [2, 2]]],
                     [[[-2]], [[1, 1], [2, 2]]])
                """
                return tuple(self.highest_weight_vectors_iterator())

            def highest_weight_vectors_iterator(self):
                r"""
                Iterate over the highest weight vectors of ``self``.

                This works by using a backtracing algorithm since if
                `b_2 \otimes b_1` is highest weight then `b_1` is
                highest weight.

                EXAMPLES::

                    sage: C = crystals.Tableaux(['D',4], shape=[2,2])
                    sage: D = crystals.Tableaux(['D',4], shape=[1])
                    sage: T = crystals.TensorProduct(D, C)
                    sage: tuple(T.highest_weight_vectors_iterator())
                    ([[[1]], [[1, 1], [2, 2]]],
                     [[[3]], [[1, 1], [2, 2]]],
                     [[[-2]], [[1, 1], [2, 2]]])
                    sage: L = filter(lambda x: x.is_highest_weight(), T)
                    sage: tuple(L) == tuple(T.highest_weight_vectors_iterator())
                    True

                TESTS:

                We check this works with Kashiwara's convention for
                tensor products::

                    sage: C = crystals.Tableaux(['B',3], shape=[2,2])
                    sage: D = crystals.Tableaux(['B',3], shape=[1])
                    sage: T = crystals.TensorProduct(D, C)
                    sage: T.options(convention='Kashiwara')
                    sage: tuple(T.highest_weight_vectors_iterator())
                    ([[[1, 1], [2, 2]], [[1]]],
                     [[[1, 1], [2, 2]], [[3]]],
                     [[[1, 1], [2, 2]], [[-2]]])
                    sage: T.options._reset()
                    sage: tuple(T.highest_weight_vectors_iterator())
                    ([[[1]], [[1, 1], [2, 2]]],
                     [[[3]], [[1, 1], [2, 2]]],
                     [[[-2]], [[1, 1], [2, 2]]])

                This currently is not implemented for infinite crystals::

                    sage: P = RootSystem(['A',3,1]).weight_lattice(extended=True)
                    sage: M = crystals.NakajimaMonomials(P.fundamental_weight(0))
                    sage: T = tensor([M, M])
                    sage: list(T.highest_weight_vectors_iterator())
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: not implemented for infinite crystals

                Check that :issue:`30493` is fixed::

                    sage: CW = CartanType("G", 2)
                    sage: C = crystals.Letters(CW)
                    sage: C.highest_weight_vectors()
                    (1,)
                    sage: T = crystals.TensorProduct(C)
                    sage: T.highest_weight_vectors()
                    ([1],)
                """
                if len(self.crystals) == 1:
                    for b in self.crystals[0].highest_weight_vectors():
                        yield self.element_class(self, [b])
                    return
                I = self.index_set()
                try:
                    T_elts = [C.list() for C in self.crystals[:-1]]
                except (TypeError, NotImplementedError, AttributeError):
                    raise NotImplementedError("not implemented for infinite crystals")
                from sage.categories.regular_crystals import RegularCrystals
                if self in RegularCrystals:
                    def hw_test(b2, i, d):
                        return d < 0
                else:
                    def hw_test(b2, i, d):
                        return d < 0 and b2.e(i) is not None
                T_len = [len(elts) for elts in T_elts]
                m = len(self.crystals) - 1
                for b in self.crystals[-1].highest_weight_vectors():
                    T_pos = m - 1  # current tensor position
                    T_cur = [0]*m  # index of current element for each tensor position
                    path = [None]*m + [b]
                    # cache phi for path up to current tensor position
                    T_phi = [None]*(m-1) + [{i: b.phi(i) for i in I}]
                    while T_pos < m:
                        if T_cur[T_pos] == T_len[T_pos]:
                            T_cur[T_pos] = 0
                            T_pos += 1
                            continue

                        b2 = T_elts[T_pos][T_cur[T_pos]]
                        T_cur[T_pos] += 1
                        b1_phi = T_phi[T_pos]
                        b1_phi_minus_b2_epsilon = {}
                        # break if (b2, b1) is not highest weight
                        for i in I:
                            d = b1_phi[i] - b2.epsilon(i)
                            # In the non-regular case, d may be nan.
                            # In this case b2.e(i) is None,
                            # and we may rely on max(0, nan) == 0.
                            # In the regular case, the next line is simply
                            #   if d < 0:
                            if hw_test(b2, i, d):
                                break
                            b1_phi_minus_b2_epsilon[i] = d
                        else:
                            path[T_pos] = b2
                            if T_pos:
                                T_pos -= 1
                                # In the regular case, the next line is simply
                                #   T_phi[T_pos] = {i: b2.phi(i) + b1_phi_minus_b2_epsilon[i] for i in I}
                                T_phi[T_pos] = {i: b2.phi(i) + max(0, b1_phi_minus_b2_epsilon[i])
                                                for i in I}
                            else:
                                yield self.element_class(self, path)

###############################################################################
## Morphisms


class HighestWeightCrystalMorphism(CrystalMorphismByGenerators):
    r"""
    A virtual crystal morphism whose domain is a highest weight crystal.

    INPUT:

    - ``parent`` -- a homset
    - ``on_gens`` -- a function or list that determines the image of the
      generators (if given a list, then this uses the order of the
      generators of the domain) of the domain under ``self``
    - ``cartan_type`` -- (optional) a Cartan type; the default is the
      Cartan type of the domain
    - ``virtualization`` -- (optional) a dictionary whose keys are
      in the index set of the domain and whose values are lists of
      entries in the index set of the codomain
    - ``scaling_factors`` -- (optional) a dictionary whose keys are in
      the index set of the domain and whose values are scaling factors
      for the weight, `\varepsilon` and `\varphi`
    - ``gens`` -- (optional) a list of generators to define the morphism;
      the default is to use the highest weight vectors of the crystal
    - ``check`` -- boolean (default: ``True``); check if the crystal morphism
      is valid
    """
    def __init__(self, parent, on_gens, cartan_type=None,
                 virtualization=None, scaling_factors=None,
                 gens=None, check=True):
        """
        Construct a crystal morphism.

        TESTS::

            sage: B = crystals.infinity.Tableaux(['B',2])
            sage: C = crystals.infinity.NakajimaMonomials(['B',2])
            sage: psi = B.crystal_morphism(C.module_generators)

            sage: B = crystals.Tableaux(['B',3], shape=[1])
            sage: C = crystals.Tableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
        """
        if cartan_type is None:
            cartan_type = parent.domain().cartan_type()
        if isinstance(on_gens, dict):
            gens = on_gens.keys()
        I = cartan_type.index_set()
        if gens is None:
            if cartan_type == parent.domain().cartan_type():
                gens = parent.domain().highest_weight_vectors()
            else:
                gens = tuple(x for x in parent.domain() if x.is_highest_weight(I))
            self._hw_gens = True
        elif check:
            self._hw_gens = all(x.is_highest_weight(I) for x in gens)
        else:
            self._hw_gens = False
        CrystalMorphismByGenerators.__init__(self, parent, on_gens, cartan_type,
                                             virtualization, scaling_factors,
                                             gens, check)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        TESTS::

            sage: B = crystals.infinity.Tableaux(['B',2])
            sage: C = crystals.infinity.NakajimaMonomials(['B',2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: b = B.highest_weight_vector()
            sage: psi(b)
            1
            sage: c = psi(b.f_string([1,1,1,2,2,1,2,2])); c
            Y(1,0)^-4 Y(2,0)^4 Y(2,1)^-4
            sage: c == C.highest_weight_vector().f_string([1,1,1,2,2,1,2,2])
            True

            sage: B = crystals.Tableaux(['B',3], shape=[1])
            sage: C = crystals.Tableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
            sage: psi(B.module_generators[0])
            [[1, 1]]

        We check with the morphism defined on the lowest weight vector::

            sage: B = crystals.Tableaux(['A',2], shape=[1])
            sage: La = RootSystem(['A',2]).weight_lattice().fundamental_weights()
            sage: T = crystals.elementary.T(['A',2], La[2])
            sage: Bp = T.tensor(B)
            sage: C = crystals.Tableaux(['A',2], shape=[2,1])
            sage: H = Hom(Bp, C)
            sage: x = C.module_generators[0].f_string([1,2])
            sage: psi = H({Bp.lowest_weight_vectors()[0]: x})
            sage: psi
            ['A', 2] Crystal morphism:
              From: Full tensor product of the crystals
                [The T crystal of type ['A', 2] and weight Lambda[2],
                 The crystal of tableaux of type ['A', 2] and shape(s) [[1]]]
              To:   The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
              Defn: [Lambda[2], [[3]]] |--> [[1, 3], [2]]
            sage: psi(Bp.highest_weight_vector())
            [[1, 1], [2]]
        """
        if not self._hw_gens:
            return CrystalMorphismByGenerators._call_(self, x)
        mg, path = x.to_highest_weight(self._cartan_type.index_set())
        cur = self._on_gens(mg)
        for i in reversed(path):
            if cur is None:
                return None
            s = []
            sf = self._scaling_factors[i]
            for j in self._virtualization[i]:
                s += [j]*sf
            cur = cur.f_string(s)
        return cur


class HighestWeightCrystalHomset(CrystalHomset):
    """
    The set of crystal morphisms from a highest weight crystal to
    another crystal.

    .. SEEALSO::

        See :class:`sage.categories.crystals.CrystalHomset` for more
        information.
    """
    def __init__(self, X, Y, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: B = crystals.Tableaux(['A', 2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: B = crystals.infinity.Tableaux(['B',2])
            sage: H = Hom(B, B)
        """
        if category is None:
            category = HighestWeightCrystals()
        CrystalHomset.__init__(self, X, Y, category)

    Element = HighestWeightCrystalMorphism
