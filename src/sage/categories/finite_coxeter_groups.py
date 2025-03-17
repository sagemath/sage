# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.combinat sage.groups
r"""
Finite Coxeter Groups
"""
# ****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#  Copyright (C) 2009    Nicolas Borie <nicolas dot borie at math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import sage.rings.abc

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.coxeter_groups import CoxeterGroups


class FiniteCoxeterGroups(CategoryWithAxiom):
    r"""
    The category of finite Coxeter groups.

    EXAMPLES::

        sage: CoxeterGroups.Finite()
        Category of finite Coxeter groups
        sage: FiniteCoxeterGroups().super_categories()
        [Category of finite generalized Coxeter groups,
         Category of Coxeter groups]

        sage: G = CoxeterGroups().Finite().example()
        sage: G.cayley_graph(side = "right").plot()
        Graphics object consisting of 40 graphics primitives

    Here are some further examples::

        sage: WeylGroups().Finite().example()
        The symmetric group on {0, ..., 3}

        sage: WeylGroup(["B", 3])
        Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space)

    Those other examples will eventually be also in this category::

        sage: SymmetricGroup(4)
        Symmetric group of order 4! as a permutation group
        sage: DihedralGroup(5)
        Dihedral group of order 10 as a permutation group
    """
    def extra_super_categories(self):
        r"""
        EXAMPLES::

            sage: CoxeterGroups().Finite().super_categories()
            [Category of finite generalized Coxeter groups,
             Category of Coxeter groups]
        """
        from sage.categories.complex_reflection_groups import ComplexReflectionGroups
        return [ComplexReflectionGroups().Finite().WellGenerated()]

    class ParentMethods:
        """
        Ambiguity resolution: the implementation of ``some_elements``
        is preferable to that of :class:`FiniteGroups`. The same holds
        for ``__iter__``, although a breadth first search would be more
        natural; at least this maintains backward compatibility after
        :issue:`13589`.

        TESTS::

            sage: W = FiniteCoxeterGroups().example(3)

            sage: W.some_elements.__module__
            'sage.categories.complex_reflection_or_generalized_coxeter_groups'
            sage: W.__iter__.__module__
            'sage.categories.coxeter_groups'

            sage: W.some_elements()
            [(1,), (2,), (), (1, 2)]
            sage: list(W)
            [(), (1,), (2,), (1, 2), (2, 1), (1, 2, 1)]
        """
        __iter__ = CoxeterGroups.ParentMethods.__dict__["__iter__"]

        @lazy_attribute
        def w0(self):
            r"""
            Return the longest element of ``self``.

            This attribute is deprecated, use :meth:`long_element` instead.

            EXAMPLES::

                sage: D8 = FiniteCoxeterGroups().example(8)
                sage: D8.w0
                (1, 2, 1, 2, 1, 2, 1, 2)
                sage: D3 = FiniteCoxeterGroups().example(3)
                sage: D3.w0
                (1, 2, 1)
            """
            return self.long_element()

        def long_element(self, index_set=None, as_word=False):
            r"""
            Return the longest element of ``self``, or of the
            parabolic subgroup corresponding to the given ``index_set``.

            INPUT:

            - ``index_set`` -- a subset (as a list or iterable) of the
              nodes of the Dynkin diagram; (default: all of them)

            - ``as_word`` -- boolean (default: ``False``); if ``True``, then
              return instead a reduced decomposition of the longest element

            Should this method be called maximal_element? longest_element?

            EXAMPLES::

                sage: D10 = FiniteCoxeterGroups().example(10)
                sage: D10.long_element()
                (1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
                sage: D10.long_element([1])
                (1,)
                sage: D10.long_element([2])
                (2,)
                sage: D10.long_element([])
                ()

                sage: D7 = FiniteCoxeterGroups().example(7)
                sage: D7.long_element()
                (1, 2, 1, 2, 1, 2, 1)

            One can require instead a reduced word for w0::

                sage: A3 = CoxeterGroup(['A', 3])
                sage: A3.long_element(as_word=True)
                [1, 2, 1, 3, 2, 1]
            """
            if index_set is None:
                index_set = self.index_set()
            w = self.one()
            if as_word:
                word = []
            while True:
                i = w.first_descent(index_set=index_set, positive=True)
                if i is None:
                    if as_word:
                        return word
                    else:
                        return w
                else:
                    if as_word:
                        word.append(i)
                    w = w.apply_simple_reflection(i)

        @cached_method
        def bruhat_poset(self, facade=False):
            """
            Return the Bruhat poset of ``self``.

            .. SEEALSO::

                :meth:`bhz_poset`, :meth:`shard_poset`, :meth:`weak_poset`

            EXAMPLES::

                sage: W = WeylGroup(["A", 2])
                sage: P = W.bruhat_poset()
                sage: P
                Finite poset containing 6 elements
                sage: P.show()

            Here are some typical operations on this poset::

                sage: W = WeylGroup(["A", 3])
                sage: P = W.bruhat_poset()
                sage: u = W.from_reduced_word([3,1])
                sage: v = W.from_reduced_word([3,2,1,2,3])
                sage: P(u) <= P(v)
                True
                sage: len(P.interval(P(u), P(v)))
                10
                sage: P.is_join_semilattice()
                False

            By default, the elements of `P` are aware that they belong
            to `P`::

                sage: P.an_element().parent()
                Finite poset containing 24 elements

            If instead one wants the elements to be plain elements of
            the Coxeter group, one can use the ``facade`` option::

                sage: P = W.bruhat_poset(facade = True)
                sage: P.an_element().parent()
                Weyl Group of type ['A', 3] (as a matrix group acting on the ambient space)

            .. SEEALSO:: :func:`Poset` for more on posets and facade posets.

            TESTS::

                sage: [len(WeylGroup(["A", n]).bruhat_poset().cover_relations()) for n in [1,2,3]]
                [1, 8, 58]

            .. TODO::

                - Use the symmetric group in the examples (for nicer
                  output), and print the edges for a stronger test.
                - The constructed poset should be lazy, in order to
                  handle large / infinite Coxeter groups.
            """
            from sage.combinat.posets.posets import Poset
            covers = tuple([u, v] for v in self for u in v.bruhat_lower_covers())
            return Poset((self, covers), cover_relations=True, facade=facade)

        def shard_poset(self, side='right'):
            """
            Return the shard intersection order attached to `W`.

            This is a lattice structure on `W`, introduced in [Rea2009]_. It
            contains the noncrossing partition lattice, as the induced lattice
            on the subset of `c`-sortable elements.

            The partial order is given by simultaneous inclusion of inversion sets
            and subgroups attached to every element.

            The precise description used here can be found in [STW2018]_.

            Another implementation for the symmetric groups is
            available as :func:`~sage.combinat.shard_order.shard_poset`.

            .. SEEALSO::

                :meth:`bhz_poset`, :meth:`bruhat_poset`, :meth:`weak_poset`

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
                sage: SH = W.shard_poset(); SH
                Finite lattice containing 24 elements
                sage: SH.is_graded()
                True
                sage: SH.characteristic_polynomial()
                q^3 - 11*q^2 + 23*q - 13
                sage: SH.f_polynomial()
                34*q^3 + 22*q^2 + q
            """
            from sage.combinat.posets.lattices import LatticePoset
            data = {w: (frozenset(u.lift()
                                  for u in w.covered_reflections_subgroup()),
                        frozenset((~w).inversions_as_reflections()))
                    for w in self}

            def shard_comparison(u, v):
                Gu, Nu = data[u]
                Gv, Nv = data[v]
                return Gu.issubset(Gv) and Nu.issubset(Nv)

            return LatticePoset([self, shard_comparison])

        def bhz_poset(self):
            r"""
            Return the Bergeron-Hohlweg-Zabrocki partial order on the Coxeter
            group.

            This is a partial order on the elements of a finite
            Coxeter group `W`, which is distinct from the Bruhat
            order, the weak order and the shard intersection order. It
            was defined in [BHZ2005]_.

            This partial order is not a lattice, as there is no unique
            maximal element. It can be succinctly defined as follows.

            Let `u` and `v` be two elements of the Coxeter group `W`. Let
            `S(u)` be the support of `u`. Then `u \leq v` if and only
            if `v_{S(u)} = u` (here `v = v^I v_I` denotes the usual
            parabolic decomposition with respect to the standard parabolic
            subgroup `W_I`).

            .. SEEALSO::

                :meth:`bruhat_poset`, :meth:`shard_poset`, :meth:`weak_poset`

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
                sage: P = W.bhz_poset(); P
                Finite poset containing 24 elements
                sage: P.relations_number()
                103
                sage: P.chain_polynomial()
                34*q^4 + 90*q^3 + 79*q^2 + 24*q + 1
                sage: len(P.maximal_elements())
                13
            """
            from sage.graphs.digraph import DiGraph
            from sage.combinat.posets.posets import Poset

            def covered_by(ux, vy):
                u, iu, Su = ux
                v, iv, Sv = vy
                if len(Sv) != len(Su) + 1:
                    return False
                if not all(u in Sv for u in Su):
                    return False
                return all((v * iu).has_descent(x, positive=True) for x in Su)

            vertices = [(u, u.inverse(),
                         tuple(set(u.reduced_word_reverse_iterator())))
                        for u in self]
            dg = DiGraph([vertices, covered_by])
            dg.relabel(lambda x: x[0])
            return Poset(dg, cover_relations=True)

        def degrees(self):
            """
            Return the degrees of the Coxeter group.

            The output is an increasing list of integers.

            EXAMPLES::

                sage: CoxeterGroup(['A', 4]).degrees()
                (2, 3, 4, 5)
                sage: CoxeterGroup(['B', 4]).degrees()
                (2, 4, 6, 8)
                sage: CoxeterGroup(['D', 4]).degrees()
                (2, 4, 4, 6)
                sage: CoxeterGroup(['F', 4]).degrees()
                (2, 6, 8, 12)
                sage: CoxeterGroup(['E', 8]).degrees()
                (2, 8, 12, 14, 18, 20, 24, 30)
                sage: CoxeterGroup(['H', 3]).degrees()
                (2, 6, 10)

                sage: WeylGroup([["A",3], ["A",3], ["B",2]]).degrees()
                (2, 3, 4, 2, 3, 4, 2, 4)

            TESTS::

                sage: CoxeterGroup(['A', 4]).degrees()
                (2, 3, 4, 5)
                sage: SymmetricGroup(3).degrees()
                (2, 3)
            """
            from sage.rings.qqbar import QQbar
            from sage.rings.integer_ring import ZZ

            def degrees_of_irreducible_component(I):
                """Return the degrees for the irreducible component indexed by I"""
                # A Coxeter element
                s = self.simple_reflections()
                c = self.prod(s[i] for i in I)
                roots = c.matrix().change_ring(QQbar).charpoly().roots()
                args = ((z.rational_argument(), m) for z, m in roots)
                args = [(z if z >= 0 else 1 + z, m) for z, m in args]
                h = max(z.denominator() for z, m in args)
                return tuple(sorted(ZZ(z * h + 1)
                                    for z, m in args if z
                                    for i in range(m)))

            return sum((degrees_of_irreducible_component(I)
                        for I in self.irreducible_component_index_sets()), ())

        def codegrees(self):
            """
            Return the codegrees of the Coxeter group.

            These are just the degrees minus 2.

            EXAMPLES::

                sage: CoxeterGroup(['A', 4]).codegrees()
                (0, 1, 2, 3)
                sage: CoxeterGroup(['B', 4]).codegrees()
                (0, 2, 4, 6)
                sage: CoxeterGroup(['D', 4]).codegrees()
                (0, 2, 2, 4)
                sage: CoxeterGroup(['F', 4]).codegrees()
                (0, 4, 6, 10)
                sage: CoxeterGroup(['E', 8]).codegrees()
                (0, 6, 10, 12, 16, 18, 22, 28)
                sage: CoxeterGroup(['H', 3]).codegrees()
                (0, 4, 8)

                sage: WeylGroup([["A",3], ["A",3], ["B",2]]).codegrees()
                (0, 1, 2, 0, 1, 2, 0, 2)
            """
            return tuple(d - 2 for d in self.degrees())

        @cached_method
        def weak_poset(self, side='right', facade=False):
            """
            INPUT:

            - ``side`` -- "left", "right", or "twosided" (default: ``'right'``)
            - ``facade`` -- boolean (default: ``False``)

            Returns the left (resp. right) poset for weak order.  In
            this poset, `u` is smaller than `v` if some reduced word
            of `u` is a right (resp. left) factor of some reduced word
            of `v`.

            .. SEEALSO::

                :meth:`bhz_poset`, :meth:`bruhat_poset`, :meth:`shard_poset`

            EXAMPLES::

                sage: W = WeylGroup(["A", 2])
                sage: P = W.weak_poset()
                sage: P
                Finite lattice containing 6 elements
                sage: P.show()

            This poset is in fact a lattice::

                sage: W = WeylGroup(["B", 3])
                sage: P = W.weak_poset(side = "left")
                sage: P.is_lattice()
                True

            so this method has an alias :meth:`weak_lattice`::

                sage: W.weak_lattice(side = "left") is W.weak_poset(side = "left")
                True

            As a bonus feature, one can create the left-right weak
            poset::

                sage: W = WeylGroup(["A",2])
                sage: P = W.weak_poset(side = "twosided")
                sage: P.show()
                sage: len(P.hasse_diagram().edges(sort=False))
                8

            This is the transitive closure of the union of left and
            right order. In this poset, `u` is smaller than `v` if
            some reduced word of `u` is a factor of some reduced word
            of `v`. Note that this is not a lattice::

                sage: P.is_lattice()
                False

            By default, the elements of `P` are aware of that they
            belong to `P`::

                sage: P.an_element().parent()
                Finite poset containing 6 elements

            If instead one wants the elements to be plain elements of
            the Coxeter group, one can use the ``facade`` option::

                sage: P = W.weak_poset(facade = True)
                sage: P.an_element().parent()
                Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)

            .. SEEALSO:: :func:`Poset` for more on posets and facade posets.

            TESTS::

                sage: [len(WeylGroup(["A", n]).weak_poset(side = "right").cover_relations()) for n in [1,2,3]]
                [1, 6, 36]
                sage: [len(WeylGroup(["A", n]).weak_poset(side = "left" ).cover_relations()) for n in [1,2,3]]
                [1, 6, 36]

            .. TODO::

                - Use the symmetric group in the examples (for nicer
                  output), and print the edges for a stronger test.
                - The constructed poset should be lazy, in order to
                  handle large / infinite Coxeter groups.
            """
            from sage.combinat.posets.posets import Poset
            from sage.combinat.posets.lattices import LatticePoset
            if side == "twosided":
                covers = tuple([u, v] for u in self for v in u.upper_covers(side='left') + u.upper_covers(side='right'))
                return Poset((self, covers), cover_relations=True,
                             facade=facade)
            covers = tuple([u, v] for u in self for v in u.upper_covers(side=side))
            return LatticePoset((self, covers), cover_relations=True,
                                facade=facade)

        weak_lattice = weak_poset

        def inversion_sequence(self, word):
            r"""
            Return the inversion sequence corresponding to the ``word``
            in indices of simple generators of ``self``.

            If ``word`` corresponds to `[w_0,w_1,...w_k]`, the output is
            `[w_0,w_0w_1w_0,\ldots,w_0w_1\cdots w_k \cdots w_1 w_0]`.

            INPUT:

            - ``word`` -- a word in the indices of the simple
              generators of ``self``

            EXAMPLES::

                sage: CoxeterGroup(["A", 2]).inversion_sequence([1,2,1])
                [
                [-1  1]  [ 0 -1]  [ 1  0]
                [ 0  1], [-1  0], [ 1 -1]
                ]

                sage: [t.reduced_word() for t in CoxeterGroup(["A",3]).inversion_sequence([2,1,3,2,1,3])]
                [[2], [1, 2, 1], [2, 3, 2], [1, 2, 3, 2, 1], [3], [1]]
            """
            return [self.from_reduced_word(word[:i+1]+list(reversed(word[:i])))
                    for i in range(len(word))]

        def reflections_from_w0(self):
            """
            Return the reflections of ``self`` using the inversion set
            of ``w_0``.

            EXAMPLES::

                sage: WeylGroup(['A',2]).reflections_from_w0()
                [
                [0 1 0]  [0 0 1]  [1 0 0]
                [1 0 0]  [0 1 0]  [0 0 1]
                [0 0 1], [1 0 0], [0 1 0]
                ]

                sage: WeylGroup(['A',3]).reflections_from_w0()
                [
                [0 1 0 0]  [0 0 1 0]  [1 0 0 0]  [0 0 0 1]  [1 0 0 0]  [1 0 0 0]
                [1 0 0 0]  [0 1 0 0]  [0 0 1 0]  [0 1 0 0]  [0 0 0 1]  [0 1 0 0]
                [0 0 1 0]  [1 0 0 0]  [0 1 0 0]  [0 0 1 0]  [0 0 1 0]  [0 0 0 1]
                [0 0 0 1], [0 0 0 1], [0 0 0 1], [1 0 0 0], [0 1 0 0], [0 0 1 0]
                ]
            """
            return self.long_element().inversions_as_reflections()

        @cached_method
        def m_cambrian_lattice(self, c, m=1, on_roots=False):
            """
            Return the `m`-Cambrian lattice on `m`-delta sequences.

            See :arxiv:`1503.00710` and :arxiv:`math/0611106`.

            The `m`-delta sequences are certain `m`-colored minimal
            factorizations of `c` into reflections.

            INPUT:

            - ``c`` -- a Coxeter element of ``self`` (as a tuple, or
              as an element of ``self``)

            - ``m`` -- positive integer (default: 1)

            - ``on_roots`` -- boolean (default: ``False``); if
              ``on_roots`` is ``True``, the lattice is realized on
              roots rather than on reflections. In order for this to
              work, the ElementMethod ``reflection_to_root`` must be
              available.

            EXAMPLES::

                sage: CoxeterGroup(["A",2]).m_cambrian_lattice((1,2))
                Finite lattice containing 5 elements

                sage: CoxeterGroup(["A",2]).m_cambrian_lattice((1,2),2)
                Finite lattice containing 12 elements
            """
            from sage.combinat.posets.lattices import LatticePoset
            if hasattr(c, "reduced_word"):
                c = c.reduced_word()
            c = list(c)

            sorting_word = self.long_element().coxeter_sorting_word(c)

            if on_roots:
                if not hasattr(self.long_element(), "reflection_to_root"):
                    raise ValueError("The parameter 'on_root=True' needs "
                                     "the ElementMethod 'reflection_to_root'")

                inv_woc = [t.reflection_to_root()
                           for t in self.inversion_sequence(sorting_word)]
                S = [s.reflection_to_root() for s in self.simple_reflections()]
                PhiP = [t.reflection_to_root() for t in self.reflections()]
            else:
                inv_woc = self.inversion_sequence(sorting_word)
                S = self.simple_reflections()
                T = self.reflections_from_w0()
                Twords = {t: t.reduced_word() for t in T}

            elements = set()
            covers = []

            bottom_elt = frozenset((s, 0) for s in S)
            new = {bottom_elt}
            while new:
                new_element = new.pop()
                elements.add(new_element)
                for t in new_element:
                    if t[1] < m:
                        cov_element = [s for s in new_element if s != t]
                        cov_element.append((t[0], t[1] + 1))
                        idx_t0 = inv_woc.index(t[0])
                        for t_conj in [(i, t[1]) for i in inv_woc[idx_t0:]] + [(i, t[1] + 1) for i in inv_woc[:idx_t0]]:
                            if t_conj in cov_element:
                                cov_element.remove(t_conj)
                                if on_roots:
                                    tmp = t_conj[0].weyl_action(t[0].associated_reflection())
                                    if tmp in PhiP:
                                        cov_element.append((tmp, t_conj[1]))
                                    else:
                                        cov_element.append((-tmp, t_conj[1] - 1))
                                else:
                                    tmp = t[0] * t_conj[0] * t[0]
                                    invs = self.inversion_sequence(Twords[t[0]]+Twords[t_conj[0]])
                                    plus_or_minus = invs.count(tmp)
                                    if plus_or_minus % 2:
                                        cov_element.append((tmp, t_conj[1]))
                                    else:
                                        cov_element.append((tmp, t_conj[1] - 1))

                        cov_element = frozenset(cov_element)
                        if cov_element not in elements:
                            new.add(cov_element)
                        covers.append((new_element, cov_element))
            return LatticePoset([elements, covers], cover_relations=True)

        def cambrian_lattice(self, c, on_roots=False):
            """
            Return the `c`-Cambrian lattice on delta sequences.

            See :arxiv:`1503.00710` and :arxiv:`math/0611106`.

            Delta sequences are certain 2-colored minimal factorizations
            of `c` into reflections.

            INPUT:

            - ``c`` -- a standard Coxeter element in ``self``
              (as a tuple, or as an element of ``self``)

            - ``on_roots`` -- boolean (default: ``False``); if
              ``on_roots`` is ``True``, the lattice is realized on
              roots rather than on reflections. In order for this to
              work, the ElementMethod ``reflection_to_root`` must be
              available.

            EXAMPLES::

                sage: CoxeterGroup(["A", 2]).cambrian_lattice((1,2))
                Finite lattice containing 5 elements

                sage: CoxeterGroup(["B", 2]).cambrian_lattice((1,2))
                Finite lattice containing 6 elements

                sage: CoxeterGroup(["G", 2]).cambrian_lattice((1,2))
                Finite lattice containing 8 elements
            """
            return self.m_cambrian_lattice(c=c, m=1, on_roots=on_roots)

        def is_real(self):
            """
            Return ``True`` since ``self`` is a real reflection group.

            EXAMPLES::

                sage: CoxeterGroup(['F',4]).is_real()
                True
                sage: CoxeterGroup(['H',4]).is_real()
                True
            """
            return True

        def permutahedron(self, point=None, base_ring=None):
            r"""
            Return the permutahedron of ``self``.

            This is the convex hull of the point ``point`` in the weight
            basis under the action of ``self`` on the underlying vector
            space `V`.

            .. SEEALSO::

                :meth:`~sage.combinat.root_system.reflection_group_real.permutahedron`

            INPUT:

            - ``point`` -- (optional) a point given by its coordinates in
              the weight basis (default: `(1, 1, 1, \ldots)`)
            - ``base_ring`` -- (optional) the base ring of the polytope

            .. NOTE::

                The result is expressed in the root basis coordinates.

            .. NOTE::

                If function is too slow, switching the base ring to
                :class:`RDF` will almost certainly speed things up.

            EXAMPLES::

                sage: W = CoxeterGroup(['H',3], base_ring=RDF)
                sage: W.permutahedron()
                doctest:warning
                ...
                UserWarning: This polyhedron data is numerically complicated; cdd could not convert between the inexact V and H representation without loss of data. The resulting object might show inconsistencies.
                A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 120 vertices

                sage: W = CoxeterGroup(['I',7])
                sage: W.permutahedron()
                A 2-dimensional polyhedron in AA^2 defined as the convex hull of 14 vertices
                sage: W.permutahedron(base_ring=RDF)
                A 2-dimensional polyhedron in RDF^2 defined as the convex hull of 14 vertices

                sage: W = ReflectionGroup(['A',3])                          # optional - gap3
                sage: W.permutahedron()                                     # optional - gap3
                A 3-dimensional polyhedron in QQ^3 defined as the convex hull
                of 24 vertices

                sage: W = ReflectionGroup(['A',3],['B',2])                  # optional - gap3
                sage: W.permutahedron()                                     # optional - gap3
                A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 192 vertices

            TESTS::

                sage: W = ReflectionGroup(['A',3])                          # optional - gap3
                sage: W.permutahedron([3,5,8])                              # optional - gap3
                A 3-dimensional polyhedron in QQ^3 defined as the convex hull
                of 24 vertices


            .. PLOT::
                :width: 300 px

                W = CoxeterGroup(['I',7])
                p = W.permutahedron()
                sphinx_plot(p)
            """
            n = self.one().canonical_matrix().rank()
            weights = self.fundamental_weights()
            if point is None:
                from sage.rings.integer_ring import ZZ
                point = [ZZ.one()] * n
            v = sum(point[i-1] * weights[i] for i in weights.keys())
            vertices = [v*w for w in self]
            if base_ring is None:
                if isinstance(v.base_ring(), (sage.rings.abc.UniversalCyclotomicField,
                                              sage.rings.abc.AlgebraicField_common)):
                    from sage.rings.qqbar import AA
                    vertices = [v.change_ring(AA) for v in vertices]
                    base_ring = AA
            from sage.geometry.polyhedron.constructor import Polyhedron
            return Polyhedron(vertices=vertices, base_ring=base_ring)

        def coxeter_poset(self):
            r"""
            Return the Coxeter poset of ``self``.

            Let `W` be a Coxeter group. The *Coxeter poset* is defined as
            the set of (right) standard cosets `gW_J`, where `J` is a
            subset of the index set `I` of `W`, ordered by reverse inclusion.

            This is equal to the face poset of the :meth:`Coxeter complex
            <coxeter_complex()>`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3])
                sage: P = W.coxeter_poset()
                sage: P
                Finite meet-semilattice containing 75 elements
                sage: P.rank()
                3

                sage: W = WeylGroup(['B', 3])
                sage: P = W.coxeter_poset()
                sage: P
                Finite meet-semilattice containing 147 elements
                sage: P.rank()
                3

                sage: W = CoxeterGroup(['I', 7])
                sage: P = W.coxeter_poset()
                sage: P
                Finite meet-semilattice containing 29 elements
                sage: P.rank()
                2

                sage: W = CoxeterGroup(['H', 3])
                sage: P = W.coxeter_poset()
                sage: P
                Finite meet-semilattice containing 363 elements
                sage: P.rank()
                3

                sage: # optional - gap3
                sage: W = CoxeterGroup(['H', 3], implementation='permutation')
                sage: P = W.coxeter_poset()
                sage: P
                Finite meet-semilattice containing 363 elements
                sage: P.rank()
                3
            """
            I = self.index_set()
            data = {}
            next_level = {(g, ()) for g in self}
            while next_level:
                cur = next_level
                next_level = set()
                for Y in cur:
                    g, J = Y
                    for i in I:
                        if i in J:
                            continue
                        Jp = tuple(sorted(J + (i,)))
                        gp = g.coset_representative(Jp, side='right')
                        X = (gp, Jp)
                        if X in data:
                            data[X].append(Y)
                        else:
                            data[X] = [Y]
                        next_level.add(X)
            from sage.combinat.posets.lattices import MeetSemilattice
            return MeetSemilattice(data)

        def coxeter_complex(self):
            r"""
            Return the Coxeter complex of ``self``.

            Let `W` be a Coxeter group, and let `X` be the corresponding Tits
            cone, which is constructed as the `W` orbit of the fundamental
            chamber in the reflection representation. The *Coxeter complex*
            of `W` is the simplicial complex `(X \setminus \{0\}) / \RR_{>0}`.
            The face poset of this simplicial complex is given by the
            :meth:`coxeter_poset()`. When `W` is a finite group, then the
            Coxeter complex is homeomorphic to an `(n-1)`-dimensional sphere,
            where `n` is the rank of `W`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A', 3])
                sage: C = W.coxeter_complex()
                sage: C
                Simplicial complex with 14 vertices and 24 facets
                sage: C.homology()
                {0: 0, 1: 0, 2: Z}

                sage: W = WeylGroup(['B', 3])
                sage: C = W.coxeter_complex()
                sage: C
                Simplicial complex with 26 vertices and 48 facets
                sage: C.homology()
                {0: 0, 1: 0, 2: Z}

                sage: W = CoxeterGroup(['I', 7])
                sage: C = W.coxeter_complex()
                sage: C
                Simplicial complex with 14 vertices and 14 facets
                sage: C.homology()
                {0: 0, 1: Z}

                sage: W = CoxeterGroup(['H', 3])
                sage: C = W.coxeter_complex()
                sage: C
                Simplicial complex with 62 vertices and 120 facets
                sage: C.homology()
                {0: 0, 1: 0, 2: Z}

                sage: # optional - gap3
                sage: W = CoxeterGroup(['H', 3], implementation='permutation')
                sage: C = W.coxeter_complex()
                sage: C
                Simplicial complex with 62 vertices and 120 facets
                sage: C.homology()
                {0: 0, 1: 0, 2: Z}
            """
            I = self.index_set()
            facets = {}
            Ip = {i: tuple([j for j in I if j != i]) for i in I}
            for g in self:
                V = []
                for i in I:
                    gp = g
                    D = gp.descents(side='right')
                    if D and D[0] == i:
                        D.pop(0)
                    while D:
                        gp = gp.apply_simple_reflection(D[0])
                        D = gp.descents(side='right')
                        if D and D[0] == i:
                            D.pop(0)
                    V.append((gp, Ip[i]))
                facets[g] = V
            verts = set()
            for F in facets.values():
                verts.update(F)
            labels = {x: i for i, x in enumerate(verts)}
            result = [[labels[v] for v in F] for F in facets.values()]
            from sage.topology.simplicial_complex import SimplicialComplex
            return SimplicialComplex(result)

    class ElementMethods:

        def absolute_length(self):
            """
            Return the absolute length of ``self``.

            The absolute length is the length of the shortest expression
            of the element as a product of reflections. For finite Coxeter
            groups, the absolute length is the codimension of the
            1-eigenspace of the element (Lemmas 1-3 in [Car1972a]_).

            For permutations in the symmetric groups, the absolute
            length is the size minus the number of its disjoint
            cycles.

            .. SEEALSO::

                :meth:`~sage.categories.coxeter_groups.absolute_le`

            EXAMPLES::

                sage: W = WeylGroup(["A", 3])                                           # needs sage.combinat sage.groups
                sage: s = W.simple_reflections()                                        # needs sage.combinat sage.groups
                sage: (s[1]*s[2]*s[3]).absolute_length()                                # needs sage.combinat sage.groups
                3

                sage: W = SymmetricGroup(4)                                             # needs sage.groups
                sage: s = W.simple_reflections()                                        # needs sage.groups
                sage: (s[3]*s[2]*s[1]).absolute_length()                                # needs sage.combinat sage.groups
                3
            """
            M = self.canonical_matrix()
            return (M - 1).image().dimension()

        @cached_in_parent_method
        def bruhat_upper_covers(self):
            r"""
            Return all the elements that cover ``self`` in Bruhat order.

            EXAMPLES::

                sage: W = WeylGroup(["A",4])
                sage: w = W.from_reduced_word([3,2])
                sage: print([v.reduced_word() for v in w.bruhat_upper_covers()])
                [[4, 3, 2], [3, 4, 2], [2, 3, 2], [3, 1, 2], [3, 2, 1]]

                sage: W = WeylGroup(["B",6])
                sage: w = W.from_reduced_word([1,2,1,4,5])
                sage: C = w.bruhat_upper_covers()
                sage: len(C)
                9
                sage: print([v.reduced_word() for v in C])
                [[6, 4, 5, 1, 2, 1], [4, 5, 6, 1, 2, 1], [3, 4, 5, 1, 2, 1], [2, 3, 4, 5, 1, 2],
                [1, 2, 3, 4, 5, 1], [4, 5, 4, 1, 2, 1], [4, 5, 3, 1, 2, 1], [4, 5, 2, 3, 1, 2],
                [4, 5, 1, 2, 3, 1]]
                sage: ww = W.from_reduced_word([5,6,5])
                sage: CC = ww.bruhat_upper_covers()
                sage: print([v.reduced_word() for v in CC])
                [[6, 5, 6, 5], [4, 5, 6, 5], [5, 6, 4, 5], [5, 6, 5, 4], [5, 6, 5, 3], [5, 6, 5, 2],
                [5, 6, 5, 1]]

            Recursive algorithm: write `w` for ``self``. If `i` is a
            non-descent of `w`, then the covers of `w` are exactly
            `\{ws_i, u_1s_i, u_2s_i,..., u_js_i\}`, where the `u_k`
            are those covers of `ws_i` that have a descent at `i`.
            """

            i = self.first_descent(positive=True, side='right')
            if i is not None:
                wsi = self.apply_simple_reflection(i, side='right')
                return [u.apply_simple_reflection(i, side='right')
                        for u in wsi.bruhat_upper_covers()
                        if u.has_descent(i, side='right')] + [wsi]
            return []

        def coxeter_knuth_neighbor(self, w):
            r"""
            Return the Coxeter-Knuth (oriented) neighbors of the reduced word `w` of ``self``.

            INPUT:

            - ``w`` -- reduced word of ``self``

            The Coxeter-Knuth relations are given by `a a+1 a \sim a+1 a a+1`, `abc \sim acb`
            if `b<a<c` and `abc \sim bac` if `a<c<b`. This method returns all neighbors of
            ``w`` under the Coxeter-Knuth relations oriented from left to right.

            EXAMPLES::

                sage: W = WeylGroup(['A',4], prefix='s')
                sage: word = [1,2,1,3,2]
                sage: w = W.from_reduced_word(word)
                sage: w.coxeter_knuth_neighbor(word)
                {(1, 2, 3, 1, 2), (2, 1, 2, 3, 2)}

                sage: word = [1,2,1,3,2,4,3]
                sage: w = W.from_reduced_word(word)
                sage: w.coxeter_knuth_neighbor(word)
                {(1, 2, 1, 3, 4, 2, 3), (1, 2, 3, 1, 2, 4, 3), (2, 1, 2, 3, 2, 4, 3)}

            TESTS::

                sage: W = WeylGroup(['B',4], prefix='s')
                sage: word = [1,2]
                sage: w = W.from_reduced_word(word)
                sage: w.coxeter_knuth_neighbor(word)
                Traceback (most recent call last):
                ...
                NotImplementedError: this has only been implemented in finite type A so far
            """
            C = self.parent().cartan_type()
            if not C[0] == 'A':
                raise NotImplementedError("this has only been implemented in finite type A so far")
            d = []
            for i in range(2, len(w)):
                v = list(w)
                if w[i-2] == w[i]:
                    if w[i] == w[i-1] - 1:
                        v[i-2] = w[i-1]
                        v[i] = w[i-1]
                        v[i-1] = w[i]
                        d += [tuple(v)]
                elif w[i-1] < w[i-2] and w[i-2] < w[i]:
                    v[i] = w[i-1]
                    v[i-1] = w[i]
                    d += [tuple(v)]
                elif w[i-2] < w[i] and w[i] < w[i-1]:
                    v[i-2] = w[i-1]
                    v[i-1] = w[i-2]
                    d += [tuple(v)]
            return set(d)

        def coxeter_knuth_graph(self):
            r"""
            Return the Coxeter-Knuth graph of type `A`.

            The Coxeter-Knuth graph of type `A` is generated by the Coxeter-Knuth relations which are
            given by `a a+1 a \sim a+1 a a+1`, `abc \sim acb` if `b<a<c` and `abc \sim bac` if `a<c<b`.

            EXAMPLES::

                sage: W = WeylGroup(['A',4], prefix='s')
                sage: w = W.from_reduced_word([1,2,1,3,2])
                sage: D = w.coxeter_knuth_graph()
                sage: D.vertices(sort=True)
                [(1, 2, 1, 3, 2),
                (1, 2, 3, 1, 2),
                (2, 1, 2, 3, 2),
                (2, 1, 3, 2, 3),
                (2, 3, 1, 2, 3)]
                sage: D.edges(sort=True)
                [((1, 2, 1, 3, 2), (1, 2, 3, 1, 2), None),
                ((1, 2, 1, 3, 2), (2, 1, 2, 3, 2), None),
                ((2, 1, 2, 3, 2), (2, 1, 3, 2, 3), None),
                ((2, 1, 3, 2, 3), (2, 3, 1, 2, 3), None)]

                sage: w = W.from_reduced_word([1,3])
                sage: D = w.coxeter_knuth_graph()
                sage: D.vertices(sort=True)
                [(1, 3), (3, 1)]
                sage: D.edges(sort=False)
                []

            TESTS::

                sage: W = WeylGroup(['B',4], prefix='s')
                sage: w = W.from_reduced_word([1,2])
                sage: w.coxeter_knuth_graph()
                Traceback (most recent call last):
                ...
                NotImplementedError: this has only been implemented in finite type A so far
            """
            from sage.graphs.graph import Graph
            R = [tuple(v) for v in self.reduced_words()]
            G = Graph()
            G.add_vertices(R)
            G.add_edges([v, vp] for v in R for vp in self.coxeter_knuth_neighbor(v))
            return G

        def is_coxeter_element(self):
            r"""
            Return whether this is a Coxeter element.

            This is, whether ``self`` has an eigenvalue `e^{2\pi i/h}`
            where `h` is the Coxeter number.

            .. SEEALSO:: :meth:`~sage.categories.finite_complex_reflection_groups.coxeter_elements`

            EXAMPLES::

                sage: W = CoxeterGroup(['A',2])
                sage: c = prod(W.gens())
                sage: c.is_coxeter_element()
                True
                sage: W.one().is_coxeter_element()
                False

                sage: W = WeylGroup(['G', 2])
                sage: c = prod(W.gens())
                sage: c.is_coxeter_element()
                True
                sage: W.one().is_coxeter_element()
                False
            """
            return self in self.parent().coxeter_elements()

        def covered_reflections_subgroup(self):
            """
            Return the subgroup of `W` generated by the conjugates by `w`
            of the simple reflections indexed by right descents of `w`.

            This is used to compute the shard intersection order on `W`.

            EXAMPLES::

                sage: W = CoxeterGroup(['A',3], base_ring=ZZ)
                sage: len(W.long_element().covered_reflections_subgroup())
                24
                sage: s = W.simple_reflection(1)
                sage: Gs = s.covered_reflections_subgroup()
                sage: len(Gs)
                2
                sage: s in [u.lift() for u in Gs]
                True
                sage: len(W.one().covered_reflections_subgroup())
                1
            """
            W = self.parent()
            winv = ~self
            cov_down = [self * W.simple_reflection(i) * winv
                        for i in self.descents(side='right')]
            return W.submonoid(cov_down)
