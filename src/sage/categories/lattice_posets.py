r"""
Lattice posets
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiéry <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.posets import Posets
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport


class LatticePosets(Category):
    r"""
    The category of lattices, i.e. partially ordered sets in which any
    two elements have a unique supremum (the elements' least upper
    bound; called their *join*) and a unique infimum (greatest lower bound;
    called their *meet*).

    EXAMPLES::

        sage: LatticePosets()
        Category of lattice posets
        sage: LatticePosets().super_categories()
        [Category of posets]
        sage: LatticePosets().example()
        NotImplemented

    .. SEEALSO::

        - :class:`~sage.categories.posets.Posets`
        - :class:`FiniteLatticePosets`, :func:`LatticePoset`

    TESTS::

        sage: C = LatticePosets()
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self) -> list:
        r"""
        Return a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: LatticePosets().super_categories()
            [Category of posets]
        """
        return [Posets()]

    class ParentMethods:

        @abstract_method
        def meet(self, x, y):
            """
            Return the meet of `x` and `y` in this lattice.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: D = LatticePoset((divisors(30), attrcall("divides")))             # needs sage.graphs sage.modules
                sage: D.meet( D(6), D(15) )                                             # needs sage.graphs sage.modules
                3
            """

        @abstract_method
        def join(self, x, y):
            """
            Return the join of `x` and `y` in this lattice.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: D = LatticePoset((divisors(60), attrcall("divides")))             # needs sage.graphs sage.modules
                sage: D.join( D(6), D(10) )                                             # needs sage.graphs sage.modules
                30
            """

    class SubcategoryMethods:
        def ChainGraded(self):
            r"""
            A lattice is graded if all maximal chains have the same length.

            To avoid possible confusion, the name of the axiom
            is ``ChainGraded``.

            EXAMPLES::

                sage: P = posets.DivisorLattice(24)
                sage: P in FiniteLatticePosets().ChainGraded()
                True
            """
            return self._with_axiom("ChainGraded")

        def Stone(self):
            r"""
            A Stone lattice `(L, \vee, \wedge)` is a pseudo-complemented
            distributive lattice such that `a^* \vee a^{**} = 1`.

            See :wikipedia:`Stone algebra`.

            EXAMPLES::

                sage: P = posets.DivisorLattice(24)
                sage: P in FiniteLatticePosets().Stone()
                True
            """
            return self._with_axiom("Stone")

        def Distributive(self):
            r"""
            A lattice `(L, \vee, \wedge)` is distributive if meet
            distributes over join: `x \wedge (y \vee z) = (x \wedge y)
            \vee (x \wedge z)` for every `x,y,z \in L`.

            From duality in lattices, it follows that then also join
            distributes over meet.

            A distributive lattice is always graded.

            See :wikipedia:`Distributive lattice`.

            EXAMPLES::

                sage: P = posets.ChainPoset(2).order_ideals_lattice()
                sage: P in FiniteLatticePosets().Distributive()
                True
            """
            return self._with_axiom("Trim")._with_axiom("ChainGraded")

        def CongruenceUniform(self):
            r"""
            A finite lattice `(L, \vee, \wedge)` is congruence uniform if it
            can be constructed by a sequence of interval doublings
            starting with the lattice with one element.

            EXAMPLES::

                sage: P = posets.TamariLattice(2)
                sage: P in FiniteLatticePosets().CongruenceUniform()
                True
            """
            return self._with_axiom("CongruenceUniform")

        def Semidistributive(self):
            r"""
            A finite lattice `(L, \vee, \wedge)` is semidistributive if
            it is both join-semidistributive and meet-semidistributive.

            A finite lattice is join-semidistributive if
            for all elements `e, x, y` in the lattice we have

            .. MATH::

                e \vee x = e \vee y \implies e \vee x = e \vee (x \wedge y)

            Meet-semidistributivity is the dual property.

            EXAMPLES::

                sage: P = posets.TamariLattice(2)
                sage: P in FiniteLatticePosets().Semidistributive()
                True
            """
            return self._with_axiom("Semidistributive")

        def Trim(self):
            r"""
            A finite lattice `(L, \vee, \wedge)` is trim if it is extremal
            and left modular.

            This notion is defined in [Thom2006]_.

            EXAMPLES::

                sage: P = posets.TamariLattice(2)
                sage: P in FiniteLatticePosets().Trim()
                True
            """
            return self._with_axiom("Trim")

        def Extremal(self):
            r"""
            A finite lattice `(L, \vee, \wedge)` is extremal if
            if it has a chain of length `n` (containing `n+1` elements)
            and exactly `n` join-irreducibles and `n` meet-irreducibles.

            This notion was defined by George Markowsky.

            EXAMPLES::

                sage: P = posets.TamariLattice(2)
                sage: P in FiniteLatticePosets().Extremal()
                True
            """
            return self._with_axiom("Extremal")

    Finite = LazyImport('sage.categories.finite_lattice_posets',
                        'FiniteLatticePosets')

    class Extremal(CategoryWithAxiom):
        """
        The category of extremal uniform lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().Extremal(); cat
            Category of finite extremal lattice posets

            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of extremal lattice posets]
        """
        class ParentMethods:
            def is_extremal(self) -> bool:
                """
                Return whether ``self`` is an extremal lattice.

                EXAMPLES::

                    sage: posets.TamariLattice(4).is_extremal()
                    True
                """
                return True

    class Trim(CategoryWithAxiom):
        """
        The category of trim uniform lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().Trim(); cat
            Category of finite trim lattice posets
            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of trim lattice posets]
        """
        @cached_method
        def extra_super_categories(self) -> list:
            r"""
            Return a list of the super categories of ``self``.

            These encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().Trim().super_categories()
                [Category of finite lattice posets,
                 Category of trim lattice posets]
            """
            return [LatticePosets().Extremal()]

        class ParentMethods:
            def is_trim(self) -> bool:
                """
                Return whether ``self`` is a trim lattice.

                EXAMPLES::

                    sage: posets.TamariLattice(4).is_trim()
                    True
                """
                return True

    class Semidistributive(CategoryWithAxiom):
        """
        The category of semidistributive lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().Semidistributive(); cat
            Category of finite semidistributive lattice posets

            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of semidistributive lattice posets]
        """
        class ParentMethods:
            def is_semidistributive(self) -> bool:
                """
                Return whether ``self`` is a semidistributive lattice.

                EXAMPLES::

                    sage: posets.TamariLattice(4).is_semidistributive()
                    True
                """
                return True

            def kappa(self, a):
                r"""
                Return the maximum element greater than the element covered
                by ``a`` but not greater than ``a``.

                Define `\kappa(a)` as the maximum element of
                `(\uparrow a_*) \setminus (\uparrow a)`, where `a_*` is the element
                covered by `a`. It is always a meet-irreducible element, if it exists.

                INPUT:

                - ``a`` -- a join-irreducible element of the lattice

                .. WARNING::

                    Element ``a`` is expected to be join-irreducible, and
                    this is *not* checked.

                OUTPUT:

                the element `\kappa(a)`

                This will raise a :exc:`ValueError` if there is not
                a unique greatest element with given constraints.

                EXAMPLES::

                    sage: V = ['b', 0, 1, 2, 3, 4, 't']
                    sage: C = [['b', 0], ['b', 1], [0, 2], [2, 3], [2, 4], [3, 't'], [1, 4], [4, 't']]
                    sage: L = LatticePoset([V, C], category=LatticePosets().Finite().Semidistributive())
                    sage: [(a, L.kappa(a)) for a in L.join_irreducibles()]
                    [(0, 1), (2, 0), (3, 4), (1, 3)]
                """
                H = self._hasse_diagram
                k = H.kappa(self._element_to_vertex(a))
                return self._vertex_to_element(k)

            def kappa_dual(self, a):
                r"""
                Return the minimum element smaller than the element covering
                ``a`` but not smaller than ``a``.

                Define `\kappa^*(a)` as the minimum element of
                `(\downarrow a_*) \setminus (\downarrow a)`, where `a_*` is the element
                covering `a`. It is always a join-irreducible element, if it exists.

                INPUT:

                - ``a`` -- a meet-irreducible element of the lattice

                .. WARNING::

                    Element ``a`` is expected to be meet-irreducible, and
                    this is *not* checked.

                OUTPUT:

                the element `\kappa^*(a)`

                This will raise a :exc:`ValueError` if there is not
                a unique greatest element with given constraints.

                EXAMPLES::

                    sage: V = ['b', 0, 1, 2, 3, 4, 't']
                    sage: C = [['b', 0], ['b', 1], [0, 2], [2, 3], [2, 4], [3, 't'], [1, 4], [4, 't']]
                    sage: L = LatticePoset([V, C], category=LatticePosets().Finite().Semidistributive())
                    sage: [(a, L.kappa_dual(a)) for a in L.meet_irreducibles()]
                    [(0, 2), (3, 1), (1, 0), (4, 3)]
                """
                H = self._hasse_diagram
                k = H.kappa_dual(self._element_to_vertex(a))
                return self._vertex_to_element(k)

            def rowmotion_semidistributive(self, a):
                r"""
                Return the image of the element ``a`` under
                semidistributive rowmotion in ``self``.

                Classical rowmotion is usually defined as an
                automorphism on the set of order ideals `J(P)` of a
                finite poset `P`.  It is a special case of
                semidistributive rowmotion because every distributive
                lattice is isomorphic to `J(P)` for some `P` by
                Birkhoff's representation theorem.

                .. SEEALSO::

                    If the image of rowmotion of several elements is
                    needed,
                    :class:`~sage.dynamics.finite_dynamical_system_catalog.semidistributive_rowmotion`
                    is much more efficient.

                EXAMPLES::

                    sage: V = ['b', 0, 1, 2, 3, 4, 't']
                    sage: C = [['b', 0], ['b', 1], [0, 2], [2, 3], [2, 4], [3, 't'], [1, 4], [4, 't']]
                    sage: L = LatticePoset([V, C], category=LatticePosets().Finite().Semidistributive())
                    sage: L.rowmotion_semidistributive(0)
                    2

                    sage: L = posets.TamariLattice(3)
                    sage: row = L.rowmotion_semidistributive
                    sage: DS = DiscreteDynamicalSystem(L, row)
                    sage: sorted([sorted([DyckWord(x[:-1]) for x in c]) for c in DS.cycles()])
                    [[[1, 0, 1, 0, 1, 0], [1, 1, 1, 0, 0, 0]],
                     [[1, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 0], [1, 1, 0, 1, 0, 0]]]
                    sage: L = posets.TamariLattice(4)
                    sage: L.rowmotion_semidistributive((1,1,0,1,1,0,0,0,0))
                    (1, 0, 1, 1, 0, 0, 1, 0, 0)

                Check that classical rowmotion is a special case of
                semidistributive rowmotion::

                    sage: T = posets.TamariLattice(3)
                    sage: L = T.order_ideals_lattice()
                    sage: all(L.rowmotion_semidistributive(a) == T.rowmotion(a) for a in L)
                    True

                    sage: P = posets.UpDownPoset(10)
                    sage: L = T.order_ideals_lattice()
                    sage: all(L.rowmotion_semidistributive(a) == T.rowmotion(a) for a in L)
                    True
                """
                kd = [self.kappa_dual(e) for e in self.canonical_meetands(a)]
                return self.join(kd)

            def spine(self):
                """
                Return the spine of ``self``.

                For a semidistributive lattice `L`, the *spine* of `L` is
                the distributive lattice constructed as the subposet
                on the union of longest maximal chains.

                EXAMPLES::

                    sage: P = posets.TamariLattice(4)
                    sage: S = P.spine(); S
                    Finite lattice containing 8 elements
                    sage: S.category()
                    Category of facade finite enumerated distributive lattices
                """
                from sage.combinat.posets.lattices import LatticePoset
                subset_H, _ = self._hasse_diagram.spine()
                subset = [self._vertex_to_element(v) for v in subset_H]

                H = self.hasse_diagram()
                covers = [(x, y) for x in subset for y in H.neighbors_in(x)
                          if y in subset]
                cat = LatticePosets().Finite().Distributive()
                return LatticePoset([subset, covers],
                                    cover_relations=True, category=cat)

    class CongruenceUniform(CategoryWithAxiom):
        """
        The category of congruence uniform lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().CongruenceUniform(); cat
            Category of finite congruence uniform lattice posets
            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of congruence uniform lattice posets]
        """
        @cached_method
        def extra_super_categories(self) -> list:
            r"""
            Return a list of the super categories of ``self``.

            These encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().CongruenceUniform().super_categories()
                [Category of finite lattice posets,
                 Category of congruence uniform lattice posets]
            """
            return [LatticePosets().Semidistributive()]

        class ParentMethods:
            def is_congruence_uniform(self) -> bool:
                """
                Return whether ``self`` is a congruence uniform lattice.

                EXAMPLES::

                    sage: posets.TamariLattice(4).is_congruence_uniform()
                    True
                """
                return True

    class Stone(CategoryWithAxiom):
        """
        The category of Stone lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().Stone(); cat
            Category of finite stone distributive lattices

            sage: cat.super_categories()
            [Category of finite distributive lattices,
             Category of stone lattice posets]
        """
        @cached_method
        def extra_super_categories(self) -> list:
            r"""
            Return a list of the super categories of ``self``.

            These encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().Stone().super_categories()
                [Category of finite distributive lattices,
                 Category of stone lattice posets]
            """
            return [LatticePosets().Trim().ChainGraded()]

        class ParentMethods:
            def is_stone(self) -> bool:
                """
                Return whether ``self`` is a Stone lattice.

                EXAMPLES::

                    sage: posets.DivisorLattice(12).is_stone()
                    True
                """
                return True

    class ChainGraded(CategoryWithAxiom):
        """
        The category of graded lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().ChainGraded(); cat
            Category of finite chain graded lattice posets

            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of chain graded lattice posets]
        """
        class ParentMethods:
            def is_graded(self) -> bool:
                """
                Return whether ``self`` is a graded lattice.

                EXAMPLES::

                    sage: posets.DivisorLattice(12).is_graded()
                    True
                """
                return True


# the following was moved out of the main class

class DistributiveLattices(CategoryWithAxiom):
    """
    The category of distributive lattices.

    EXAMPLES::

        sage: cat = FiniteLatticePosets().Distributive(); cat
        Category of finite distributive lattices

        sage: cat.super_categories()
        [Category of finite lattice posets,
         Category of distributive lattices]

    TESTS::

        sage: from sage.categories.lattice_posets import DistributiveLattices
        sage: LatticePosets().Distributive() is DistributiveLattices()
        True
    """
    _base_category_class_and_axiom = (LatticePosets.Trim,
                                      "ChainGraded")

    @cached_method
    def extra_super_categories(self) -> list:
        r"""
        Return a list of the super categories of ``self``.

        These encode implications between properties.

        EXAMPLES::

            sage: LatticePosets().Distributive().super_categories()
            [Category of congruence uniform lattice posets,
             Category of trim lattice posets,
             Category of chain graded lattice posets]
        """
        return [LatticePosets().CongruenceUniform()]

    class Finite(CategoryWithAxiom):
        pass

    class ParentMethods:
        def is_distributive(self) -> bool:
            """
            Return whether ``self`` is a distributive lattice.

            EXAMPLES::

                sage: P = posets.Crown(4).order_ideals_lattice()
                sage: P.is_distributive()
                True
            """
            return True


LatticePosets.Trim.ChainGraded = DistributiveLattices
