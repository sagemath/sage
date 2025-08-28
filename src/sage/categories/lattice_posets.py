# sage_setup: distribution = sagemath-categories
r"""
Lattice posets
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
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
    def super_categories(self):
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

            See :wikipedia:`Distributive lattice`.

            EXAMPLES::

                sage: P = posets.ChainPoset(2).order_ideals_lattice()
                sage: P in FiniteLatticePosets().Distributive()
                True
            """
            return self._with_axiom("Distributive")

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
            def is_extremal(self):
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
        def extra_super_categories(self):
            r"""
            Return a list of the super categories of ``self``.

            This encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().Trim().super_categories()
                [Category of finite lattice posets,
                 Category of trim lattice posets]
            """
            return [LatticePosets().Extremal()]

        class ParentMethods:
            def is_trim(self):
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
            def is_semidistributive(self):
                """
                Return whether ``self`` is a semidistributive lattice.

                EXAMPLES::

                    sage: posets.TamariLattice(4).is_semidistributive()
                    True
                """
                return True

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
        def extra_super_categories(self):
            r"""
            Return a list of the super categories of ``self``.

            This encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().CongruenceUniform().super_categories()
                [Category of finite lattice posets,
                 Category of congruence uniform lattice posets]
            """
            return [LatticePosets().Semidistributive()]

        class ParentMethods:
            def is_congruence_uniform(self):
                """
                Return whether ``self`` is a congruence uniform lattice.

                EXAMPLES::

                    sage: posets.TamariLattice(4).is_congruence_uniform()
                    True
                """
                return True

    class Distributive(CategoryWithAxiom):
        """
        The category of distributive lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().Distributive(); cat
            Category of finite distributive lattice posets

            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of distributive lattice posets]
        """
        @cached_method
        def extra_super_categories(self):
            r"""
            Return a list of the super categories of ``self``.

            This encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().Distributive().super_categories()
                [Category of finite lattice posets,
                 Category of distributive lattice posets]
            """
            return [LatticePosets().Trim(),
                    LatticePosets().CongruenceUniform()]

        class ParentMethods:
            def is_distributive(self):
                """
                Return whether ``self`` is a distributive lattice.

                EXAMPLES::

                    sage: P = posets.Crown(4).order_ideals_lattice()
                    sage: P.is_distributive()
                    True
                """
                return True

    class Stone(CategoryWithAxiom):
        """
        The category of Stone lattices.

        EXAMPLES::

            sage: cat = FiniteLatticePosets().Stone(); cat
            Category of finite stone lattice posets

            sage: cat.super_categories()
            [Category of finite lattice posets,
             Category of stone lattice posets]
        """
        @cached_method
        def extra_super_categories(self):
            r"""
            Return a list of the super categories of ``self``.

            This encode implications between properties.

            EXAMPLES::

                sage: FiniteLatticePosets().Stone().super_categories()
                [Category of finite lattice posets,
                 Category of stone lattice posets]
            """
            return [LatticePosets().Distributive()]

        class ParentMethods:
            def is_stone(self):
                """
                Return whether ``self`` is a Stone lattice.

                EXAMPLES::

                    sage: posets.DivisorLattice(12).is_stone()
                    True
                """
                return True
