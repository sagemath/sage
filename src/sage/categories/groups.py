# sage_setup: distribution = sagemath-categories
r"""
Groups
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.monoids import Monoids
from sage.categories.cartesian_product import CartesianProductsCategory, cartesian_product
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.topological_spaces import TopologicalSpacesCategory


class Groups(CategoryWithAxiom):
    """
    The category of (multiplicative) groups, i.e. monoids with
    inverses.

    EXAMPLES::

        sage: Groups()
        Category of groups
        sage: Groups().super_categories()
        [Category of monoids, Category of inverse unital magmas]

    TESTS::

        sage: TestSuite(Groups()).run()
    """
    _base_category_class_and_axiom = (Monoids, "Inverse")

    def example(self):
        """
        EXAMPLES::

            sage: Groups().example()                                                    # needs sage.modules
            General Linear Group of degree 4 over Rational Field
        """
        from sage.rings.rational_field import QQ
        from sage.groups.matrix_gps.linear import GL
        return GL(4,QQ)

    @staticmethod
    def free(index_set=None, names=None, **kwds):
        r"""
        Return the free group.

        INPUT:

        - ``index_set`` -- (optional) an index set for the generators; if
          an integer, then this represents `\{0, 1, \ldots, n-1\}`

        - ``names`` -- string or list/tuple/iterable of strings
          (default: ``'x'``); the generator names or name prefix

        When the index set is an integer or only variable names are given,
        this returns :class:`~sage.groups.free_group.FreeGroup_class`, which
        currently has more features due to the interface with GAP than
        :class:`~sage.groups.indexed_free_group.IndexedFreeGroup`.

        EXAMPLES::

            sage: # needs sage.combinat sage.groups
            sage: Groups.free(index_set=ZZ)
            Free group indexed by Integer Ring
            sage: Groups().free(ZZ)
            Free group indexed by Integer Ring
            sage: Groups().free(5)
            Free Group on generators {x0, x1, x2, x3, x4}
            sage: F.<x,y,z> = Groups().free(); F
            Free Group on generators {x, y, z}
        """
        from sage.rings.integer_ring import ZZ
        if index_set in ZZ or (index_set is None and names is not None):
            from sage.groups.free_group import FreeGroup
            if names is None:
                return FreeGroup(index_set, **kwds)
            return FreeGroup(index_set, names, **kwds)

        from sage.groups.indexed_free_group import IndexedFreeGroup
        return IndexedFreeGroup(index_set, **kwds)

    class ParentMethods:

        def group_generators(self):
            """
            Return group generators for ``self``.

            This default implementation calls :meth:`gens`, for
            backward compatibility.

            EXAMPLES::

                sage: A = AlternatingGroup(4)                                           # needs sage.groups
                sage: A.group_generators()                                              # needs sage.groups
                Family ((1,2,3), (2,3,4))
            """
            from sage.sets.family import Family
            try:
                return Family(self.gens())
            except AttributeError:
                raise NotImplementedError("no generators are implemented for this group")

        def monoid_generators(self):
            r"""
            Return the generators of ``self`` as a monoid.

            Let `G` be a group with generating set `X`. In general, the
            generating set of `G` as a monoid is given by `X \cup X^{-1}`,
            where `X^{-1}` is the set of inverses of `X`. If `G` is a finite
            group, then the generating set as a monoid is `X`.

            EXAMPLES::

                sage: # needs sage.groups
                sage: A = AlternatingGroup(4)
                sage: A.monoid_generators()
                Family ((1,2,3), (2,3,4))
                sage: F.<x,y> = FreeGroup()
                sage: F.monoid_generators()
                Family (x, y, x^-1, y^-1)
            """
            G = self.group_generators()
            from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
            if G not in FiniteEnumeratedSets():
                raise NotImplementedError("currently only implemented for finitely generated groups")
            from sage.sets.family import Family
            return Family(tuple(G) + tuple(~x for x in G))

        def _test_inverse(self, **options):
            """
            Run generic tests on the method :meth:`.__invert__`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: G = SymmetricGroup(3)                                             # needs sage.groups
                sage: G._test_inverse()                                                 # needs sage.groups
            """
            tester = self._tester(**options)
            for x in tester.some_elements():
                tester.assertEqual(x * ~x, self.one())
                tester.assertEqual(~x * x, self.one())

        def semidirect_product(self, N, mapping, check=True):
            r"""
            The semi-direct product of two groups.

            EXAMPLES::

                sage: G = Groups().example()                                            # needs sage.modules
                sage: G.semidirect_product(G, Morphism(G, G))                           # needs sage.modules
                Traceback (most recent call last):
                ...
                NotImplementedError: semidirect product of General Linear Group of degree 4
                over Rational Field and General Linear Group of degree 4 over Rational Field
                not yet implemented
            """
            raise NotImplementedError("semidirect product of %s and %s not yet implemented" % (self, N))

        def holomorph(self):
            r"""
            The holomorph of a group.

            The holomorph of a group `G` is the semidirect product
            `G \rtimes_{id} Aut(G)`, where `id` is the identity function
            on `Aut(G)`, the automorphism group of `G`.

            See :wikipedia:`Holomorph (mathematics)`

            EXAMPLES::

                sage: G = Groups().example()                                            # needs sage.modules
                sage: G.holomorph()                                                     # needs sage.modules
                Traceback (most recent call last):
                ...
                NotImplementedError: holomorph of General Linear Group of degree 4
                over Rational Field not yet implemented
            """
            raise NotImplementedError("holomorph of %s not yet implemented" % self)

        def cayley_table(self, names='letters', elements=None):
            r"""
            Return the "multiplication" table of this multiplicative group,
            which is also known as the "Cayley table".

            .. NOTE:: The order of the elements in the row and column
              headings is equal to the order given by the table's
              :meth:`~sage.matrix.operation_table.OperationTable.column_keys`
              method.  The association between the actual elements and the
              names/symbols used in the table can also be retrieved as
              a dictionary with the
              :meth:`~sage.matrix.operation_table.OperationTable.translation`
              method.

            For groups, this routine should behave identically to the
            :meth:`~sage.categories.magmas.Magmas.ParentMethods.multiplication_table`
            method for magmas, which applies in greater generality.

            INPUT:

            - ``names`` -- the type of names used, values are:

              * ``'letters'`` -- lowercase ASCII letters are used
                for a base 26 representation of the elements'
                positions in the list given by :meth:`list`,
                padded to a common width with leading 'a's.
              * ``'digits'`` -- base 10 representation of the
                elements' positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading zeros.
              * ``'elements'`` -- the string representations
                of the elements themselves.
              * a list - a list of strings, where the length
                of the list equals the number of elements.

            - ``elements`` -- (default: ``None``) a list of
              elements of the group, in forms that can be
              coerced into the structure, eg. their string
              representations. This may be used to impose an
              alternate ordering on the elements, perhaps when
              this is used in the context of a particular structure.
              The default is to use whatever ordering is provided by the
              the group, which is reported by the
              :meth:`~sage.matrix.operation_table.OperationTable.column_keys`
              method.  Or the ``elements`` can be a subset
              which is closed under the operation. In particular,
              this can be used when the base set is infinite.

            OUTPUT:

            An object representing the multiplication table.  This is
            an :class:`~sage.matrix.operation_table.OperationTable` object
            and even more documentation can be found there.

            EXAMPLES:

            Permutation groups, matrix groups and abelian groups
            can all compute their multiplication tables.  ::

                sage: # needs sage.groups
                sage: G = DiCyclicGroup(3)
                sage: T = G.cayley_table()
                sage: T.column_keys()
                ((), (5,6,7), ..., (1,4,2,3)(5,7))
                sage: T
                *  a b c d e f g h i j k l
                 +------------------------
                a| a b c d e f g h i j k l
                b| b c a e f d i g h l j k
                c| c a b f d e h i g k l j
                d| d e f a b c j k l g h i
                e| e f d b c a l j k i g h
                f| f d e c a b k l j h i g
                g| g h i j k l d e f a b c
                h| h i g k l j f d e c a b
                i| i g h l j k e f d b c a
                j| j k l g h i a b c d e f
                k| k l j h i g c a b f d e
                l| l j k i g h b c a e f d

            ::

                sage: M = SL(2, 2)                                                      # needs sage.modules
                sage: M.cayley_table()                                                  # needs sage.libs.gap sage.modules
                *  a b c d e f
                 +------------
                a| a b c d e f
                b| b a d c f e
                c| c e a f b d
                d| d f b e a c
                e| e c f a d b
                f| f d e b c a
                <BLANKLINE>

            ::

                sage: A = AbelianGroup([2, 3])                                          # needs sage.groups
                sage: A.cayley_table()                                                  # needs sage.groups
                *  a b c d e f
                 +------------
                a| a b c d e f
                b| b c a e f d
                c| c a b f d e
                d| d e f a b c
                e| e f d b c a
                f| f d e c a b

            Lowercase ASCII letters are the default symbols used
            for the table, but you can also specify the use of
            decimal digit strings, or provide your own strings
            (in the proper order if they have meaning).
            Also, if the elements themselves are not too complex,
            you can choose to just use the string representations
            of the elements themselves.  ::

                sage: C = CyclicPermutationGroup(11)                                    # needs sage.groups
                sage: C.cayley_table(names='digits')                                    # needs sage.groups
                 *  00 01 02 03 04 05 06 07 08 09 10
                  +---------------------------------
                00| 00 01 02 03 04 05 06 07 08 09 10
                01| 01 02 03 04 05 06 07 08 09 10 00
                02| 02 03 04 05 06 07 08 09 10 00 01
                03| 03 04 05 06 07 08 09 10 00 01 02
                04| 04 05 06 07 08 09 10 00 01 02 03
                05| 05 06 07 08 09 10 00 01 02 03 04
                06| 06 07 08 09 10 00 01 02 03 04 05
                07| 07 08 09 10 00 01 02 03 04 05 06
                08| 08 09 10 00 01 02 03 04 05 06 07
                09| 09 10 00 01 02 03 04 05 06 07 08
                10| 10 00 01 02 03 04 05 06 07 08 09

            ::

                sage: G = QuaternionGroup()                                             # needs sage.groups
                sage: names = ['1', 'I', '-1', '-I', 'J', '-K', '-J', 'K']
                sage: G.cayley_table(names=names)                                       # needs sage.groups
                 *   1  I -1 -I  J -K -J  K
                  +------------------------
                 1|  1  I -1 -I  J -K -J  K
                 I|  I -1 -I  1  K  J -K -J
                -1| -1 -I  1  I -J  K  J -K
                -I| -I  1  I -1 -K -J  K  J
                 J|  J -K -J  K -1 -I  1  I
                -K| -K -J  K  J  I -1 -I  1
                -J| -J  K  J -K  1  I -1 -I
                 K|  K  J -K -J -I  1  I -1

            ::

                sage: A = AbelianGroup([2, 2])                                          # needs sage.groups
                sage: A.cayley_table(names='elements')                                  # needs sage.groups
                    *      1    f1    f0 f0*f1
                     +------------------------
                    1|     1    f1    f0 f0*f1
                   f1|    f1     1 f0*f1    f0
                   f0|    f0 f0*f1     1    f1
                f0*f1| f0*f1    f0    f1     1

            The :meth:`~sage.matrix.operation_table.OperationTable.change_names`
            routine behaves similarly, but changes an existing table "in-place."
            ::

                sage: # needs sage.groups
                sage: G = AlternatingGroup(3)
                sage: T = G.cayley_table()
                sage: T.change_names('digits')
                sage: T
                *  0 1 2
                 +------
                0| 0 1 2
                1| 1 2 0
                2| 2 0 1

            For an infinite group, you can still work with finite sets of
            elements, provided the set is closed under multiplication.
            Elements will be coerced into the group as part of setting
            up the table.  ::

                sage: # needs sage.modules
                sage: G = SL(2,ZZ); G
                Special Linear Group of degree 2 over Integer Ring
                sage: identity = matrix(ZZ, [[1,0], [0,1]])
                sage: G.cayley_table(elements=[identity, -identity])
                *  a b
                 +----
                a| a b
                b| b a

            The
            :class:`~sage.matrix.operation_table.OperationTable`
            class provides even greater flexibility, including changing
            the operation.  Here is one such example, illustrating the
            computation of commutators.  ``commutator`` is defined as
            a function of two variables, before being used to build
            the table. From this, the commutator subgroup seems obvious,
            and creating a Cayley table with just these three elements
            confirms that they form a closed subset in the group.
            ::

                sage: # needs sage.groups sage.modules
                sage: from sage.matrix.operation_table import OperationTable
                sage: G = DiCyclicGroup(3)
                sage: commutator = lambda x, y: x*y*x^-1*y^-1
                sage: T = OperationTable(G, commutator); T
                .  a b c d e f g h i j k l
                 +------------------------
                a| a a a a a a a a a a a a
                b| a a a a a a c c c c c c
                c| a a a a a a b b b b b b
                d| a a a a a a a a a a a a
                e| a a a a a a c c c c c c
                f| a a a a a a b b b b b b
                g| a b c a b c a c b a c b
                h| a b c a b c b a c b a c
                i| a b c a b c c b a c b a
                j| a b c a b c a c b a c b
                k| a b c a b c b a c b a c
                l| a b c a b c c b a c b a
                sage: trans = T.translation()
                sage: comm = [trans['a'], trans['b'], trans['c']]
                sage: comm
                [(), (5,6,7), (5,7,6)]
                sage: P = G.cayley_table(elements=comm)
                sage: P
                *  a b c
                 +------
                a| a b c
                b| b c a
                c| c a b

            .. TODO::

                Arrange an ordering of elements into cosets of a normal
                subgroup close to size `\sqrt{n}`.  Then the quotient
                group structure is often apparent in the table.  See
                comments on :issue:`7555`.

            AUTHOR:

            - Rob Beezer (2010-03-15)
            """
            from sage.matrix.operation_table import OperationTable
            import operator
            return OperationTable(self, operation=operator.mul, names=names, elements=elements)

        def conjugacy_class(self, g):
            r"""
            Return the conjugacy class of the element ``g``.

            This is a fall-back method for groups not defined over GAP.

            EXAMPLES::

                sage: A = AbelianGroup([2, 2])                                          # needs sage.groups
                sage: c = A.conjugacy_class(A.an_element())                             # needs sage.groups
                sage: type(c)                                                           # needs sage.groups
                <class 'sage.groups.conjugacy_classes.ConjugacyClass_with_category'>
            """
            from sage.groups.conjugacy_classes import ConjugacyClass
            return ConjugacyClass(self, g)

        def minimum_generating_set(self):
            """
            Return a list of the minimum generating set of this group.

            EXAMPLES::

                sage: G = GL(2,GF(3))
                sage: G.minimum_generating_set()
                [
                [1 2]  [1 2]
                [0 1], [1 1]
                ]

                sage: G = SymmetricGroup(3)
                sage: s = G.minimum_generating_set(); s
                [(2,3), (1,2,3)]
                sage: s[0].parent()
                Symmetric group of order 3! as a permutation group

                sage: A5 = AlternatingGroup(5)
                sage: A5.minimum_generating_set()
                [(1,2,3,4,5), (3,4,5)]

                sage: H = groups.matrix.Heisenberg(1,3); H
                Heisenberg group of degree 1 over Ring of integers modulo 3
                sage: H.minimum_generating_set()
                [
                [1 1 0]  [1 0 0]
                [0 1 0]  [0 1 1]
                [0 0 1], [0 0 1]
                ]

            TESTS:

            Test that function gives an error for infinite groups::

                sage: G = GL(2, ZZ)
                sage: G.minimum_generating_set()
                Traceback (most recent call last):
                ...
                NotImplementedError: only implemented for finite groups
            """
            from sage.groups.libgap_mixin import minimum_generating_set
            return [self(x) for x in minimum_generating_set(self)]

    class ElementMethods:
        def conjugacy_class(self):
            r"""
            Return the conjugacy class of ``self``.

            EXAMPLES::

                sage: D = DihedralGroup(5)                                              # needs sage.groups
                sage: g = D((1,3,5,2,4))                                                # needs sage.groups
                sage: g.conjugacy_class()                                               # needs sage.groups
                Conjugacy class of (1,3,5,2,4)
                 in Dihedral group of order 10 as a permutation group

                sage: H = MatrixGroup([matrix(GF(5), 2, [1,2, -1,1]),                   # needs sage.modules
                ....:                  matrix(GF(5), 2, [1,1, 0,1])])
                sage: h = H(matrix(GF(5), 2, [1,2, -1,1]))                              # needs sage.modules
                sage: h.conjugacy_class()                                               # needs sage.groups sage.modules
                Conjugacy class of [1 2]
                                   [4 1]
                in Matrix group over Finite Field of size 5 with 2 generators (
                [1 2]  [1 1]
                [4 1], [0 1]
                )

                sage: G = SL(2, GF(2))                                                  # needs sage.modules
                sage: g = G.gens()[0]                                                   # needs sage.groups sage.modules
                sage: g.conjugacy_class()                                               # needs sage.groups sage.modules
                Conjugacy class of [1 1]
                [0 1] in Special Linear Group of degree 2 over Finite Field of size 2

                sage: G = SL(2, QQ)                                                     # needs sage.modules
                sage: g = G([[1,1], [0,1]])                                             # needs sage.modules
                sage: g.conjugacy_class()                                               # needs sage.groups sage.modules
                Conjugacy class of [1 1]
                [0 1] in Special Linear Group of degree 2 over Rational Field
            """
            return self.parent().conjugacy_class(self)

    Finite = LazyImport('sage.categories.finite_groups', 'FiniteGroups', at_startup=True)
    Lie = LazyImport('sage.categories.lie_groups', 'LieGroups', 'Lie')
    Algebras = LazyImport('sage.categories.group_algebras', 'GroupAlgebras', at_startup=True)

    class Commutative(CategoryWithAxiom):
        r"""
        Category of commutative (abelian) groups.

        A group `G` is *commutative* if `xy = yx` for all `x,y \in G`.
        """
        @staticmethod
        def free(index_set=None, names=None, **kwds):
            r"""
            Return the free commutative group.

            INPUT:

            - ``index_set`` -- (optional) an index set for the generators; if
              an integer, then this represents `\{0, 1, \ldots, n-1\}`

            - ``names`` -- string or list/tuple/iterable of strings
              (default: ``'x'``); the generator names or name prefix

            EXAMPLES::

                sage: # needs sage.combinat sage.groups
                sage: Groups.Commutative.free(index_set=ZZ)
                Free abelian group indexed by Integer Ring
                sage: Groups().Commutative().free(ZZ)
                Free abelian group indexed by Integer Ring
                sage: Groups().Commutative().free(5)
                Multiplicative Abelian group isomorphic to Z x Z x Z x Z x Z
                sage: F.<x,y,z> = Groups().Commutative().free(); F
                Multiplicative Abelian group isomorphic to Z x Z x Z
            """
            from sage.rings.integer_ring import ZZ
            if names is not None:
                if isinstance(names, str):
                    if ',' not in names and index_set in ZZ:
                        names = [names + repr(i) for i in range(index_set)]
                    else:
                        names = names.split(',')
                names = tuple(names)
                if index_set is None:
                    index_set = ZZ(len(names))
                if index_set in ZZ:
                    from sage.groups.abelian_gps.abelian_group import AbelianGroup
                    return AbelianGroup(index_set, names=names, **kwds)

            if index_set in ZZ:
                from sage.groups.abelian_gps.abelian_group import AbelianGroup
                return AbelianGroup(index_set, **kwds)

            from sage.groups.indexed_free_group import IndexedFreeAbelianGroup
            return IndexedFreeAbelianGroup(index_set, names=names, **kwds)

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of groups constructed as Cartesian products of groups.

        This construction gives the direct product of groups. See
        :wikipedia:`Direct_product` and :wikipedia:`Direct_product_of_groups`
        for more information.
        """
        def extra_super_categories(self):
            """
            A Cartesian product of groups is endowed with a natural
            group structure.

            EXAMPLES::

                sage: C = Groups().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of groups]
                sage: sorted(C.super_categories(), key=str)
                [Category of Cartesian products of inverse unital magmas,
                 Category of Cartesian products of monoids,
                 Category of groups]
            """
            return [self.base_category()]

        class ParentMethods:
            @cached_method
            def group_generators(self):
                """
                Return the group generators of ``self``.

                EXAMPLES::

                    sage: # needs sage.groups
                    sage: C5 = CyclicPermutationGroup(5)
                    sage: C4 = CyclicPermutationGroup(4)
                    sage: S4 = SymmetricGroup(3)
                    sage: C = cartesian_product([C5, C4, S4])
                    sage: C.group_generators()
                    Family (((1,2,3,4,5), (), ()),
                            ((), (1,2,3,4), ()),
                            ((), (), (1,2)),
                            ((), (), (2,3)))

                We check the other portion of :issue:`16718` is fixed::

                    sage: len(C.j_classes())                                            # needs sage.graphs sage.groups
                    1

                An example with an infinitely generated group (a better output
                is needed)::

                    sage: # needs sage.combinat sage.groups
                    sage: G = Groups.free([1,2])
                    sage: H = Groups.free(ZZ)
                    sage: C = cartesian_product([G, H])
                    sage: C.monoid_generators()
                    Lazy family (gen(i))_{i in The Cartesian product of (...)}
                """
                F = self.cartesian_factors()
                ids = tuple(G.one() for G in F)

                def lift(i, gen):
                    cur = list(ids)
                    cur[i] = gen
                    return self._cartesian_product_of_elements(cur)
                from sage.sets.family import Family

                # Finitely generated
                cat = FiniteEnumeratedSets()
                if all(G.group_generators() in cat
                       or isinstance(G.group_generators(), (tuple, list)) for G in F):
                    ret = [lift(i, gen) for i, G in enumerate(F) for gen in G.group_generators()]
                    return Family(ret)

                # Infinitely generated
                # This does not return a good output, but it is "correct"
                # TODO: Figure out a better way to do things
                gens_prod = cartesian_product([Family(G.group_generators(),
                                                      lambda g: (i, g))
                                               for i, G in enumerate(F)])
                return Family(gens_prod, lift, name='gen')

            def order(self):
                r"""
                Return the cardinality of ``self``.

                EXAMPLES::

                    sage: C = cartesian_product([SymmetricGroup(10), SL(2, GF(3))])     # needs sage.groups sage.modules
                    sage: C.order()                                                     # needs sage.groups sage.modules
                    87091200

                TESTS::

                    sage: C.order.__module__                                            # needs sage.groups sage.modules
                    'sage.categories.groups'

                .. TODO::

                    this method is just here to prevent
                    ``FiniteGroups.ParentMethods`` to call
                    ``_cardinality_from_iterator``.
                """
                from sage.misc.misc_c import prod
                return prod(c.cardinality() for c in self.cartesian_factors())

    class Topological(TopologicalSpacesCategory):
        """
        Category of topological groups.

        A topological group `G` is a group which has a topology such that
        multiplication and taking inverses are continuous functions.

        REFERENCES:

        - :wikipedia:`Topological_group`
        """
