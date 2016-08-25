r"""
Groups
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.monoids import Monoids
from sage.categories.algebra_functor import AlgebrasCategory
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

            sage: Groups().example()
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

        - ``names`` -- a string or list/tuple/iterable of strings
          (default: ``'x'``); the generator names or name prefix

        When the index set is an integer or only variable names are given,
        this returns :class:`~sage.groups.free_group.FreeGroup_class`, which
        currently has more features due to the interface with GAP than
        :class:`~sage.groups.indexed_free_group.IndexedFreeGroup`.

        EXAMPLES::

            sage: Groups.free(index_set=ZZ)
            Free group indexed by Integer Ring
            sage: Groups().free(ZZ)
            Free group indexed by Integer Ring
            sage: Groups().free(5)
            Free Group on generators {x0, x1, x2, x3, x4}
            sage: F.<x,y,z> = Groups().free(); F
            Free Group on generators {x, y, z}
        """
        from sage.rings.all import ZZ
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
            Returns group generators for self.

            This default implementation calls :meth:`.gens`, for
            backward compatibility.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.group_generators()
                Family ((2,3,4), (1,2,3))
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

                sage: A = AlternatingGroup(4)
                sage: A.monoid_generators()
                Family ((2,3,4), (1,2,3))
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

                sage: G = SymmetricGroup(3)
                sage: G._test_inverse()
            """
            tester = self._tester(**options)
            for x in tester.some_elements():
                tester.assertEquals(x * ~x, self.one())
                tester.assertEquals(~x * x, self.one())

        def semidirect_product(self, N, mapping, check = True):
            r"""
            The semi-direct product of two groups

            EXAMPLES::

                sage: G = Groups().example()
                sage: G.semidirect_product(G,Morphism(G,G))
                Traceback (most recent call last):
                ...
                NotImplementedError: semidirect product of General Linear Group of degree 4 over Rational Field and General Linear Group of degree 4 over Rational Field not yet implemented
            """
            raise NotImplementedError("semidirect product of %s and %s not yet implemented"%(self, N))

        def holomorph(self):
            r"""
            The holomorph of a group

            The holomorph of a group `G` is the semidirect product
            `G \rtimes_{id} Aut(G)`, where `id` is the identity function
            on `Aut(G)`, the automorphism group of `G`.

            See :wikipedia:`Holomorph (mathematics)`

            EXAMPLES::

                sage: G = Groups().example()
                sage: G.holomorph()
                Traceback (most recent call last):
                ...
                NotImplementedError: holomorph of General Linear Group of degree 4 over Rational Field not yet implemented
            """
            raise NotImplementedError("holomorph of %s not yet implemented"%self)

        def cayley_table(self, names='letters', elements=None):
            r"""
            Returns the "multiplication" table of this multiplicative group,
            which is also known as the "Cayley table".

            .. note:: The order of the elements in the row and column
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

            - ``names`` - the type of names used, values are:

              * ``'letters'`` - lowercase ASCII letters are used
                for a base 26 representation of the elements'
                positions in the list given by :meth:`list`,
                padded to a common width with leading 'a's.
              * ``'digits'`` - base 10 representation of the
                elements' positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading zeros.
              * ``'elements'`` - the string representations
                of the elements themselves.
              * a list - a list of strings, where the length
                of the list equals the number of elements.

            - ``elements`` - default = ``None``.  A list of
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

                sage: G = DiCyclicGroup(3)
                sage: T = G.cayley_table()
                sage: T.column_keys()
                ((), (1,3,2,4)(5,7), ..., (1,2)(3,4)(5,7,6))
                sage: T
                *  a b c d e f g h i j k l
                 +------------------------
                a| a b c d e f g h i j k l
                b| b e f j i h d k a l c g
                c| c g d e h b k l j f i a
                d| d k e h l g i a f b j c
                e| e i h l a k j c b g f d
                f| f d j i k e c g l h a b
                g| g h b f j l e i c a d k
                h| h j l a c i f d g k b e
                i| i a k g b c l f e d h j
                j| j c i k g d a b h e l f
                k| k l g b f a h j d c e i
                l| l f a c d j b e k i g h

            ::

                sage: M=SL(2,2)
                sage: M.cayley_table()
                *  a b c d e f
                 +------------
                a| a b c d e f
                b| b a d c f e
                c| c f e b a d
                d| d e f a b c
                e| e d a f c b
                f| f c b e d a
                <BLANKLINE>

            ::

                sage: A=AbelianGroup([2,3])
                sage: A.cayley_table()
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

                sage: C=CyclicPermutationGroup(11)
                sage: C.cayley_table(names='digits')
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

                sage: G=QuaternionGroup()
                sage: names=['1', 'I', 'J', '-1', '-K', 'K', '-I', '-J']
                sage: G.cayley_table(names=names)
                 *   1  I  J -1 -K  K -I -J
                  +------------------------
                 1|  1  I  J -1 -K  K -I -J
                 I|  I -1  K -I  J -J  1 -K
                 J|  J -K -1 -J -I  I  K  1
                -1| -1 -I -J  1  K -K  I  J
                -K| -K -J  I  K -1  1  J -I
                 K|  K  J -I -K  1 -1 -J  I
                -I| -I  1 -K  I -J  J -1  K
                -J| -J  K  1  J  I -I -K -1

            ::

                sage: A=AbelianGroup([2,2])
                sage: A.cayley_table(names='elements')
                    *      1    f1    f0 f0*f1
                     +------------------------
                    1|     1    f1    f0 f0*f1
                   f1|    f1     1 f0*f1    f0
                   f0|    f0 f0*f1     1    f1
                f0*f1| f0*f1    f0    f1     1

            The :meth:`~sage.matrix.operation_table.OperationTable.change_names`
            routine behaves similarly, but changes an existing table "in-place."
            ::

                sage: G=AlternatingGroup(3)
                sage: T=G.cayley_table()
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

                sage: G=SL(2,ZZ)
                sage: G
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

                sage: from sage.matrix.operation_table import OperationTable
                sage: G=DiCyclicGroup(3)
                sage: commutator = lambda x, y: x*y*x^-1*y^-1
                sage: T=OperationTable(G, commutator)
                sage: T
                .  a b c d e f g h i j k l
                 +------------------------
                a| a a a a a a a a a a a a
                b| a a h d a d h h a h d d
                c| a d a a a d d a d d d a
                d| a h a a a h h a h h h a
                e| a a a a a a a a a a a a
                f| a h h d a a d h h d a d
                g| a d h d a h a h d a h d
                h| a d a a a d d a d d d a
                i| a a h d a d h h a h d d
                j| a d h d a h a h d a h d
                k| a h h d a a d h h d a d
                l| a h a a a h h a h h h a
                sage: trans = T.translation()
                sage: comm = [trans['a'], trans['d'],trans['h']]
                sage: comm
                [(), (5,7,6), (5,6,7)]
                sage: P=G.cayley_table(elements=comm)
                sage: P
                *  a b c
                 +------
                a| a b c
                b| b c a
                c| c a b

            TODO:

            Arrange an ordering of elements into cosets of a normal
            subgroup close to size `\sqrt{n}`.  Then the quotient
            group structure is often apparent in the table.  See
            comments on :trac:`7555`.

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

                sage: A = AbelianGroup([2,2])
                sage: c = A.conjugacy_class(A.an_element())
                sage: type(c)
                <class 'sage.groups.conjugacy_classes.ConjugacyClass_with_category'>
            """
            from sage.groups.conjugacy_classes import ConjugacyClass
            return ConjugacyClass(self, g)

    class ElementMethods:
        def conjugacy_class(self):
            r"""
            Return the conjugacy class of ``self``.

            EXAMPLES::

                sage: D = DihedralGroup(5)
                sage: g = D((1,3,5,2,4))
                sage: g.conjugacy_class()
                Conjugacy class of (1,3,5,2,4) in Dihedral group of order 10 as a permutation group

                sage: H = MatrixGroup([matrix(GF(5),2,[1,2, -1, 1]), matrix(GF(5),2, [1,1, 0,1])])
                sage: h = H(matrix(GF(5),2,[1,2, -1, 1]))
                sage: h.conjugacy_class()
                Conjugacy class of [1 2]
                [4 1] in Matrix group over Finite Field of size 5 with 2 generators (
                [1 2]  [1 1]
                [4 1], [0 1]
                )

                sage: G = SL(2, GF(2))
                sage: g = G.gens()[0]
                sage: g.conjugacy_class()
                Conjugacy class of [1 1]
                [0 1] in Special Linear Group of degree 2 over Finite Field of size 2

                sage: G = SL(2, QQ)
                sage: g = G([[1,1],[0,1]])
                sage: g.conjugacy_class()
                Conjugacy class of [1 1]
                [0 1] in Special Linear Group of degree 2 over Rational Field
            """
            return self.parent().conjugacy_class(self)

    Finite = LazyImport('sage.categories.finite_groups', 'FiniteGroups')
    Lie = LazyImport('sage.categories.lie_groups', 'LieGroups', 'Lie')
    #Algebras = LazyImport('sage.categories.group_algebras', 'GroupAlgebras')

    class Commutative(CategoryWithAxiom):
        """
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

            - ``names`` -- a string or list/tuple/iterable of strings
              (default: ``'x'``); the generator names or name prefix

            EXAMPLES::

                sage: Groups.Commutative.free(index_set=ZZ)
                Free abelian group indexed by Integer Ring
                sage: Groups().Commutative().free(ZZ)
                Free abelian group indexed by Integer Ring
                sage: Groups().Commutative().free(5)
                Multiplicative Abelian group isomorphic to Z x Z x Z x Z x Z
                sage: F.<x,y,z> = Groups().Commutative().free(); F
                Multiplicative Abelian group isomorphic to Z x Z x Z
            """
            from sage.rings.all import ZZ
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

    class Algebras(AlgebrasCategory):
        r"""
        The category of group algebras over a given base ring.

        EXAMPLES::

            sage: GroupAlgebras(IntegerRing())
            Category of group algebras over Integer Ring
            sage: GroupAlgebras(IntegerRing()).super_categories()
            [Category of hopf algebras with basis over Integer Ring,
             Category of monoid algebras over Integer Ring]

        Here is how to create the group algebra of a group `G`::

            sage: G = DihedralGroup(5)
            sage: QG = G.algebra(QQ); QG
            Group algebra of Dihedral group of order 10 as a permutation group over Rational Field

        and an example of computation::

            sage: g = G.an_element(); g
            (1,2,3,4,5)
            sage: (QG.term(g) + 1)**3
            B[()] + 3*B[(1,2,3,4,5)] + 3*B[(1,3,5,2,4)] + B[(1,4,2,5,3)]

        .. TODO::

            - Check which methods would be better located in
              ``Monoid.Algebras`` or ``Groups.Finite.Algebras``.

        TESTS::

            sage: A = GroupAlgebras(QQ).example(GL(3, GF(11)))
            sage: A.one_basis()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: A = SymmetricGroupAlgebra(QQ,4)
            sage: x = Permutation([4,3,2,1])
            sage: A.product_on_basis(x,x)
            [1, 2, 3, 4]

            sage: C = GroupAlgebras(ZZ)
            sage: TestSuite(C).run()
        """

        def extra_super_categories(self):
            """
            Implement the fact that the algebra of a group is a Hopf
            algebra.

            EXAMPLES::

                sage: C = Groups().Algebras(QQ)
                sage: C.extra_super_categories()
                [Category of hopf algebras over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of hopf algebras with basis over Rational Field,
                 Category of monoid algebras over Rational Field]
            """
            from sage.categories.hopf_algebras import HopfAlgebras
            return [HopfAlgebras(self.base_ring())]

        def example(self, G = None):
            """
            Return an example of group algebra.

            EXAMPLES::

                sage: GroupAlgebras(QQ['x']).example()
                Group algebra of Dihedral group of order 8 as a permutation group over Univariate Polynomial Ring in x over Rational Field

            An other group can be specified as optional argument::

                sage: GroupAlgebras(QQ).example(AlternatingGroup(4))
                Group algebra of Alternating group of order 4!/2 as a permutation group over Rational Field
            """
            from sage.groups.perm_gps.permgroup_named import DihedralGroup
            if G is None:
                G = DihedralGroup(4)
            return G.algebra(self.base_ring())

        class ParentMethods:

            def _repr_(self):
                r"""
                Return the string representation of `self`.

                EXAMPLES::

                    sage: A = Groups().example().algebra(QQ); A
                    Group algebra of General Linear Group of degree 4 over Rational Field over Rational Field
                    sage: A._name= "foo"
                    sage: A
                    foo over Rational Field
                """
                if hasattr(self, "_name"):
                    return self._name + " over {}".format(self.base_ring())
                else:
                    return 'Group algebra of {} over {}'.format(self.basis().keys(),
                                                                self.base_ring())

            def group(self):
                r"""
                Return the underlying group of the group algebra.

                EXAMPLES::

                    sage: GroupAlgebras(QQ).example(GL(3, GF(11))).group()
                    General Linear Group of degree 3 over Finite Field of size 11
                    sage: SymmetricGroup(10).algebra(QQ).group()
                    Symmetric group of order 10! as a permutation group
                """
                return self.basis().keys()

            def algebra_generators(self):
                r"""
                Return generators of this group algebra (as an algebra).

                EXAMPLES::

                    sage: GroupAlgebras(QQ).example(AlternatingGroup(10)).algebra_generators()
                    Finite family {(8,9,10): B[(8,9,10)], (1,2,3,4,5,6,7,8,9): B[(1,2,3,4,5,6,7,8,9)]}
                """
                from sage.sets.family import Family
                return Family(self.group().gens(), self.term)

            def center_basis(self):
                r"""
                Return a basis of the center of the group algebra.

                The canonical basis of the center of the group algebra
                is the family `(f_\sigma)_{\sigma\in C}`, where `C` is
                any collection of representatives of the conjugacy
                classes of the group, and `f_\sigma` is the sum of the
                elements in the conjugacy class of `\sigma`.

                OUTPUT:

                - ``list`` of elements of ``self``

                .. WARNING::

                    - This method requires the underlying group to
                      have a method ``conjugacy_classes``
                      (every permutation group has one, thanks GAP!).

                EXAMPLES::

                    sage: SymmetricGroup(3).algebra(QQ).center_basis()
                    [(), (2,3) + (1,2) + (1,3), (1,2,3) + (1,3,2)]

                .. SEEALSO::

                    - :meth:`Groups.Algebras.ElementMethods.central_form`
                    - :meth:`Monoids.Algebras.ElementMethods.is_central`
                """
                return [self.sum_of_monomials(conj) for conj  in
                        self.basis().keys().conjugacy_classes()]

            # Coalgebra structure

            def coproduct_on_basis(self, g):
                r"""
                Return the coproduct of the element ``g`` of the basis.

                Each basis element ``g`` is group-like. This method is
                used to compute the coproduct of any element.

                EXAMPLES::

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: g=CyclicPermutationGroup(6).an_element();g
                    (1,2,3,4,5,6)
                    sage: A.coproduct_on_basis(g)
                    B[(1,2,3,4,5,6)] # B[(1,2,3,4,5,6)]
                    sage: a=A.an_element();a
                    B[()] + 3*B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)]
                    sage: a.coproduct()
                    B[()] # B[()] + 3*B[(1,2,3,4,5,6)] # B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)] # B[(1,3,5)(2,4,6)]
                """
                from sage.categories.tensor import tensor
                g = self.term(g)
                return tensor([g, g])

            def antipode_on_basis(self,g):
                r"""
                Return the antipode of the element ``g`` of the basis.

                Each basis element ``g`` is group-like, and so has
                antipode `g^{-1}`. This method is used to compute the
                antipode of any element.

                EXAMPLES::

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: g=CyclicPermutationGroup(6).an_element();g
                    (1,2,3,4,5,6)
                    sage: A.antipode_on_basis(g)
                    B[(1,6,5,4,3,2)]
                    sage: a=A.an_element();a
                    B[()] + 3*B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)]
                    sage: a.antipode()
                    B[()] + 3*B[(1,5,3)(2,6,4)] + 3*B[(1,6,5,4,3,2)]
                """
                return self.term(~g)

            def counit_on_basis(self,g):
                r"""
                Return the counit of the element ``g`` of the basis.

                Each basis element ``g`` is group-like, and so has
                counit `1`. This method is used to compute the
                counit of any element.

                EXAMPLES::

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: g=CyclicPermutationGroup(6).an_element();g
                    (1,2,3,4,5,6)
                    sage: A.counit_on_basis(g)
                    1
                """
                return self.base_ring().one()

            def counit(self,x):
                r"""
                Return the counit of the element ``x`` of the group
                algebra.

                This is the sum of all coefficients of ``x`` with respect
                to the standard basis of the group algebra.

                EXAMPLES::

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: a=A.an_element();a
                    B[()] + 3*B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)]
                    sage: a.counit()
                    7
                """
                return self.base_ring().sum(x.coefficients())

        class ElementMethods:

            def central_form(self):
                r"""
                Return ``self`` expressed in the canonical basis of the center
                of the group algebra.

                INPUT:

                - ``self`` -- an element of the center of the group algebra

                OUTPUT:

                - A formal linear combination of the conjugacy class
                  representatives representing its coordinates in the
                  canonical basis of the center. See
                  :meth:`Groups.Algebras.ParentMethods.center_basis` for
                  details.

                .. WARNING::

                    - This method requires the underlying group to
                      have a method ``conjugacy_classes_representatives``
                      (every permutation group has one, thanks GAP!).
                    - This method does not check that the element is
                      indeed central. Use the method
                      :meth:`Monoids.Algebras.ElementMethods.is_central`
                      for this purpose.
                    - This function has a complexity linear in the
                      number of conjugacy classes of the group. One
                      could easily implement a function whose
                      complexity is linear in the size of the support
                      of ``self``.

                EXAMPLES::

                    sage: QS3 = SymmetricGroup(3).algebra(QQ)
                    sage: A = QS3([2,3,1]) + QS3([3,1,2])
                    sage: A.central_form()
                    B[(1,2,3)]
                    sage: QS4 = SymmetricGroup(4).algebra(QQ)
                    sage: B = sum(len(s.cycle_type())*QS4(s) for s in Permutations(4))
                    sage: B.central_form()
                    4*B[()] + 3*B[(1,2)] + 2*B[(1,2)(3,4)] + 2*B[(1,2,3)] + B[(1,2,3,4)]

                    sage: QG = GroupAlgebras(QQ).example(PermutationGroup([[(1,2,3),(4,5)],[(3,4)]]))
                    sage: sum(i for i in QG.basis()).central_form()
                    B[()] + B[(4,5)] + B[(3,4,5)] + B[(2,3)(4,5)] + B[(2,3,4,5)] + B[(1,2)(3,4,5)] + B[(1,2,3,4,5)]

                .. SEEALSO::

                    - :meth:`Groups.Algebras.ParentMethods.center_basis`
                    - :meth:`Monoids.Algebras.ElementMethods.is_central`
                """
                from sage.combinat.free_module import CombinatorialFreeModule
                conj_classes_reps = self.parent().basis().keys().conjugacy_classes_representatives()
                Z = CombinatorialFreeModule(self.base_ring(), conj_classes_reps)
                return sum(self[i] * Z.basis()[i] for i in Z.basis().keys())

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

                    sage: C5 = CyclicPermutationGroup(5)
                    sage: C4 = CyclicPermutationGroup(4)
                    sage: S4 = SymmetricGroup(3)
                    sage: C = cartesian_product([C5, C4, S4])
                    sage: C.group_generators()
                    Family (((1,2,3,4,5), (), ()),
                            ((), (1,2,3,4), ()),
                            ((), (), (1,2)),
                            ((), (), (2,3)))

                We check the other portion of :trac:`16718` is fixed::

                    sage: len(C.j_classes())
                    1

                An example with an infinitely generated group (a better output
                is needed)::

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
                    ret = [lift(i, gen) for i,G in enumerate(F) for gen in G.group_generators()]
                    return Family(ret)

                # Infinitely generated
                # This does not return a good output, but it is "correct"
                # TODO: Figure out a better way to do things
                from sage.categories.cartesian_product import cartesian_product
                gens_prod = cartesian_product([Family(G.group_generators(),
                                                      lambda g: (i, g))
                                               for i,G in enumerate(F)])
                return Family(gens_prod, lift, name="gen")

            def order(self):
                r"""
                Return the cardinality of self.

                EXAMPLES::

                    sage: C = cartesian_product([SymmetricGroup(10), SL(2,GF(3))])
                    sage: C.order()
                    87091200

                TESTS::

                    sage: C.order.__module__
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
