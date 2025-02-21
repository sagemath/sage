# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.groups
r"""
Finite groups
"""
# ****************************************************************************
#  Copyright (C) 2010-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.algebra_functor import AlgebrasCategory


class FiniteGroups(CategoryWithAxiom):
    r"""
    The category of finite (multiplicative) groups.

    EXAMPLES::

        sage: C = FiniteGroups(); C
        Category of finite groups
        sage: C.super_categories()
        [Category of finite monoids, Category of groups]
        sage: C.example()
        General Linear Group of degree 2 over Finite Field of size 3

    TESTS::

        sage: TestSuite(C).run()
    """

    def example(self):
        """
        Return an example of finite group, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: G = FiniteGroups().example(); G
            General Linear Group of degree 2 over Finite Field of size 3
        """
        from sage.groups.matrix_gps.linear import GL
        return GL(2,3)

    class ParentMethods:

        def semigroup_generators(self):
            """
            Return semigroup generators for ``self``.

            For finite groups, the group generators are also semigroup
            generators. Hence, this default implementation calls
            :meth:`~sage.categories.groups.Groups.ParentMethods.group_generators`.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.semigroup_generators()
                Family ((1,2,3), (2,3,4))
            """
            return self.group_generators()

        def monoid_generators(self):
            """
            Return monoid generators for ``self``.

            For finite groups, the group generators are also monoid
            generators. Hence, this default implementation calls
            :meth:`~sage.categories.groups.Groups.ParentMethods.group_generators`.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.monoid_generators()
                Family ((1,2,3), (2,3,4))
            """
            return self.group_generators()

        def cardinality(self):
            """
            Return the cardinality of ``self``, as per
            :meth:`EnumeratedSets.ParentMethods.cardinality`.

            This default implementation calls :meth:`.order` if
            available, and otherwise resorts to
            :meth:`._cardinality_from_iterator`. This is for backward
            compatibility only. Finite groups should override this
            method instead of :meth:`.order`.

            EXAMPLES:

            We need to use a finite group which uses this default
            implementation of cardinality::

                sage: G = groups.misc.SemimonomialTransformation(GF(5), 3); G           # needs sage.rings.number_field
                Semimonomial transformation group over Finite Field of size 5 of degree 3
                sage: G.cardinality.__module__                                          # needs sage.rings.number_field
                'sage.categories.finite_groups'
                sage: G.cardinality()                                                   # needs sage.rings.number_field
                384
            """
            if hasattr(self, 'order'):
                return self.order()
            return self._cardinality_from_iterator()

        def some_elements(self):
            """
            Return some elements of ``self``.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.some_elements()
                Family ((1,2,3), (2,3,4))
            """
            return self.group_generators()

        # TODO: merge with that of finite semigroups
        def cayley_graph_disabled(self, connecting_set=None):
            """

            AUTHORS:

            - Bobby Moretti (2007-08-10)

            - Robert Miller (2008-05-01): editing
            """
            if connecting_set is None:
                connecting_set = self.group_generators()
            else:
                for g in connecting_set:
                    if g not in self:
                        raise RuntimeError("each element of the connecting set must be in the group")
                connecting_set = [self(g) for g in connecting_set]
            from sage.graphs.digraph import DiGraph
            arrows = {}
            for x in self:
                arrows[x] = {}
                for g in connecting_set:
                    xg = x*g # cache the multiplication
                    if not xg == x:
                        arrows[x][xg] = g

            return DiGraph(arrows, implementation='networkx')

        def conjugacy_classes(self):
            r"""
            Return a list with all the conjugacy classes of the group.

            This will eventually be a fall-back method for groups not defined
            over GAP. Right now, it just raises a
            :exc:`NotImplementedError`, until we include a non-GAP
            way of listing the conjugacy classes representatives.

            EXAMPLES::

                sage: from sage.groups.group import FiniteGroup
                sage: G = FiniteGroup()
                sage: G.conjugacy_classes()
                Traceback (most recent call last):
                ...
                NotImplementedError: Listing the conjugacy classes for group <sage.groups.group.FiniteGroup object at ...> is not implemented
            """
            raise NotImplementedError("Listing the conjugacy classes for group %s is not implemented" % self)

        def conjugacy_classes_representatives(self):
            r"""
            Return a list of the conjugacy classes representatives of the group.

            EXAMPLES::

                sage: G = SymmetricGroup(3)
                sage: G.conjugacy_classes_representatives()                             # needs sage.combinat
                [(), (1,2), (1,2,3)]
            """
            return [C.representative() for C in self.conjugacy_classes()]

    class ElementMethods:
        pass

    class Algebras(AlgebrasCategory):
        def extra_super_categories(self):
            r"""
            Implement Maschke's theorem.

            In characteristic 0 all finite group algebras are semisimple.

            EXAMPLES::

                sage: FiniteGroups().Algebras(QQ).is_subcategory(Algebras(QQ).Semisimple())
                True
                sage: FiniteGroups().Algebras(FiniteField(7)).is_subcategory(Algebras(FiniteField(7)).Semisimple())
                False
                sage: FiniteGroups().Algebras(ZZ).is_subcategory(Algebras(ZZ).Semisimple())
                False
                sage: FiniteGroups().Algebras(Fields()).is_subcategory(Algebras(Fields()).Semisimple())
                False

                sage: Cat = CommutativeAdditiveGroups().Finite()
                sage: Cat.Algebras(QQ).is_subcategory(Algebras(QQ).Semisimple())
                True
                sage: Cat.Algebras(GF(7)).is_subcategory(Algebras(GF(7)).Semisimple())
                False
                sage: Cat.Algebras(ZZ).is_subcategory(Algebras(ZZ).Semisimple())
                False
                sage: Cat.Algebras(Fields()).is_subcategory(Algebras(Fields()).Semisimple())
                False
            """
            from sage.categories.fields import Fields
            K = self.base_ring()
            if K in Fields() and K.characteristic() == 0:
                from sage.categories.algebras import Algebras
                return [Algebras(self.base_ring()).Semisimple()]
            else:
                return []

        class ParentMethods:
            def __init_extra__(self):
                r"""
                Implement Maschke's theorem.

                EXAMPLES::

                    sage: G = groups.permutation.Dihedral(8)
                    sage: A = G.algebra(GF(5))
                    sage: A in Algebras.Semisimple
                    True
                    sage: A = G.algebra(Zmod(4))
                    sage: A in Algebras.Semisimple
                    False

                    sage: G = groups.misc.AdditiveCyclic(4)                             # needs sage.rings.number_field
                    sage: Cat = CommutativeAdditiveGroups().Finite()
                    sage: A = G.algebra(GF(5), category=Cat)
                    sage: A in Algebras.Semisimple
                    True
                    sage: A = G.algebra(GF(2), category=Cat)
                    sage: A in Algebras.Semisimple
                    False
                """
                base_ring = self.base_ring()
                group = self.group()
                from sage.categories.fields import Fields
                # If base_ring is of characteristic 0, this is handled
                #    in the FiniteGroups.Algebras category
                # Maschke's theorem: under some conditions, the algebra is semisimple.
                if (base_ring in Fields
                    and base_ring.characteristic() > 0
                    and hasattr(group, "cardinality")
                    and group.cardinality() % base_ring.characteristic() != 0):
                    self._refine_category_(self.category().Semisimple())
