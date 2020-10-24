r"""
Affine Weyl groups
"""
# *****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.weyl_groups import WeylGroups
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets



class AffineWeylGroups(Category_singleton):
    """
    The category of affine Weyl groups

    .. todo:: add a description of this category

    .. SEEALSO::

        - :wikipedia:`Affine_weyl_group`
        - :class:`WeylGroups`, :class:`WeylGroup`

    EXAMPLES::

        sage: C = AffineWeylGroups(); C
        Category of affine weyl groups
        sage: C.super_categories()
        [Category of infinite weyl groups]

        sage: C.example()
        NotImplemented
        sage: W = WeylGroup(["A",4,1]); W
        Weyl Group of type ['A', 4, 1] (as a matrix group acting on the root space)
        sage: W.category()
        Category of irreducible affine weyl groups

    TESTS::

        sage: TestSuite(C).run()
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: AffineWeylGroups().super_categories()
            [Category of infinite weyl groups]
        """
        return [WeylGroups().Infinite()]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of affine Weyl groups defines no
        additional structure: affine Weyl groups are a special class
        of Weyl groups.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`CategoryWithAxiom`?

        EXAMPLES::

            sage: AffineWeylGroups().additional_structure()
        """
        return None

    class ParentMethods:

        @cached_method
        def special_node(self):
            """
            Return the distinguished special node of the underlying
            Dynkin diagram.

            EXAMPLES::

                sage: W = WeylGroup(['A',3,1])
                sage: W.special_node()
                0
            """
            return self.cartan_type().special_node()

        def affine_grassmannian_elements_of_given_length(self, k):
            """
            Return the affine Grassmannian elements of length `k`.

            This is returned as a finite enumerated set.

            EXAMPLES::

                sage: W = WeylGroup(['A',3,1])
                sage: [x.reduced_word() for x in W.affine_grassmannian_elements_of_given_length(3)]
                [[2, 1, 0], [3, 1, 0], [2, 3, 0]]

            .. SEEALSO::

                :meth:`AffineWeylGroups.ElementMethods.is_affine_grassmannian`
            """
            from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest

            def select_length(pair):
                u, length = pair
                if length == k:
                    return u

            def succ(pair):
                u, length = pair
                for i in u.descents(positive=True, side="left"):
                    u1 = u.apply_simple_reflection(i, "left")
                    if (length < k and i == u1.first_descent(side="left") and
                            u1.is_affine_grassmannian()):
                        yield (u1, length + 1)
                return
            return RecursivelyEnumeratedSet_forest(((self.one(), 0),), succ, algorithm='breadth',
                                                   category=FiniteEnumeratedSets(),
                                                   post_process=select_length)

    class ElementMethods:
        def is_affine_grassmannian(self):
            """
            Test whether ``self`` is affine Grassmannian.

            An element of an affine Weyl group is *affine Grassmannian*
            if any of the following equivalent properties holds:

            - all reduced words for ``self`` end with 0.
            - ``self`` is the identity, or 0 is its single right descent.
            - ``self`` is a minimal coset representative for W / cl W.

            EXAMPLES::

                sage: W = WeylGroup(['A',3,1])
                sage: w = W.from_reduced_word([2,1,0])
                sage: w.is_affine_grassmannian()
                True
                sage: w = W.from_reduced_word([2,0])
                sage: w.is_affine_grassmannian()
                False
                sage: W.one().is_affine_grassmannian()
                True
            """
            D = self.descents()
            return (not D) or D == [self.parent().special_node()]

        def affine_grassmannian_to_core(self):
            r"""
            Bijection between affine Grassmannian elements of type `A_k^{(1)}` and `(k+1)`-cores.

            INPUT:

            - ``self`` -- an affine Grassmannian element of some affine Weyl group of type `A_k^{(1)}`

            Recall that an element `w` of an affine Weyl group is
            affine Grassmannian if all its all reduced words end in 0, see :meth:`is_affine_grassmannian`.

            OUTPUT:

            - a `(k+1)`-core

            .. SEEALSO:: :meth:`affine_grassmannian_to_partition`

            EXAMPLES::

                sage: W = WeylGroup(['A',2,1])
                sage: w = W.from_reduced_word([0,2,1,0])
                sage: la = w.affine_grassmannian_to_core(); la
                [4, 2]
                sage: type(la)
                <class 'sage.combinat.core.Cores_length_with_category.element_class'>
                sage: la.to_grassmannian() == w
                True

                sage: w = W.from_reduced_word([0,2,1])
                sage: w.affine_grassmannian_to_core()
                Traceback (most recent call last):
                ...
                ValueError: this only works on type 'A' affine Grassmannian elements
            """
            from sage.combinat.partition import Partition
            from sage.combinat.core import Core
            if not self.is_affine_grassmannian() or not self.parent().cartan_type().letter == 'A':
                raise ValueError("this only works on type 'A' affine Grassmannian elements")
            out = Partition([])
            rword = self.reduced_word()
            kp1 = self.parent().n
            for i in range(len(rword)):
                for c in out.outside_corners():
                    if (c[1] - c[0]) % kp1 == rword[-i - 1]:
                        out = out.add_cell(c[0], c[1])
            return Core(out._list, kp1)

        def affine_grassmannian_to_partition(self):
            r"""
            Bijection between affine Grassmannian elements of type `A_k^{(1)}`
            and `k`-bounded partitions.

            INPUT:

            - ``self`` is affine Grassmannian element of the affine Weyl group of type `A_k^{(1)}` (i.e. all reduced words end in 0)

            OUTPUT:

            - `k`-bounded partition

            .. SEEALSO:: :meth:`affine_grassmannian_to_core`

            EXAMPLES::

                sage: k = 2
                sage: W = WeylGroup(['A',k,1])
                sage: w = W.from_reduced_word([0,2,1,0])
                sage: la = w.affine_grassmannian_to_partition(); la
                [2, 2]
                sage: la.from_kbounded_to_grassmannian(k) == w
                True
            """
            return self.affine_grassmannian_to_core().to_bounded_partition()
