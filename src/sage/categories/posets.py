# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.graphs, sage.combinat
r"""
Posets
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category import Category
from sage.categories.sets_cat import Sets


class Posets(Category):
    r"""
    The category of posets i.e. sets with a partial order structure.

    EXAMPLES::

        sage: Posets()
        Category of posets
        sage: Posets().super_categories()
        [Category of sets]
        sage: P = Posets().example(); P
        An example of a poset: sets ordered by inclusion

    The partial order is implemented by the mandatory method
    :meth:`~Posets.ParentMethods.le`::

        sage: x = P(Set([1,3])); y = P(Set([1,2,3]))
        sage: x, y
        ({1, 3}, {1, 2, 3})
        sage: P.le(x, y)
        True
        sage: P.le(x, x)
        True
        sage: P.le(y, x)
        False

    The other comparison methods are called
    :meth:`~Posets.ParentMethods.lt`, :meth:`~Posets.ParentMethods.ge`,
    :meth:`~Posets.ParentMethods.gt`, following Python's naming
    convention in :mod:`operator`. Default implementations are
    provided::

        sage: P.lt(x, x)
        False
        sage: P.ge(y, x)
        True

    Unless the poset is a facade (see :class:`Sets.Facade`), one can
    compare directly its elements using the usual Python operators::

        sage: D = Poset((divisors(30), attrcall("divides")), facade = False)
        sage: D(3) <= D(6)
        True
        sage: D(3) <= D(3)
        True
        sage: D(3) <= D(5)
        False
        sage: D(3) < D(3)
        False
        sage: D(10) >= D(5)
        True

    At this point, this has to be implemented by hand. Once
    :issue:`10130` will be resolved, this will be automatically
    provided by this category::

        sage: # not implemented
        sage: x < y
        True
        sage: x < x
        False
        sage: x <= x
        True
        sage: y >= x
        True

    .. SEEALSO:: :func:`Poset`, :class:`FinitePosets`, :class:`LatticePosets`

    TESTS::

        sage: C = Posets()
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        r"""
        Return a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        EXAMPLES::

            sage: Posets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    def example(self, choice=None):
        r"""
        Return examples of objects of ``Posets()``, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: Posets().example()
            An example of a poset: sets ordered by inclusion

            sage: Posets().example("facade")
            An example of a facade poset:
             the positive integers ordered by divisibility
        """
        from sage.categories.examples.posets import FiniteSetsOrderedByInclusion, PositiveIntegersOrderedByDivisibilityFacade
        if choice == "facade":
            return PositiveIntegersOrderedByDivisibilityFacade()
        else:
            return FiniteSetsOrderedByInclusion()

    def __iter__(self):
        r"""
        Iterator over representatives of the isomorphism classes of
        posets with finitely many vertices.

        .. warning:: this feature may become deprecated, since it does
           of course not iterate through all posets.

        EXAMPLES::

            sage: P = Posets()
            sage: it = iter(P)
            sage: for _ in range(10): print(next(it))
            Finite poset containing 0 elements
            Finite poset containing 1 elements
            Finite poset containing 2 elements
            Finite poset containing 2 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 3 elements
            Finite poset containing 4 elements
        """
        from sage.combinat.posets.posets import FinitePosets_n
        n = 0
        while True:
            yield from FinitePosets_n(n)
            n += 1

    Finite = LazyImport('sage.categories.finite_posets', 'FinitePosets')

    class ParentMethods:

        @abstract_method
        def le(self, x, y):
            r"""
            Return whether `x \le y` in the poset ``self``.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.le( 3, 6 )
                True
                sage: D.le( 3, 3 )
                True
                sage: D.le( 3, 5 )
                False
            """

        def lt(self, x, y):
            r"""
            Return whether `x < y` in the poset ``self``.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            This default implementation delegates the work to :meth:`le`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.lt( 3, 6 )
                True
                sage: D.lt( 3, 3 )
                False
                sage: D.lt( 3, 5 )
                False
            """
            return self.le(x,y) and x != y

        def ge(self, x, y):
            r"""
            Return whether `x \ge y` in the poset ``self``.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            This default implementation delegates the work to :meth:`le`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.ge( 6, 3 )
                True
                sage: D.ge( 3, 3 )
                True
                sage: D.ge( 3, 5 )
                False
            """
            return self.le(y,x)

        def gt(self, x, y):
            r"""
            Return whether `x > y` in the poset ``self``.

            INPUT:

            - ``x``, ``y`` -- elements of ``self``

            This default implementation delegates the work to :meth:`lt`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.gt( 3, 6 )
                False
                sage: D.gt( 3, 3 )
                False
                sage: D.gt( 3, 5 )
                False
            """
            return self.lt(y,x)

        @abstract_method(optional=True)
        def upper_covers(self, x):
            r"""
            Return the upper covers of `x`, that is, the elements `y`
            such that `x<y` and there exists no `z` such that `x<z<y`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.upper_covers(3)
                [6, 15]
            """

        @abstract_method(optional=True)
        def lower_covers(self, x):
            r"""
            Return the lower covers of `x`, that is, the elements `y`
            such that `y<x` and there exists no `z` such that `y<z<x`.

            EXAMPLES::

                sage: D = Poset((divisors(30), attrcall("divides")))
                sage: D.lower_covers(15)
                [3, 5]
            """

        @abstract_method(optional=True)
        def order_ideal(self, elements):
            r"""
            Return the order ideal in ``self`` generated by the elements
            of an iterable ``elements``.

            A subset `I` of a poset is said to be an order ideal if, for
            any `x` in `I` and `y` such that `y \le x`, then `y` is in `I`.

            This is also called the lower set generated by these elements.

            EXAMPLES::

                sage: B = posets.BooleanLattice(4)
                sage: B.order_ideal([7,10])
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
            """

        @abstract_method(optional=True)
        def order_filter(self, elements):
            r"""
            Return the order filter generated by a list of elements.

            A subset `I` of a poset is said to be an order filter if, for
            any `x` in `I` and `y` such that `y \ge x`, then `y` is in `I`.

            This is also called the upper set generated by these elements.

            EXAMPLES::

                sage: B = posets.BooleanLattice(4)
                sage: B.order_filter([3,8])
                [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            """

        def directed_subset(self, elements, direction):
            r"""
            Return the order filter or the order ideal generated by a
            list of elements.

            If ``direction`` is ``'up'``, the order filter (upper set) is
            being returned.

            If ``direction`` is ``'down'``, the order ideal (lower set) is
            being returned.

            INPUT:

            - ``elements`` -- list of elements

            - ``direction`` -- ``'up'`` or ``'down'``

            EXAMPLES::

                sage: B = posets.BooleanLattice(4)
                sage: B.directed_subset([3, 8], 'up')
                [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
                sage: B.directed_subset([7, 10], 'down')
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]

            TESTS::

                sage: B = posets.BooleanLattice(3)
                sage: B.directed_subset([3, 1], 'banana')
                Traceback (most recent call last):
                ...
                ValueError: direction must be either 'up' or 'down'
            """
            if direction == 'up':
                return self.order_filter(elements)
            if direction == 'down':
                return self.order_ideal(elements)
            raise ValueError("direction must be either 'up' or 'down'")

        def principal_order_ideal(self, x):
            r"""
            Return the order ideal generated by an element ``x``.

            This is also called the lower set generated by this element.

            EXAMPLES::

                sage: B = posets.BooleanLattice(4)
                sage: B.principal_order_ideal(6)
                [0, 2, 4, 6]
            """
            return self.order_ideal([x])

        principal_lower_set = principal_order_ideal

        def principal_order_filter(self, x):
            r"""
            Return the order filter generated by an element ``x``.

            This is also called the upper set generated by this element.

            EXAMPLES::

                sage: B = posets.BooleanLattice(4)
                sage: B.principal_order_filter(2)
                [2, 3, 6, 7, 10, 11, 14, 15]
            """
            return self.order_filter([x])

        principal_upper_set = principal_order_filter

        def order_ideal_toggle(self, I, v):
            r"""
            Return the result of toggling the element ``v`` in the
            order ideal ``I``.

            If `v` is an element of a poset `P`, then toggling the
            element `v` is an automorphism of the set `J(P)` of all
            order ideals of `P`. It is defined as follows: If `I`
            is an order ideal of `P`, then the image of `I` under
            toggling the element `v` is

            - the set `I \cup \{ v \}`, if `v \not\in I` but
              every element of `P` smaller than `v` is in `I`;

            - the set `I \setminus \{ v \}`, if `v \in I` but
              no element of `P` greater than `v` is in `I`;

            - `I` otherwise.

            This image always is an order ideal of `P`.

            EXAMPLES::

                sage: P = Poset({1: [2,3], 2: [4], 3: []})
                sage: I = Set({1, 2})
                sage: I in P.order_ideals_lattice()                                     # needs sage.modules
                True
                sage: P.order_ideal_toggle(I, 1)
                {1, 2}
                sage: P.order_ideal_toggle(I, 2)
                {1}
                sage: P.order_ideal_toggle(I, 3)
                {1, 2, 3}
                sage: P.order_ideal_toggle(I, 4)
                {1, 2, 4}
                sage: P4 = Posets(4)
                sage: all(all(all(P.order_ideal_toggle(P.order_ideal_toggle(I, i), i) == I          # needs sage.modules
                ....:               for i in range(4))
                ....:          for I in P.order_ideals_lattice(facade=True))
                ....:     for P in P4)
                True
            """
            if v not in I:
                if all(u in I for u in self.lower_covers(v)):
                    from sage.sets.set import Set
                    return I.union(Set({v}))
            else:
                if all(u not in I for u in self.upper_covers(v)):
                    from sage.sets.set import Set
                    return I.difference(Set({v}))
            return I

        def order_ideal_toggles(self, I, vs):
            r"""
            Return the result of toggling the elements of the list (or
            iterable) ``vs`` (one by one, from left to right) in the order
            ideal ``I``.

            See :meth:`order_ideal_toggle` for a definition of toggling.

            EXAMPLES::

                sage: P = Poset({1: [2,3], 2: [4], 3: []})
                sage: I = Set({1, 2})
                sage: P.order_ideal_toggles(I, [1,2,3,4])
                {1, 3}
                sage: P.order_ideal_toggles(I, (1,2,3,4))
                {1, 3}
            """
            for v in vs:
                I = self.order_ideal_toggle(I, v)
            return I

        def is_order_ideal(self, o):
            """
            Return whether ``o`` is an order ideal of ``self``, assuming ``self``
            has no infinite descending path.

            INPUT:

            - ``o`` -- list (or set, or tuple) containing some elements of ``self``

            EXAMPLES::

                sage: P = Poset((divisors(12), attrcall("divides")),
                ....:           facade=True, linear_extension=True)
                sage: sorted(P.list())
                [1, 2, 3, 4, 6, 12]
                sage: P.is_order_ideal([1, 3])
                True
                sage: P.is_order_ideal([])
                True
                sage: P.is_order_ideal({1, 3})
                True
                sage: P.is_order_ideal([1, 3, 4])
                False
            """
            return all((u in self and all(x in o for x in self.lower_covers(u))) for u in o)

        def is_order_filter(self, o):
            """
            Return whether ``o`` is an order filter of ``self``, assuming ``self``
            has no infinite ascending path.

            INPUT:

            - ``o`` -- list (or set, or tuple) containing some elements of ``self``

            EXAMPLES::

                sage: P = Poset((divisors(12), attrcall("divides")),
                ....:           facade=True, linear_extension=True)
                sage: sorted(P.list())
                [1, 2, 3, 4, 6, 12]
                sage: P.is_order_filter([4, 12])
                True
                sage: P.is_order_filter([])
                True
                sage: P.is_order_filter({3, 4, 12})
                False
                sage: P.is_order_filter({3, 6, 12})
                True
            """
            return all((u in self and all(x in o for x in self.upper_covers(u))) for u in o)

        def is_chain_of_poset(self, o, ordered=False):
            """
            Return whether an iterable ``o`` is a chain of ``self``,
            including a check for ``o`` being ordered from smallest
            to largest element if the keyword ``ordered`` is set to
            ``True``.

            INPUT:

            - ``o`` -- an iterable (e. g., list, set, or tuple)
              containing some elements of ``self``

            - ``ordered`` -- a Boolean (default: ``False``); decides
              whether the notion of a chain includes being ordered

            OUTPUT:

            If ``ordered`` is set to ``False``, the truth value of
            the following assertion is returned: The subset of ``self``
            formed by the elements of ``o`` is a chain in ``self``.

            If ``ordered`` is set to ``True``, the truth value of
            the following assertion is returned: Every element of the
            list ``o`` is (strictly!) smaller than its successor in
            ``self``. (This makes no sense if ``ordered`` is a set.)

            EXAMPLES::

                sage: P = Poset((divisors(12), attrcall("divides")),
                ....:           facade=True, linear_extension=True)
                sage: sorted(P.list())
                [1, 2, 3, 4, 6, 12]
                sage: P.is_chain_of_poset([1, 3])
                True
                sage: P.is_chain_of_poset([3, 1])
                True
                sage: P.is_chain_of_poset([1, 3], ordered=True)
                True
                sage: P.is_chain_of_poset([3, 1], ordered=True)
                False
                sage: P.is_chain_of_poset([])
                True
                sage: P.is_chain_of_poset([], ordered=True)
                True
                sage: P.is_chain_of_poset((2, 12, 6))
                True
                sage: P.is_chain_of_poset((2, 6, 12), ordered=True)
                True
                sage: P.is_chain_of_poset((2, 12, 6), ordered=True)
                False
                sage: P.is_chain_of_poset((2, 12, 6, 3))
                False
                sage: P.is_chain_of_poset((2, 3))
                False

                sage: Q = Poset({2: [3, 1], 3: [4], 1: [4]})
                sage: Q.is_chain_of_poset([1, 2], ordered=True)
                False
                sage: Q.is_chain_of_poset([1, 2])
                True
                sage: Q.is_chain_of_poset([2, 1], ordered=True)
                True
                sage: Q.is_chain_of_poset([2, 1, 1], ordered=True)
                False
                sage: Q.is_chain_of_poset([3])
                True
                sage: Q.is_chain_of_poset([4, 2, 3])
                True
                sage: Q.is_chain_of_poset([4, 2, 3], ordered=True)
                False
                sage: Q.is_chain_of_poset([2, 3, 4], ordered=True)
                True

            Examples with infinite posets::

                sage: from sage.categories.examples.posets import FiniteSetsOrderedByInclusion
                sage: R = FiniteSetsOrderedByInclusion()
                sage: R.is_chain_of_poset([R(set([3, 1, 2])),
                ....:                      R(set([1, 4])),
                ....:                      R(set([4, 5]))])
                False
                sage: R.is_chain_of_poset([R(set([3, 1, 2])),
                ....:                      R(set([1, 2])),
                ....:                      R(set([1]))], ordered=True)
                False
                sage: R.is_chain_of_poset([R(set([3, 1, 2])),
                ....:                      R(set([1, 2])), R(set([1]))])
                True

                sage: from sage.categories.examples.posets import PositiveIntegersOrderedByDivisibilityFacade
                sage: T = PositiveIntegersOrderedByDivisibilityFacade()
                sage: T.is_chain_of_poset((T(3), T(4), T(7)))
                False
                sage: T.is_chain_of_poset((T(3), T(6), T(3)))
                True
                sage: T.is_chain_of_poset((T(3), T(6), T(3)), ordered=True)
                False
                sage: T.is_chain_of_poset((T(3), T(3), T(6)))
                True
                sage: T.is_chain_of_poset((T(3), T(3), T(6)), ordered=True)
                False
                sage: T.is_chain_of_poset((T(3), T(6)), ordered=True)
                True
                sage: T.is_chain_of_poset((), ordered=True)
                True
                sage: T.is_chain_of_poset((T(3),), ordered=True)
                True
                sage: T.is_chain_of_poset((T(q) for q in divisors(27)))
                True
                sage: T.is_chain_of_poset((T(q) for q in divisors(18)))
                False
            """
            list_o = list(o)
            if ordered:
                return all(self.lt(a, b) for a, b in zip(list_o, list_o[1:]))
            else:
                for (i, x) in enumerate(list_o):
                    for y in list_o[:i]:
                        if (not self.le(x, y)) and (not self.gt(x, y)):
                            return False
                return True

        def is_antichain_of_poset(self, o):
            """
            Return whether an iterable ``o`` is an antichain of
            ``self``.

            INPUT:

            - ``o`` -- an iterable (e. g., list, set, or tuple)
              containing some elements of ``self``

            OUTPUT:

            ``True`` if the subset of ``self`` consisting of the entries
            of ``o`` is an antichain of ``self``, and ``False`` otherwise.

            EXAMPLES::

                sage: P = Poset((divisors(12), attrcall("divides")),
                ....:           facade=True, linear_extension=True)
                sage: sorted(P.list())
                [1, 2, 3, 4, 6, 12]
                sage: P.is_antichain_of_poset([1, 3])
                False
                sage: P.is_antichain_of_poset([3, 1])
                False
                sage: P.is_antichain_of_poset([1, 1, 3])
                False
                sage: P.is_antichain_of_poset([])
                True
                sage: P.is_antichain_of_poset([1])
                True
                sage: P.is_antichain_of_poset([1, 1])
                True
                sage: P.is_antichain_of_poset([3, 4])
                True
                sage: P.is_antichain_of_poset([3, 4, 12])
                False
                sage: P.is_antichain_of_poset([6, 4])
                True
                sage: P.is_antichain_of_poset(i for i in divisors(12)
                ....:                         if (2 < i and i < 6))
                True
                sage: P.is_antichain_of_poset(i for i in divisors(12)
                ....:                         if (2 <= i and i < 6))
                False

                sage: Q = Poset({2: [3, 1], 3: [4], 1: [4]})
                sage: Q.is_antichain_of_poset((1, 2))
                False
                sage: Q.is_antichain_of_poset((2, 4))
                False
                sage: Q.is_antichain_of_poset((4, 2))
                False
                sage: Q.is_antichain_of_poset((2, 2))
                True
                sage: Q.is_antichain_of_poset((3, 4))
                False
                sage: Q.is_antichain_of_poset((3, 1))
                True
                sage: Q.is_antichain_of_poset((1, ))
                True
                sage: Q.is_antichain_of_poset(())
                True

            An infinite poset::

                sage: from sage.categories.examples.posets import FiniteSetsOrderedByInclusion
                sage: R = FiniteSetsOrderedByInclusion()
                sage: R.is_antichain_of_poset([R(set([3, 1, 2])),
                ....:                          R(set([1, 4])), R(set([4, 5]))])
                True
                sage: R.is_antichain_of_poset([R(set([3, 1, 2, 4])),
                ....:                          R(set([1, 4])), R(set([4, 5]))])
                False
            """
            return all(not self.lt(x,y) for x in o for y in o)

        CartesianProduct = LazyImport(
            'sage.combinat.posets.cartesian_product', 'CartesianProductPoset')

    class ElementMethods:
        pass
        # TODO: implement x<y, x<=y, x>y, x>=y appropriately once #10130 is resolved
        #
        # def __le__(self, other):
        #     r"""
        #     Return whether ``self`` is smaller or equal to ``other``
        #     in the poset.
        #
        #     EXAMPLES::
        #
        #         sage: P = Posets().example(); P
        #         An example of poset: sets ordered by inclusion
        #         sage: x = P(Set([1,3])); y = P(Set([1,2,3]))
        #         sage: x.__le__(y)
        #         sage: x <= y
        #     """
        #     return self.parent().le(self, other)
