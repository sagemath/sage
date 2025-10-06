"""
Cartesian products of posets

AUTHORS:

- Daniel Krenn (2015)
"""
# ****************************************************************************
#  Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                https://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.cartesian_product import CartesianProduct


class CartesianProductPoset(CartesianProduct):
    r"""
    A class implementing Cartesian products of posets (and elements
    thereof). Compared to :class:`CartesianProduct` you are able to
    specify an order for comparison of the elements.

    INPUT:

    - ``sets`` -- tuple of parents

    - ``category`` -- a subcategory of ``Sets().CartesianProducts() & Posets()``

    - ``order`` -- string or function specifying an order less or equal;
      it can be one of the following:

      - ``'native'`` -- elements are ordered by their native ordering,
        i.e., the order the wrapped elements (tuples) provide

      - ``'lex'`` -- elements are ordered lexicographically

      - ``'product'`` -- an element is less or equal to another
        element, if less or equal is true for all its components
        (Cartesian projections)

      - a function which performs the comparison `\leq`; it takes two
        input arguments and outputs a boolean

    Other keyword arguments (``kwargs``) are passed to the constructor
    of :class:`CartesianProduct`.

    EXAMPLES::

        sage: P = Poset((srange(3), lambda left, right: left <= right))
        sage: Cl = cartesian_product((P, P), order='lex')
        sage: Cl((1, 1)) <= Cl((2, 0))
        True
        sage: Cp = cartesian_product((P, P), order='product')
        sage: Cp((1, 1)) <= Cp((2, 0))
        False
        sage: def le_sum(left, right):
        ....:     return (sum(left) < sum(right) or
        ....:             sum(left) == sum(right) and left[0] <= right[0])
        sage: Cs = cartesian_product((P, P), order=le_sum)
        sage: Cs((1, 1)) <= Cs((2, 0))
        True

    TESTS::

        sage: Cl.category()
        Join of Category of finite posets and
        Category of Cartesian products of finite enumerated sets
        sage: TestSuite(Cl).run(skip=['_test_construction'])
        sage: Cp.category()
        Join of Category of finite posets and
        Category of Cartesian products of finite enumerated sets
        sage: TestSuite(Cp).run(skip=['_test_construction'])

    .. SEEALSO::

        :class:`CartesianProduct`
    """

    def __init__(self, sets, category, order=None, **kwargs) -> None:
        r"""
        See :class:`CartesianProductPoset` for details.

        TESTS::

            sage: P = Poset((srange(3), lambda left, right: left <= right))
            sage: C = cartesian_product((P, P), order='notexisting')
            Traceback (most recent call last):
            ...
            ValueError: no order 'notexisting' known
            sage: C = cartesian_product((P, P), category=(Groups(),))
            sage: C.category()
            Join of Category of groups and Category of posets
        """
        if order is None:
            self._le_ = self.le_product
        elif isinstance(order, str):
            try:
                self._le_ = getattr(self, 'le_' + order)
            except AttributeError:
                raise ValueError(f"no order '{order}' known")
        else:
            self._le_ = order

        from sage.categories.category import Category
        from sage.categories.posets import Posets
        if not isinstance(category, tuple):
            category = (category,)
        category = Category.join(category + (Posets(),))
        super().__init__(sets, category, **kwargs)

    def le(self, left, right):
        r"""
        Test whether ``left`` is less than or equal to ``right``.

        INPUT:

        - ``left`` -- an element

        - ``right`` -- an element

        OUTPUT: boolean

        .. NOTE::

            This method uses the order defined on creation of this
            Cartesian product. See :class:`CartesianProductPoset`.

        EXAMPLES::

            sage: P = posets.ChainPoset(10)
            sage: def le_sum(left, right):
            ....:     return (sum(left) < sum(right) or
            ....:             sum(left) == sum(right) and left[0] <= right[0])
            sage: C = cartesian_product((P, P), order=le_sum)
            sage: C.le(C((1, 6)), C((6, 1)))
            True
            sage: C.le(C((6, 1)), C((1, 6)))
            False
            sage: C.le(C((1, 6)), C((6, 6)))
            True
            sage: C.le(C((6, 6)), C((1, 6)))
            False
        """
        return self._le_(left, right)

    def le_lex(self, left, right):
        r"""
        Test whether ``left`` is lexicographically smaller or equal
        to ``right``.

        INPUT:

        - ``left`` -- an element

        - ``right`` -- an element

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Poset((srange(2), lambda left, right: left <= right))
            sage: Q = cartesian_product((P, P), order='lex')
            sage: T = [Q((0, 0)), Q((1, 1)), Q((0, 1)), Q((1, 0))]
            sage: for a in T:
            ....:     for b in T:
            ....:         assert Q.le(a, b) == (a <= b)
            ....:         print('%s <= %s = %s' % (a, b, a <= b))
            (0, 0) <= (0, 0) = True
            (0, 0) <= (1, 1) = True
            (0, 0) <= (0, 1) = True
            (0, 0) <= (1, 0) = True
            (1, 1) <= (0, 0) = False
            (1, 1) <= (1, 1) = True
            (1, 1) <= (0, 1) = False
            (1, 1) <= (1, 0) = False
            (0, 1) <= (0, 0) = False
            (0, 1) <= (1, 1) = True
            (0, 1) <= (0, 1) = True
            (0, 1) <= (1, 0) = True
            (1, 0) <= (0, 0) = False
            (1, 0) <= (1, 1) = True
            (1, 0) <= (0, 1) = False
            (1, 0) <= (1, 0) = True

        TESTS:

        Check that :issue:`19999` is resolved::

            sage: P = Poset((srange(2), lambda left, right: left <= right))
            sage: Q = cartesian_product((P, P), order='product')
            sage: R = cartesian_product((Q, P), order='lex')
            sage: R(((1, 0), 0)) <= R(((0, 1), 0))
            False
            sage: R(((0, 1), 0)) <= R(((1, 0), 0))
            False
        """
        for l, r, S in \
                zip(left.value, right.value, self.cartesian_factors()):
            if l == r:
                continue
            if S.le(l, r):
                return True
            if S.le(r, l):
                return False
            return False  # incomparable components
        return True  # equal

    def le_product(self, left, right):
        r"""
        Test whether ``left`` is component-wise smaller or equal
        to ``right``.

        INPUT:

        - ``left`` -- an element

        - ``right`` -- an element

        OUTPUT: boolean

        The comparison is ``True`` if the result of the
        comparison in each component is ``True``.

        EXAMPLES::

            sage: P = Poset((srange(2), lambda left, right: left <= right))
            sage: Q = cartesian_product((P, P), order='product')
            sage: T = [Q((0, 0)), Q((1, 1)), Q((0, 1)), Q((1, 0))]
            sage: for a in T:
            ....:     for b in T:
            ....:         assert Q.le(a, b) == (a <= b)
            ....:         print('%s <= %s = %s' % (a, b, a <= b))
            (0, 0) <= (0, 0) = True
            (0, 0) <= (1, 1) = True
            (0, 0) <= (0, 1) = True
            (0, 0) <= (1, 0) = True
            (1, 1) <= (0, 0) = False
            (1, 1) <= (1, 1) = True
            (1, 1) <= (0, 1) = False
            (1, 1) <= (1, 0) = False
            (0, 1) <= (0, 0) = False
            (0, 1) <= (1, 1) = True
            (0, 1) <= (0, 1) = True
            (0, 1) <= (1, 0) = False
            (1, 0) <= (0, 0) = False
            (1, 0) <= (1, 1) = True
            (1, 0) <= (0, 1) = False
            (1, 0) <= (1, 0) = True
        """
        return all(
            S.le(l, r)
            for l, r, S in
            zip(left.value, right.value, self.cartesian_factors()))

    def le_native(self, left, right):
        r"""
        Test whether ``left`` is smaller or equal to ``right`` in the order
        provided by the elements themselves.

        INPUT:

        - ``left`` -- an element

        - ``right`` -- an element

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Poset((srange(2), lambda left, right: left <= right))
            sage: Q = cartesian_product((P, P), order='native')
            sage: T = [Q((0, 0)), Q((1, 1)), Q((0, 1)), Q((1, 0))]
            sage: for a in T:
            ....:     for b in T:
            ....:         assert Q.le(a, b) == (a <= b)
            ....:         print('%s <= %s = %s' % (a, b, a <= b))
            (0, 0) <= (0, 0) = True
            (0, 0) <= (1, 1) = True
            (0, 0) <= (0, 1) = True
            (0, 0) <= (1, 0) = True
            (1, 1) <= (0, 0) = False
            (1, 1) <= (1, 1) = True
            (1, 1) <= (0, 1) = False
            (1, 1) <= (1, 0) = False
            (0, 1) <= (0, 0) = False
            (0, 1) <= (1, 1) = True
            (0, 1) <= (0, 1) = True
            (0, 1) <= (1, 0) = True
            (1, 0) <= (0, 0) = False
            (1, 0) <= (1, 1) = True
            (1, 0) <= (0, 1) = False
            (1, 0) <= (1, 0) = True
        """
        return left.value <= right.value

    class Element(CartesianProduct.Element):

        def _le_(self, other):
            r"""
            Return if this element is less or equal to ``other``.

            INPUT:

            - ``other`` -- an element

            OUTPUT: boolean

            .. NOTE::

                This method calls :meth:`CartesianProductPoset.le`. Override
                it in inherited class to change this.

                It can be assumed that this element and ``other`` have
                the same parent.

            TESTS::

                sage: from sage.combinat.posets.cartesian_product import CartesianProductPoset
                sage: QQ.CartesianProduct = CartesianProductPoset  # needed until #19269 is fixed
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) <= C((2, 1/3))  # indirect doctest
                True
                sage: C((1/3, 2)) <= C((2, 2))  # indirect doctest
                True
            """
            return self.parent().le(self, other)

        def __le__(self, other):
            r"""
            Return if this element is less than or equal to ``other``.

            INPUT:

            - ``other`` -- an element

            OUTPUT: boolean

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :issue:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: from sage.combinat.posets.cartesian_product import CartesianProductPoset
                sage: QQ.CartesianProduct = CartesianProductPoset  # needed until #19269 is fixed
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) <= C((2, 1/3))
                True
                sage: C((1/3, 2)) <= C((2, 2))
                True

            The following example tests that the coercion gets involved in
            comparisons; it can be simplified once :issue:`18182` is merged.
            ::

                sage: class MyCP(CartesianProductPoset):
                ....:     def _coerce_map_from_(self, S):
                ....:         if isinstance(S, self.__class__):
                ....:             S_factors = S.cartesian_factors()
                ....:             R_factors = self.cartesian_factors()
                ....:             if len(S_factors) == len(R_factors):
                ....:                 if all(r.has_coerce_map_from(s)
                ....:                        for r,s in zip(R_factors, S_factors)):
                ....:                     return True
                sage: QQ.CartesianProduct = MyCP
                sage: A = cartesian_product((QQ, ZZ), order=le_sum)
                sage: B = cartesian_product((QQ, QQ), order=le_sum)
                sage: A((1/2, 4)) <= B((1/2, 5))
                True
            """
            from sage.structure.element import have_same_parent
            if have_same_parent(self, other):
                return self._le_(other)

            from sage.structure.element import get_coercion_model
            import operator
            try:
                return get_coercion_model().bin_op(self, other, operator.le)
            except TypeError:
                return False

        def __ge__(self, other):
            r"""
            Return if this element is greater than or equal to ``other``.

            INPUT:

            - ``other`` -- an element

            OUTPUT: boolean

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :issue:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: from sage.combinat.posets.cartesian_product import CartesianProductPoset
                sage: QQ.CartesianProduct = CartesianProductPoset  # needed until #19269 is fixed
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) >= C((2, 1/3))
                False
                sage: C((1/3, 2)) >= C((2, 2))
                False
            """
            return other <= self

        def __lt__(self, other):
            r"""
            Return if this element is less than ``other``.

            INPUT:

            - ``other`` -- an element

            OUTPUT: boolean

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :issue:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: from sage.combinat.posets.cartesian_product import CartesianProductPoset
                sage: QQ.CartesianProduct = CartesianProductPoset  # needed until #19269 is fixed
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) < C((2, 1/3))
                True
                sage: C((1/3, 2)) < C((2, 2))
                True
            """
            return not self == other and self <= other

        def __gt__(self, other):
            r"""
            Return if this element is greater than ``other``.

            INPUT:

            - ``other`` -- an element

            OUTPUT: boolean

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :issue:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: from sage.combinat.posets.cartesian_product import CartesianProductPoset
                sage: QQ.CartesianProduct = CartesianProductPoset  # needed until #19269 is fixed
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) > C((2, 1/3))
                False
                sage: C((1/3, 2)) > C((2, 2))
                False
            """
            return not self == other and other <= self
