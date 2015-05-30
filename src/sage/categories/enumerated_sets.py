r"""
Enumerated Sets
"""
#*****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_singleton import Category_singleton
from sage.categories.sets_cat import Sets
from sage.categories.sets_cat import EmptySetError
from sage.categories.cartesian_product import CartesianProductsCategory

class EnumeratedSets(Category_singleton):
    """
    The category of enumerated sets

    An *enumerated set* is a *finite* or *countable* set or multiset `S`
    together with a canonical enumeration of its elements;
    conceptually, this is very similar to an immutable list. The main
    difference lies in the names and the return type of the methods,
    and of course the fact that the list of element is not supposed to
    be expanded in memory. Whenever possible one should use one of the
    two sub-categories :class:`FiniteEnumeratedSets` or
    :class:`InfiniteEnumeratedSets`.

    The purpose of this category is threefold:

     - to fix a common interface for all these sets;
     - to provide a bunch of default implementations;
     - to provide consistency tests.

    The standard methods for an enumerated set ``S`` are:

       - ``S.cardinality()``: the number of elements of the set. This
         is the equivalent for ``len`` on a list except that the
         return value is specified to be a Sage :class:`Integer` or
         ``infinity``, instead of a Python ``int``;

       - ``iter(S)``: an iterator for the elements of the set;

       - ``S.list()``: the list of the elements of the set, when
         possible; raises a NotImplementedError if the list is
         predictably too large to be expanded in memory.

       - ``S.unrank(n)``: the  ``n-th`` element of the set when ``n`` is a sage
         ``Integer``. This is the equivanlent for ``l[n]`` on a list.

       - ``S.rank(e)``: the position of the element ``e`` in the set;
         This is equivalent to ``l.index(e)`` for a list except that
         the return value is specified to be a Sage :class:`Integer`,
         instead of a Python ``int``;

       - ``S.first()``: the first object of the set; it is equivalent to
         ``S.unrank(0)``;

       - ``S.next(e)``: the object of the set which follows ``e``; It is
         equivalent to ``S.unrank(S.rank(e)+1)``.

       - ``S.random_element()``: a random generator for an element of
         the set. Unless otherwise stated, and for finite enumerated
         sets, the probability is uniform.

    For examples and tests see:

       - ``FiniteEnumeratedSets().example()``
       - ``InfiniteEnumeratedSets().example()``


    EXAMPLES::

        sage: EnumeratedSets()
        Category of enumerated sets
        sage: EnumeratedSets().super_categories()
        [Category of sets]
        sage: EnumeratedSets().all_super_categories()
        [Category of enumerated sets, Category of sets, Category of sets with partial maps, Category of objects]

    TESTS::

        sage: C = EnumeratedSets()
        sage: TestSuite(C).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: EnumeratedSets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    def additional_structure(self):
        """
        Return ``None``.

        Indeed, morphisms of enumerated sets are not required to
        preserve the enumeration.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: EnumeratedSets().additional_structure()
        """
        return None

    def _call_(self, X):
        """
        Construct an object in this category from the data in ``X``.

        EXAMPLES::

            sage: EnumeratedSets()(Primes())
            Set of all prime numbers: 2, 3, 5, 7, ...

        For now, lists, tuples, sets, Sets are coerced into finite
        enumerated sets::

            sage: S = EnumeratedSets()([1, 2, 3]); S
            {1, 2, 3}
            sage: S.category()
            Category of facade finite enumerated sets

            sage: S = EnumeratedSets()((1, 2, 3)); S
            {1, 2, 3}
            sage: S = EnumeratedSets()(set([1, 2, 3])); S
            {1, 2, 3}
            sage: S = EnumeratedSets()(Set([1, 2, 3])); S
            {1, 2, 3}
            sage: S.category()
            Category of facade finite enumerated sets
        """
        import sage.sets.set
        if isinstance(X, (tuple, list, set, sage.sets.set.Set_object_enumerated)):
            return sage.sets.all.FiniteEnumeratedSet(X)
        raise NotImplementedError

    class ParentMethods:

        def __iter__(self):
            """
            An iterator for the enumerated set.

            ``iter(self)`` allows the combinatorial class to be treated as an
            iterable. This if the default implementation from the category
            ``EnumeratedSets()`` it just goes through the iterator of the set
            to count the number of objects.

            By decreasing order of priority, the second column of the
            following array shows which methods is used to define
            ``__iter__``, when the methods of the first column are overloaded:

            +------------------------+---------------------------------+
            | Needed methods         | Default ``__iterator`` provided |
            +========================+=================================+
            | ``first`` and ``next`` | ``_iterator_from_next``         |
            +------------------------+---------------------------------+
            | ``unrank``             | ``_iterator_from_unrank``       |
            +------------------------+---------------------------------+
            | ``list`                | ``_iterator_from_next``         |
            +------------------------+---------------------------------+

            If non of these are provided raise a ``NotImplementedError``

            EXAMPLES::

            We start with an example where nothing is implemented::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: it = iter(broken()); [next(it), next(it), next(it)]
                Traceback (most recent call last):
                ...
                NotImplementedError: iterator called but not implemented

            Here is what happends when ``first`` and ``next`` are implemeted::

                sage: class set_first_next(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...    def first(self):
                ...        return 0
                ...    def next(self, elt):
                ...        return elt+1
                ...
                sage: it = iter(set_first_next()); [next(it), next(it), next(it)]
                [0, 1, 2]

            Let us try with ``unrank``::

                sage: class set_unrank(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...    def unrank(self, i):
                ...        return i + 5
                ...
                sage: it = iter(set_unrank()); [next(it), next(it), next(it)]
                [5, 6, 7]

            Let us finally try with ``list``::

                sage: class set_list(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...    def list(self):
                ...        return [5, 6, 7]
                ...
                sage: it = iter(set_list()); [next(it), next(it), next(it)]
                [5, 6, 7]

            """
            #Check to see if .first() and .next() are overridden in the subclass
            if ( self.first != self._first_from_iterator and
                 self.next  != self._next_from_iterator ):
                return self._iterator_from_next()
            #Check to see if .unrank() is overridden in the subclass
            elif self.unrank != self._unrank_from_iterator:
                return self._iterator_from_unrank()
            #Finally, check to see if .list() is overridden in the subclass
            elif self.list != self._list_default:
                return self._iterator_from_list()
            else:
                raise NotImplementedError("iterator called but not implemented")

        def is_empty(self):
            r"""
            Return whether this set is empty.

            EXAMPLES::

                sage: F = FiniteEnumeratedSet([1,2,3])
                sage: F.is_empty()
                False
                sage: F = FiniteEnumeratedSet([])
                sage: F.is_empty()
                True

            TESTS::

                sage: F.is_empty.__module__
                'sage.categories.enumerated_sets'
            """
            try:
                next(iter(self))
            except StopIteration:
                return True
            else:
                return False

        def list(self):
            """
            Returns an error since the cardinality of self is not known.

            EXAMPLES::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: broken().list()
                Traceback (most recent call last):
                ...
                NotImplementedError: unknown cardinality
            """
            raise NotImplementedError("unknown cardinality")
        _list_default  = list # needed by the check system.

        def _first_from_iterator(self):
            """
            The "first" element of ``self``.

            ``self.first()`` returns the first element of the set
            ``self``. This is a generic implementation from the category
            ``EnumeratedSets()`` which can be used when the method ``__iter__`` is
            provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.first() # indirect doctest
                1
            """
            return next(iter(self))
        first = _first_from_iterator

        def _next_from_iterator(self, obj):
            """
            The "next" element after ``obj`` in ``self``.

            ``self.next(e)`` returns the element of the set ``self`` which
            follows ``e``. This is a generic implementation from the category
            ``EnumeratedSets()`` which can be used when the method ``__iter__``
            is provided.

            Remark: this is the default (brute force) implementation
            of the category ``EnumeratedSets()``. Its complexity is
            `O(r)`, where `r` is the rank of ``obj``.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: C._next_from_iterator(10) # indirect doctest
                11

            TODO: specify the behavior when ``obj`` is not in ``self``.
            """
            it = iter(self)
            el = next(it)
            while el != obj:
                el = next(it)
            return next(it)
        next = _next_from_iterator


        def _unrank_from_iterator(self, r):
            """
            The ``r``-th element of ``self``

            ``self.unrank(r)`` returns the ``r``-th element of self, where
            ``r`` is an integer between ``0`` and ``n-1`` where ``n`` is the
            cardinality of ``self``.

            This is the default (brute force) implementation from the
            category ``EnumeratedSets()`` which can be used when the
            method ``__iter__`` is provided. Its complexity is `O(r)`,
            where `r` is the rank of ``obj``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.unrank(2) # indirect doctest
                3
                sage: C._unrank_from_iterator(5)
                Traceback (most recent call last):
                ...
                ValueError: the value must be between 0 and 2 inclusive
            """
            counter = 0
            for u in self:
                if counter == r:
                    return u
                counter += 1
            raise ValueError("the value must be between %s and %s inclusive"%(0,counter-1))
        unrank = _unrank_from_iterator

        def _rank_from_iterator(self, x):
            """
            The rank of an element of ``self``

            ``self.rank(x)`` returns the rank of `x`, that is its
            position in the enumeration of ``self``. This is an
            integer between ``0`` and ``n-1`` where ``n`` is the
            cardinality of ``self``, or None if `x` is not in `self`.

            This is the default (brute force) implementation from the
            category ``EnumeratedSets()`` which can be used when the
            method ``__iter__`` is provided. Its complexity is `O(r)`,
            where `r` is the rank of ``obj``. For infinite enumerated
            sets, this won't terminate when `x` is not in ``self``

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: list(C)
                [1, 2, 3]
                sage: C.rank(3) # indirect doctest
                2
                sage: C.rank(5) # indirect doctest
            """
            counter = 0
            for u in self:
                if u == x:
                    return counter
                counter += 1
            return None

        rank = _rank_from_iterator

        def _iterator_from_list(self):
            """
            An iterator for the elements of ``self``.

            ``iter(self)`` returns an iterator for the elements
            of ``self``. This is a generic implementation from the
            category ``EnumeratedSets()`` which can be used when the
            method ``list`` is provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: it = C._iterator_from_list()
                sage: [next(it), next(it), next(it)]
                [1, 2, 3]
            """
            for x in self.list():
                yield x

        def _iterator_from_next(self):
            """
            An iterator for the elements of ``self``.

            ``iter(self)`` returns an iterator for the element of
            the set ``self``. This is a generic implementation from
            the category ``EnumeratedSets()`` which can be used when
            the methods ``first`` and ``next`` are provided.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: it = C._iterator_from_next()
                sage: [next(it), next(it), next(it), next(it), next(it)]
                [0, 1, 2, 3, 4]
            """
            f = self.first()
            yield f
            while True:
                try:
                    f = self.next(f)
                except (TypeError, ValueError ):
                    break

                if f is None or f is False :
                    break
                else:
                    yield f

        def _iterator_from_unrank(self):
            """
            An iterator for the elements of ``self``.

            ``iter(self)`` returns an iterator for the elements
            of the set ``self``. This is a generic implementation from
            the category ``EnumeratedSets()`` which can be used when
            the method ``unrank`` is provided.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: it = C._iterator_from_unrank()
                sage: [next(it), next(it), next(it), next(it), next(it)]
                [0, 1, 2, 3, 4]
            """
            r = 0
            try:
                u = self.unrank(r)
            except (TypeError, ValueError, IndexError):
                return
            yield u
            while True:
                r += 1
                try:
                    u = self.unrank(r)
                except (TypeError, ValueError, IndexError):
                    break

                if u is None:
                    break
                else:
                    yield u

        # This @cached_method is not really needed, since the method
        # an_element itself is cached. We leave it for the moment, so
        # that Parents that do not yet inherit properly from categories
        # (e.g. Set([1,2,3]) can use the following trick:
        #    _an_element_ = EnumeratedSets.ParentMethods._an_element_
        @cached_method
        def _an_element_from_iterator(self):
            """
            Returns the first element of ``self`` returned by :meth:`__iter__`

            If ``self`` is empty, the exception
            :class:`~sage.categories.sets_cat.EmptySetError` is raised instead.

            This provides a generic implementation of the method
            :meth:`_an_element_` for all parents in :class:`EnumeratedSets`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example(); C
                An example of a finite enumerated set: {1,2,3}
                sage: C.an_element() # indirect doctest
                1
                sage: S = Set([])
                sage: S.an_element()
                Traceback (most recent call last):
                ...
                EmptySetError

            TESTS::

                sage: super(Parent, C)._an_element_
                Cached version of <function _an_element_from_iterator at ...>
            """
            it = self.__iter__()
            try:
                return next(it)
            except StopIteration:
                raise EmptySetError

        # Should this be implemented from first instead?
        _an_element_ = _an_element_from_iterator

        #FIXME: use combinatorial_class_from_iterator once class_from_iterator.patch is in
        def _some_elements_from_iterator(self):
            """
            Returns some elements in ``self``.

            See :class:`TestSuite` for a typical use case.

            This is a generic implementation from the category
            ``EnumeratedSets()`` which can be used when the method
            ``__iter__`` is provided. It returns an iterator for up to
            the first 100 elements of ``self``

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: list(C.some_elements()) # indirect doctest
                [1, 2, 3]
            """
            nb = 0
            for i in self:
                yield i
                nb += 1
                if nb >= 100:
                    break
        some_elements = _some_elements_from_iterator

        def random_element(self):
            """
            Returns a random element in ``self``.

            Unless otherwise stated, and for finite enumerated sets,
            the probability is uniform.

            This is a generic implementation from the category
            ``EnumeratedSets()``. It raise a ``NotImplementedError``
            since one does not know whether the set is finite.

            EXAMPLES::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: broken().random_element()
                Traceback (most recent call last):
                ...
                NotImplementedError: unknown cardinality
                """
            raise NotImplementedError("unknown cardinality")

        def map(self, f, name=None):
            r"""
            Returns the image `\{f(x) | x \in \text{self}\}` of this
            enumerated set by `f`, as an enumerated set.

            `f` is supposed to be injective.

            EXAMPLES::

                sage: R = SymmetricGroup(3).map(attrcall('reduced_word')); R
                Image of Symmetric group of order 3! as a permutation group by *.reduced_word()
                sage: R.cardinality()
                6
                sage: R.list()
                [[], [1], [2, 1], [1, 2], [2], [1, 2, 1]]
                sage: [ r for r in R]
                [[], [1], [2, 1], [1, 2], [2], [1, 2, 1]]

            .. warning::

                If the function is not injective, then there may be
                repeated elements::

                    sage: P = SymmetricGroup(3)
                    sage: P.list()
                    [(), (1,2), (1,2,3), (1,3,2), (2,3), (1,3)]
                    sage: P.map(attrcall('length')).list()
                    [0, 1, 2, 2, 1, 3]

            .. warning::

                :class:`MapCombinatorialClass` needs to be refactored to use categories::

                    sage: R.category()             # todo: not implemented
                    Category of enumerated sets
                    sage: TestSuite(R).run(skip=['_test_an_element', '_test_category', '_test_some_elements'])
            """
            from sage.combinat.combinat import MapCombinatorialClass
            return MapCombinatorialClass(self, f, name)

#
#  Consistency test suite for an enumerated set:
#
        def _test_enumerated_set_contains(self, **options):
            """
            Checks that the methods :meth:`.__contains__` and :meth:`.__iter__` are consistent.

            See also :class:`TestSuite`.

            TESTS::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_contains()
                sage: TestSuite(C).run()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ...       def __contains__(self, obj):
                ...           if obj == 3:
                ...               return False
                ...           else:
                ...               return obj in C
                sage: CC = CCls()
                sage: CC._test_enumerated_set_contains()
                Traceback (most recent call last):
                ...
                AssertionError: False is not true
            """
            tester = self._tester(**options)
            i = 0
            for w in self:
                tester.assertTrue(w in self)
                i += 1
                if i > tester._max_runs:
                    return

        def _test_enumerated_set_iter_list(self, **options):
            """
            Checks that the methods :meth:`.list` and :meth:`.__iter__` are consistent.

            See also: :class:`TestSuite`.

            .. NOTE::

                This test does nothing if the cardinality of the set
                is larger than the max_runs argument.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_iter_list()
                sage: TestSuite(C).run()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ...       def list(self):
                ...           return [1,2,3,4]
                sage: CC = CCls()
                sage: CC._test_enumerated_set_iter_list()
                Traceback (most recent call last):
                ...
                AssertionError: 3 != 4

            For a large enumerated set this test does nothing:
            increase tester._max_runs if you want to actually run the
            test::

                sage: class CCls(Example):
                ...       def list(self):
                ...           return [1,2,3]
                sage: CC = CCls()
                sage: CC._test_enumerated_set_iter_list(verbose=True,max_runs=2)
                Enumerated set too big; skipping test; increase tester._max_runs
            """
            tester = self._tester(**options)
            if self.list != self._list_default:
                # TODO: if self._cardinality is self._cardinality_from_iterator
                # we could make sure to stop the counting at
                # self.max_test_enumerated_set_loop
                if self.cardinality() > tester._max_runs:
                    tester.info("Enumerated set too big; skipping test; increase tester._max_runs")
                    return
                ls = self.list()
                i = 0
                for obj in self:
                    tester.assertEqual(obj, ls[i])
                    i += 1
                tester.assertEqual(i, len(ls))

    class ElementMethods:

        def rank(self):
            """
            Returns the rank of ``self`` in its parent.

            See also :meth:`EnumeratedSets.ElementMethods.rank`

            EXAMPLES::

                sage: F = FiniteSemigroups().example(('a','b','c'))
                sage: L = list(F); L
                ['a', 'c', 'ac', 'b', 'ba', 'bc', 'cb', 'ca', 'bca', 'ab', 'bac', 'cab', 'acb', 'cba', 'abc']
                sage: L[7].rank()
                7
            """
            return self.parent().rank(self)

    Finite   = LazyImport('sage.categories.finite_enumerated_sets', 'FiniteEnumeratedSets', at_startup=True)
    Infinite = LazyImport('sage.categories.infinite_enumerated_sets', 'InfiniteEnumeratedSets', at_startup=True)

    class CartesianProducts(CartesianProductsCategory):

        class ParentMethods:
            def __iter__(self):
                r"""
                Return a lexicographic iterator for the elements of this cartesian product.

                EXAMPLES::

                    sage: A = FiniteEnumeratedSets()(["a", "b"])
                    sage: B = FiniteEnumeratedSets().example(); B
                    An example of a finite enumerated set: {1,2,3}
                    sage: C = cartesian_product([A, B, A]); C
                    The cartesian product of ({'a', 'b'}, An example of a finite enumerated set: {1,2,3}, {'a', 'b'})
                    sage: C in FiniteEnumeratedSets()
                    True
                    sage: list(C)
                    [('a', 1, 'a'), ('a', 1, 'b'), ('a', 2, 'a'), ('a', 2, 'b'), ('a', 3, 'a'), ('a', 3, 'b'),
                     ('b', 1, 'a'), ('b', 1, 'b'), ('b', 2, 'a'), ('b', 2, 'b'), ('b', 3, 'a'), ('b', 3, 'b')]
                    sage: C.__iter__.__module__
                    'sage.categories.enumerated_sets'

                    sage: F22 = GF(2).cartesian_product(GF(2))
                    sage: list(F22)
                    [(0, 0), (0, 1), (1, 0), (1, 1)]

                    sage: C = cartesian_product([Permutations(10)]*4)
                    sage: it = iter(C)
                    sage: it.next()
                    ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
                    sage: it.next()
                    ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                     [1, 2, 3, 4, 5, 6, 7, 8, 10, 9])

                .. WARNING::

                    The elements are returned in lexicographic order,
                    which gives a valid enumeration only if all
                    factors, but possibly the first one, are
                    finite. So the following one is fine::

                        sage: it = iter(cartesian_product([ZZ, GF(2)]))
                        sage: [it.next() for _ in range(10)]
                        [(0, 0), (0, 1), (1, 0), (1, 1),
                         (-1, 0), (-1, 1), (2, 0), (2, 1),
                         (-2, 0), (-2, 1)]

                    But this one is not::

                        sage: it = iter(cartesian_product([GF(2), ZZ]))
                        sage: [it.next() for _ in range(10)]
                        doctest:...: UserWarning: Sage is not able to determine
                        whether the factors of this cartesian product are
                        finite. The lexicographic ordering might not go through
                        all elements.
                        [(0, 0), (0, 1), (0, -1), (0, 2), (0, -2),
                         (0, 3), (0, -3), (0, 4), (0, -4), (0, 5)]

                .. NOTE::

                    Here it would be faster to use :func:`itertools.product` for sets
                    of small size. But the latter expands all factor in memory!
                    So we can not reasonably use it in general.

                ALGORITHM:

                Recipe 19.9 in the Python Cookbook by Alex Martelli
                and David Ascher.
                """
                if any(f not in Sets().Finite() for f in self.cartesian_factors()[1:]):
                    from warnings import warn
                    warn("Sage is not able to determine whether the factors of "
                         "this cartesian product are finite. The lexicographic "
                         "ordering might not go through all elements.")

                # visualize an odometer, with "wheels" displaying "digits"...:
                factors = list(self.cartesian_factors())
                wheels = map(iter, factors)
                digits = [next(it) for it in wheels]
                while True:
                    yield self._cartesian_product_of_elements(digits)
                    for i in range(len(digits)-1, -1, -1):
                        try:
                            digits[i] = next(wheels[i])
                            break
                        except StopIteration:
                            wheels[i] = iter(factors[i])
                            digits[i] = next(wheels[i])
                    else:
                        break

            def first(self):
                r"""
                Return the first element.

                EXAMPLES::

                    sage: cartesian_product([ZZ]*10).first()
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                """
                return self._cartesian_product_of_elements(
                        tuple(c.first() for c in self.cartesian_factors()))
