# sage_setup: distribution = sagemath-categories
r"""
Finite Enumerated Sets
"""
# ****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.sets_cat import Sets
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.isomorphic_objects import IsomorphicObjectsCategory
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.cpython.getattr import raw_getattr
lazy_import("sage.rings.integer", "Integer")


class FiniteEnumeratedSets(CategoryWithAxiom):
    """
    The category of finite enumerated sets.

    EXAMPLES::

        sage: FiniteEnumeratedSets()
        Category of finite enumerated sets
        sage: FiniteEnumeratedSets().super_categories()
        [Category of enumerated sets, Category of finite sets]
        sage: FiniteEnumeratedSets().all_super_categories()
        [Category of finite enumerated sets,
         Category of enumerated sets,
         Category of finite sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    TESTS::

        sage: C = FiniteEnumeratedSets()
        sage: TestSuite(C).run()
        sage: sorted(C.Algebras(QQ).super_categories(), key=str)
        [Category of finite dimensional vector spaces with basis over Rational Field,
         Category of set algebras over Rational Field]

    .. TODO::

        :class:`sage.combinat.debruijn_sequence.DeBruijnSequences` should
        not inherit from this class. If that is solved, then
        :class:`FiniteEnumeratedSets` shall be turned into a subclass of
        :class:`~sage.categories.category_singleton.Category_singleton`.
    """

    def _call_(self, X):
        """
        Construct an object in this category from the data in ``X``.

        EXAMPLES::

            sage: FiniteEnumeratedSets()(GF(3))
            Finite Field of size 3
            sage: Partitions(3)                                                         # needs sage.combinat
            Partitions of the integer 3

        For now, lists, tuples, sets, Sets are coerced into finite
        enumerated sets::

            sage: FiniteEnumeratedSets()([1, 2, 3])
            {1, 2, 3}
            sage: FiniteEnumeratedSets()((1, 2, 3))
            {1, 2, 3}
            sage: FiniteEnumeratedSets()(set([1, 2, 3]))
            {1, 2, 3}
            sage: FiniteEnumeratedSets()(Set([1, 2, 3]))
            {1, 2, 3}
        """
        return EnumeratedSets()._call_(X)

    class ParentMethods:

        def __len__(self):
            """
            Return the number of elements of ``self``.

            EXAMPLES::

                sage: len(GF(5))
                5
                sage: len(MatrixSpace(GF(2), 3, 3))                                     # needs sage.modules
                512
            """
            return int(self.cardinality())

        def _cardinality_from_iterator(self, *ignored_args, **ignored_kwds):
            """
            Return the cardinality of ``self``.

            This brute force implementation of :meth:`cardinality`
            iterates through the elements of ``self`` to count them.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example(); C
                An example of a finite enumerated set: {1,2,3}
                sage: C._cardinality_from_iterator()
                3

            TESTS:

            This is the default implementation of :meth:`cardinality`
            from the category ``FiniteEnumeratedSet()``. To test this,
            we need a fresh example::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class FreshExample(Example): pass
                sage: C = FreshExample(); C.rename('FreshExample')
                sage: C.cardinality
                <bound method FiniteEnumeratedSets.ParentMethods._cardinality_from_iterator of FreshExample>

            This method shall return an ``Integer``; we test this
            here, because :meth:`_test_enumerated_set_iter_cardinality`
            does not do it for us::

                sage: type(C._cardinality_from_iterator())
                <class 'sage.rings.integer.Integer'>

            We ignore additional inputs since during doctests classes which
            override ``cardinality()`` call up to the category rather than
            their own ``cardinality()`` method (see :issue:`13688`)::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._cardinality_from_iterator(algorithm='testing')
                3

            Here is a more complete example::

                sage: class TestParent(Parent):
                ....:   def __init__(self):
                ....:       Parent.__init__(self, category=FiniteEnumeratedSets())
                ....:   def __iter__(self):
                ....:       yield 1
                ....:       return
                ....:   def cardinality(self, dummy_arg):
                ....:       return 1 # we don't want to change the semantics of cardinality()
                sage: P = TestParent()
                sage: P.cardinality(-1)
                1
                sage: v = P.list(); v
                [1]
                sage: P.cardinality()
                1
                sage: P.cardinality('use alt algorithm') # Used to break here: see trac #13688
                1
                sage: P.cardinality(dummy_arg='use alg algorithm') # Used to break here: see trac #13688
                1
            """
            c = 0
            for _ in self:
                c += 1
            return Integer(c)

        #Set cardinality to the default implementation
        cardinality = _cardinality_from_iterator

        def _cardinality_from_list(self, *ignored_args, **ignored_kwds):
            """
            The cardinality of ``self``.

            This implementation of :meth:`cardinality` computes the
            cardinality from :meth:`list` (which is
            cached). Reciprocally, calling ``self.list()`` makes this
            method the default implementation of :meth:`cardinality`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._cardinality_from_list()
                3

            We ignore additional inputs since during doctests classes which
            override ``cardinality()`` call up to the category rather than
            their own ``cardinality()`` method (see :issue:`13688`)::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._cardinality_from_list(algorithm='testing')
                3
            """
            return Integer(len(self.tuple()))

        def _unrank_from_list(self, r):
            r"""
            The `r`-th element of ``self``.

            INPUT:

              - ``r`` -- integer between `0` and `n-1`,
                where `n` is the cardinality of ``self``

            OUTPUT: the `r`-th element of ``self``

            This implementation of :meth:`unrank` uses the method
            :meth:`list` (which is cached). Reciprocally, calling
            ``self.list()`` makes this method the default
            implementation of :meth:`unrank`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._unrank_from_list(1)
                2
            """
            lst = self.tuple()
            try:
                return lst[r]
            except IndexError:
                raise ValueError("the value must be in the range from %s to %s" % (0, len(lst) - 1))

        def tuple(self):
            r"""
            Return a :class:`tuple`of the elements of ``self``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.tuple()
                (1, 2, 3)
                sage: C.tuple() is C.tuple()
                True
            """
            # Simpler implementation because it does not have to check whether cardinality is finite
            try: # shortcut
                if self._list is not None:
                    return self._tuple_from_list()
            except AttributeError:
                pass

            if self.list != self._list_default:
                return tuple(self.list())

            return self._tuple_from_iterator()
        _tuple_default = tuple

        def _list_from_iterator(self):
            r"""
            Return a list of the elements of ``self`` after cached.

            It moreover overrides the following methods to use this cache:

            - ``self.__iter__()``
            - ``self.cardinality()``
            - ``self.unrank()``

            .. SEEALSO:: :meth:`_cardinality_from_list`,
                :meth:`_iterator_from_list`, and :meth:`_unrank_from_list`

            .. WARNING::

                The overriding of ``self.__iter__`` to use the cache
                is ignored upon calls such as ``for x in C:`` or
                ``list(C)`` (which essentially ruins its purpose).
                Indeed, Python looks up the ``__iter__`` method
                directly in the class of ``C``, bypassing ``C``'s
                dictionary (see the Python reference manual,
                `Special method lookup for new-style classes <http://docs.python.org/reference/datamodel.html#special-method-lookup-for-new-style-classes>`_)

                Let's take an example::

                    sage: class Example(Parent):
                    ....:     def __init__(self):
                    ....:         Parent.__init__(self, category = FiniteEnumeratedSets())
                    ....:     def __iter__(self):
                    ....:         print("hello!")
                    ....:         for x in [1,2,3]: yield x
                    sage: C = Example()
                    sage: list(C)
                    hello!
                    ...
                    [1, 2, 3]
                    sage: list(C)
                    hello!
                    ...
                    [1, 2, 3]

                Note that ``hello!`` actually gets printed more than once in
                the calls to ``list(C)``. That's because of the
                implicit calls to :meth:`__len__`, which also relies
                on :meth:`__iter__`. Let's call :meth:`list`::

                    sage: C.list()
                    hello!
                    [1, 2, 3]
                    sage: C.list()
                    [1, 2, 3]

                Now we would want the original iterator of ``C`` not
                to be called anymore, but that's not the case::

                    sage: list(C)
                    hello!
                    [1, 2, 3]

            TESTS:

            To test if the caching and overriding works, we need a
            fresh finite enumerated set example, because the caching
            mechanism has already been triggered::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class FreshExample(Example): pass
                sage: C = FreshExample(); C.rename('FreshExample')
                sage: C.list
                <bound method EnumeratedSets.ParentMethods.list of FreshExample>
                sage: C.unrank
                <bound method EnumeratedSets.ParentMethods._unrank_from_iterator of FreshExample>
                sage: C.cardinality
                <bound method FiniteEnumeratedSets.ParentMethods._cardinality_from_iterator of FreshExample>
                sage: l1 = C.list(); l1
                [1, 2, 3]
                sage: C.list
                <bound method EnumeratedSets.ParentMethods.list of FreshExample>
                sage: C.unrank
                <bound method FiniteEnumeratedSets.ParentMethods._unrank_from_list of FreshExample>
                sage: C.cardinality
                <bound method FiniteEnumeratedSets.ParentMethods._cardinality_from_list of FreshExample>
                sage: C.__iter__
                <bound method EnumeratedSets.ParentMethods._iterator_from_list of FreshExample>

            We finally check that nothing breaks before and after
            calling explicitly the method ``.list()``::

                sage: class FreshExample(Example): pass
                sage: import __main__; __main__.FreshExample = FreshExample # Fake FreshExample being defined in a python module
                sage: C = FreshExample()
                sage: TestSuite(C).run()
                sage: C.list()
                [1, 2, 3]
                sage: TestSuite(C).run()
            """
            try:
                if self._list is not None:
                    return list(self._list)
            except AttributeError:
                pass
            result = tuple(self.__iter__())
            try:
                self._list = result
                self.__iter__ = self._iterator_from_list
                self.cardinality = self._cardinality_from_list
                self.tuple = self._tuple_from_list
                self.unrank = self._unrank_from_list
            except AttributeError:
                pass
            return list(result)

        def unrank_range(self, start=None, stop=None, step=None):
            r"""
            Return the range of elements of ``self`` starting at ``start``,
            ending at ``stop``, and stepping by ``step``.

            See also ``unrank()``.

            EXAMPLES::

                sage: F = FiniteEnumeratedSet([1,2,3])
                sage: F.unrank_range(1)
                [2, 3]
                sage: F.unrank_range(stop=2)
                [1, 2]
                sage: F.unrank_range(stop=2, step=2)
                [1]
                sage: F.unrank_range(start=1, step=2)
                [2]
                sage: F.unrank_range(stop=-1)
                [1, 2]

                sage: F = FiniteEnumeratedSet([1,2,3,4])
                sage: F.unrank_range(stop=10)
                [1, 2, 3, 4]
            """
            try:
                return list(self._list[start:stop:step])
            except AttributeError:
                pass
            card = self.cardinality() # This may set the list
            try:
                return list(self._list[start:stop:step])
            except AttributeError:
                pass
            if start is None and stop is not None and stop >= 0 and step is None:
                if stop < card:
                    it = self.__iter__()
                    return [next(it) for j in range(stop)]
                return self.list()
            return list(self.tuple()[start:stop:step])

        def iterator_range(self, start=None, stop=None, step=None):
            r"""
            Iterate over the range of elements of ``self`` starting
            at ``start``, ending at ``stop``, and stepping by ``step``.

            .. SEEALSO::

                ``unrank()``, ``unrank_range()``

            EXAMPLES::

                sage: F = FiniteEnumeratedSet([1,2,3])
                sage: list(F.iterator_range(1))
                [2, 3]
                sage: list(F.iterator_range(stop=2))
                [1, 2]
                sage: list(F.iterator_range(stop=2, step=2))
                [1]
                sage: list(F.iterator_range(start=1, step=2))
                [2]
                sage: list(F.iterator_range(start=1, stop=2))
                [2]
                sage: list(F.iterator_range(start=0, stop=1))
                [1]
                sage: list(F.iterator_range(start=0, stop=3, step=2))
                [1, 3]
                sage: list(F.iterator_range(stop=-1))
                [1, 2]

                sage: F = FiniteEnumeratedSet([1,2,3,4])
                sage: list(F.iterator_range(start=1, stop=3))
                [2, 3]
                sage: list(F.iterator_range(stop=10))
                [1, 2, 3, 4]
            """
            L = None
            try:
                L = self._list
            except AttributeError:
                pass
            card = self.cardinality() # This may set the list
            try:
                L = self._list
            except AttributeError:
                pass
            if L is None and start is None and stop is not None and stop >= 0 and step is None:
                if stop < card:
                    it = self.__iter__()
                    for j in range(stop):
                        yield next(it)
                    return
                yield from self
                return
            if L is None:
                L = self.tuple()
            yield from L[start:stop:step]

        def _random_element_from_unrank(self):
            """
            A random element in ``self``.

            ``self.random_element()`` returns a random element in
            ``self`` with uniform probability.

            This is the default implementation from the category
            ``EnumeratedSet()`` which uses the method ``unrank``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: n = C.random_element()
                sage: n in C
                True

                sage: n = C._random_element_from_unrank()
                sage: n in C
                True

            TODO: implement _test_random which checks uniformness
            """
            from sage.misc.prandom import randint
            c = self.cardinality()
            r = randint(0, c-1)
            return self.unrank(r)
        # Set the default implementation of random_element
        random_element = _random_element_from_unrank

        @cached_method
        def _last_from_iterator(self):
            """
            The last element of ``self``.

            ``self.last()`` returns the last element of ``self``.

            This is the default (brute force) implementation from the
            category ``FiniteEnumeratedSet()`` which can be used when
            the method ``__iter__`` is provided. Its complexity is
            `O(n)` where `n` is the size of ``self``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.last()
                3
                sage: C._last_from_iterator()
                3
            """
            for i in self:
                pass
            return i
        last = _last_from_iterator

        def _last_from_unrank(self):
            """
            The last element of ``self``.

            ``self.last()`` returns the last element of ``self``

            This is a generic implementation from the category
            ``FiniteEnumeratedSet()`` which can be used when the
            method ``unrank`` is provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._last_from_unrank()
                3
            """
            return self.unrank(self.cardinality() - 1)

        def _test_enumerated_set_iter_cardinality(self, **options):
            """
            Check that the methods :meth:`.cardinality` and
            :meth:`.__iter__` are consistent. Also checks that
            :meth:`.cardinality` returns an ``Integer``.

            For efficiency reasons, those tests are not run if
            :meth:`.cardinality` is
            :meth:`._cardinality_from_iterator`, or if ``self`` is too
            big.

            .. SEEALSO:: :class:`TestSuite`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_iter_cardinality()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ....:     def cardinality(self):
                ....:         return 4
                sage: CC = CCls()
                sage: CC._test_enumerated_set_iter_cardinality()
                Traceback (most recent call last):
                ...
                AssertionError: 4 != 3
            """
            tester = self._tester(**options)
            if self.cardinality != self._cardinality_from_iterator:
                card = self.cardinality()
                if card <= tester._max_runs:
                    tester.assertEqual(card,
                                       self._cardinality_from_iterator())

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            A Cartesian product of finite enumerated sets is a finite
            enumerated set.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of finite enumerated sets]
            """
            return [FiniteEnumeratedSets()]

        class ParentMethods:
            r"""
            TESTS:

            Ideally, these tests should be just after the declaration of the
            associated attributes. But doing this way, Sage will not consider
            them as a doctest.

            We check that Cartesian products of finite enumerated sets
            inherit various methods from `Sets.CartesianProducts`
            and not from :class:`EnumeratedSets.Finite`::

                sage: C = cartesian_product([Partitions(10), Permutations(20)])         # needs sage.combinat
                sage: C in EnumeratedSets().Finite()                                    # needs sage.combinat
                True

                sage: C.random_element.__module__                                       # needs sage.combinat
                'sage.categories.sets_cat'

                sage: C.cardinality.__module__                                          # needs sage.combinat
                'sage.categories.sets_cat'

                sage: C.__iter__.__module__                                             # needs sage.combinat
                'sage.categories.sets_cat'
            """
            random_element = raw_getattr(Sets.CartesianProducts.ParentMethods, "random_element")
            cardinality = raw_getattr(Sets.CartesianProducts.ParentMethods, "cardinality")
            __iter__ = raw_getattr(Sets.CartesianProducts.ParentMethods, "__iter__")

            def last(self):
                r"""
                Return the last element.

                EXAMPLES::

                    sage: C = cartesian_product([Zmod(42), Partitions(10),              # needs sage.combinat
                    ....:                        IntegerRange(5)])
                    sage: C.last()                                                      # needs sage.combinat
                    (41, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 4)
                """
                return self._cartesian_product_of_elements(
                        tuple(c.last() for c in self.cartesian_factors()))

            def rank(self, x):
                r"""
                Return the rank of an element of this Cartesian product.

                The *rank* of ``x`` is its position in the
                enumeration. It is an integer between ``0`` and
                ``n-1`` where ``n`` is the cardinality of this set.

                .. SEEALSO::

                    - :meth:`EnumeratedSets.ParentMethods.rank`
                    - :meth:`unrank`

                EXAMPLES::

                    sage: C = cartesian_product([GF(2), GF(11), GF(7)])
                    sage: C.rank(C((1,2,5)))
                    96
                    sage: C.rank(C((0,0,0)))
                    0

                    sage: for c in C: print(C.rank(c))
                    0
                    1
                    2
                    3
                    4
                    5
                    ...
                    150
                    151
                    152
                    153

                    sage: # needs sage.combinat
                    sage: F1 = FiniteEnumeratedSet('abcdefgh')
                    sage: F2 = IntegerRange(250)
                    sage: F3 = Partitions(20)
                    sage: C = cartesian_product([F1, F2, F3])
                    sage: c = C(('a', 86, [7,5,4,4]))
                    sage: C.rank(c)
                    54213
                    sage: C.unrank(54213)
                    ('a', 86, [7, 5, 4, 4])
                """
                from sage.rings.integer_ring import ZZ
                x = self(x)
                b = ZZ.one()
                rank = ZZ.zero()
                for f, c in zip(reversed(x.cartesian_factors()),
                                reversed(self.cartesian_factors())):
                    rank += b * c.rank(f)
                    b *= c.cardinality()
                return rank

            def unrank(self, i):
                r"""
                Return the `i`-th element of this Cartesian product.

                INPUT:

                - ``i`` -- integer between `0` and `n-1` where
                  `n` is the cardinality of this set

                .. SEEALSO::

                    - :meth:`EnumeratedSets.ParentMethods.unrank`
                    - :meth:`rank`

                EXAMPLES::

                    sage: C = cartesian_product([GF(3), GF(11), GF(7), GF(5)])
                    sage: c = C.unrank(123); c
                    (0, 3, 3, 3)
                    sage: C.rank(c)
                    123

                    sage: c = C.unrank(857); c
                    (2, 2, 3, 2)
                    sage: C.rank(c)
                    857

                    sage: C.unrank(2500)
                    Traceback (most recent call last):
                    ...
                    IndexError: index i (=2) is greater than the cardinality
                """
                from sage.rings.integer_ring import ZZ
                i = ZZ(i)
                if i < 0:
                    raise IndexError("i (={}) must be a nonnegative integer")
                elt = []
                for c in reversed(self.cartesian_factors()):
                    card = c.cardinality()
                    elt.insert(0, c.unrank(i % card))
                    i //= card
                if i:
                    raise IndexError("index i (={}) is greater than the cardinality".format(i))
                return self._cartesian_product_of_elements(elt)

    class IsomorphicObjects(IsomorphicObjectsCategory):

        def example(self):
            """
            Return an example of isomorphic object of a finite
            enumerated set, as per :meth:`Category.example
            <sage.categories.category.Category.example>`.

            EXAMPLES::

                sage: FiniteEnumeratedSets().IsomorphicObjects().example()
                The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            """
            from sage.categories.examples.finite_enumerated_sets import IsomorphicObjectOfFiniteEnumeratedSet
            return IsomorphicObjectOfFiniteEnumeratedSet()

        class ParentMethods:

            def cardinality(self):
                r"""
                Return the cardinality of ``self`` which is the same
                as that of the ambient set ``self`` is isomorphic to.

                EXAMPLES::

                    sage: A = FiniteEnumeratedSets().IsomorphicObjects().example(); A
                    The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
                    sage: A.cardinality()
                    3
                """
                return self.ambient().cardinality()

            def __iter__(self):
                r"""
                Return an iterator over ``self``, using the bijection
                with the ambient space.

                EXAMPLES::

                    sage: A = FiniteEnumeratedSets().IsomorphicObjects().example(); A
                    The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
                    sage: list(A)                  # indirect doctest
                    [1, 4, 9]
                """
                for x in self.ambient():
                    yield self.retract(x)
