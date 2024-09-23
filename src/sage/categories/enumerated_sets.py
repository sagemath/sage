# sage_setup: distribution = sagemath-categories
r"""
Enumerated sets
"""
# ****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.sets_cat import Sets
from sage.categories.sets_cat import EmptySetError
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.misc.lazy_import import lazy_import
lazy_import("sage.rings.integer", "Integer")


class EnumeratedSets(CategoryWithAxiom):
    """
    The category of enumerated sets.

    An *enumerated set* is a *finite* or *countable* set or multiset `S`
    together with a canonical enumeration of its elements;
    conceptually, this is very similar to an immutable list. The main
    difference lies in the names and the return type of the methods,
    and of course the fact that the list of elements is not supposed to
    be expanded in memory. Whenever possible one should use one of the
    two sub-categories :class:`FiniteEnumeratedSets` or
    :class:`InfiniteEnumeratedSets`.

    The purpose of this category is threefold:

     - to fix a common interface for all these sets;
     - to provide a bunch of default implementations;
     - to provide consistency tests.

    The standard methods for an enumerated set ``S`` are:

       - ``S.cardinality()`` -- the number of elements of the set. This
         is the equivalent for ``len`` on a list except that the
         return value is specified to be a Sage :class:`Integer` or
         ``infinity``, instead of a Python ``int``.

       - ``iter(S)`` -- an iterator for the elements of the set;

       - ``S.list()`` -- a fresh list of the elements of the set, when
         possible; raises a :exc:`NotImplementedError` if the list is
         predictably too large to be expanded in memory.

       - ``S.tuple()`` -- a tuple of the elements of the set, when
         possible; raises a :exc:`NotImplementedError` if the tuple is
         predictably too large to be expanded in memory.

       - ``S.unrank(n)`` -- the  ``n``-th element of the set when ``n`` is a sage
         ``Integer``. This is the equivalent for ``l[n]`` on a list.

       - ``S.rank(e)`` -- the position of the element ``e`` in the set;
         This is equivalent to ``l.index(e)`` for a list except that
         the return value is specified to be a Sage :class:`Integer`,
         instead of a Python ``int``.

       - ``S.first()`` -- the first object of the set; it is equivalent to
         ``S.unrank(0)``.

       - ``S.next(e)`` -- the object of the set which follows ``e``; it is
         equivalent to ``S.unrank(S.rank(e) + 1)``.

       - ``S.random_element()`` -- a random generator for an element of
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
        [Category of enumerated sets, Category of sets,
         Category of sets with partial maps, Category of objects]

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
            Category of finite enumerated sets

        Also Python3 range are now accepted::

            sage: S = EnumeratedSets()(range(4)); S
            {0, 1, 2, 3}
        """
        import sage.sets.set
        if isinstance(X, (tuple, list, set, range, sage.sets.set.Set_object_enumerated)):
            return sage.sets.finite_enumerated_set.FiniteEnumeratedSet(X)
        raise NotImplementedError

    class ParentMethods:

        def __iter__(self):
            """
            An iterator for the enumerated set.

            ``iter(self)`` allows the combinatorial class to be treated as an
            iterable. This is the default implementation from the category
            ``EnumeratedSets()``; it just goes through the iterator of the set
            to count the number of objects.

            By decreasing order of priority, the second column of the
            following array shows which method is used to define
            ``__iter__``, when the methods of the first column are overloaded:

            +------------------------+---------------------------------+
            | Needed methods         | Default ``__iterator`` provided |
            +========================+=================================+
            | ``first`` and ``next`` | ``_iterator_from_next``         |
            +------------------------+---------------------------------+
            | ``unrank``             | ``_iterator_from_unrank``       |
            +------------------------+---------------------------------+
            | ``list``               | ``_iterator_from_next``         |
            +------------------------+---------------------------------+

            It is also possible to override ``__iter__`` method itself. Then
            the methods of the first column are defined using  ``__iter__``

            If none of these are provided, this raises
            a :exc:`NotImplementedError`.

            EXAMPLES:

            We start with an example where nothing is implemented::

                sage: class broken(UniqueRepresentation, Parent):
                ....:     def __init__(self):
                ....:         Parent.__init__(self, category = EnumeratedSets())
                sage: it = iter(broken()); [next(it), next(it), next(it)]
                Traceback (most recent call last):
                ...
                NotImplementedError: iterator called but not implemented

            Here is what happens when ``first`` and ``next`` are implemented::

                sage: class set_first_next(UniqueRepresentation, Parent):
                ....:     def __init__(self):
                ....:         Parent.__init__(self, category = EnumeratedSets())
                ....:     def first(self):
                ....:         return 0
                ....:     def next(self, elt):
                ....:         return elt+1
                sage: it = iter(set_first_next()); [next(it), next(it), next(it)]
                [0, 1, 2]

            Let us try with ``unrank``::

                sage: class set_unrank(UniqueRepresentation, Parent):
                ....:     def __init__(self):
                ....:         Parent.__init__(self, category = EnumeratedSets())
                ....:     def unrank(self, i):
                ....:         return i + 5
                sage: it = iter(set_unrank()); [next(it), next(it), next(it)]
                [5, 6, 7]

            Let us finally try with ``list``::

                sage: class set_list(UniqueRepresentation, Parent):
                ....:     def __init__(self):
                ....:         Parent.__init__(self, category = EnumeratedSets())
                ....:     def list(self):
                ....:         return [5, 6, 7]
                sage: it = iter(set_list()); [next(it), next(it), next(it)]
                [5, 6, 7]
            """
            # Check if .first() and .next(x) are overridden in the subclass
            if ( self.first != self._first_from_iterator and
                 self.next != self._next_from_iterator ):
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

        def iterator_range(self, start=None, stop=None, step=None):
            r"""
            Iterate over the range of elements of ``self`` starting
            at ``start``, ending at ``stop``, and stepping by ``step``.

            .. SEEALSO::

                ``unrank()``, ``unrank_range()``

            EXAMPLES::

                sage: # needs sage.combinat
                sage: P = Partitions()
                sage: list(P.iterator_range(stop=5))
                [[], [1], [2], [1, 1], [3]]
                sage: list(P.iterator_range(0, 5))
                [[], [1], [2], [1, 1], [3]]
                sage: list(P.iterator_range(3, 5))
                [[1, 1], [3]]
                sage: list(P.iterator_range(3, 10))
                [[1, 1], [3], [2, 1], [1, 1, 1], [4], [3, 1], [2, 2]]
                sage: list(P.iterator_range(3, 10, 2))
                [[1, 1], [2, 1], [4], [2, 2]]
                sage: it = P.iterator_range(3)
                sage: [next(it) for x in range(10)]
                [[1, 1],
                 [3], [2, 1], [1, 1, 1],
                 [4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1],
                 [5]]
                sage: it = P.iterator_range(3, step=2)
                sage: [next(it) for x in range(5)]
                [[1, 1],
                 [2, 1],
                 [4], [2, 2], [1, 1, 1, 1]]
                sage: next(P.iterator_range(stop=-3))
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set
                sage: next(P.iterator_range(start=-3))
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set
            """
            if stop is None:
                if start is None:
                    if step is None:
                        yield from self
                        return
                    start = 0
                elif start < 0:
                    yield from self.tuple()[start::step]
                    return
                if step is None:
                    step = 1
                while True:
                    try:
                        yield self.unrank(start)
                    except ValueError:
                        return
                    start += step

            elif stop < 0:
                yield from self.tuple()[start:stop:step]
                return

            if start is None:
                if step is None:
                    it = self.__iter__()
                    for j in range(stop):
                        yield next(it)
                    return
                start = 0
            elif start < 0:
                yield from self.tuple()[start:stop:step]
            if step is None:
                step = 1
            for j in range(start, stop, step):
                yield self.unrank(j)

        def unrank_range(self, start=None, stop=None, step=None):
            """
            Return the range of elements of ``self`` starting at ``start``,
            ending at ``stop``, and stepping by ``step``.

            .. SEEALSO::

                ``unrank()``, ``iterator_range()``

            EXAMPLES::

                sage: # needs sage.combinat
                sage: P = Partitions()
                sage: P.unrank_range(stop=5)
                [[], [1], [2], [1, 1], [3]]
                sage: P.unrank_range(0, 5)
                [[], [1], [2], [1, 1], [3]]
                sage: P.unrank_range(3, 5)
                [[1, 1], [3]]
                sage: P.unrank_range(3, 10)
                [[1, 1], [3], [2, 1], [1, 1, 1], [4], [3, 1], [2, 2]]
                sage: P.unrank_range(3, 10, 2)
                [[1, 1], [2, 1], [4], [2, 2]]
                sage: P.unrank_range(3)
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set
                sage: P.unrank_range(stop=-3)
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set
                sage: P.unrank_range(start=-3)
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set
            """
            if stop is None:
                return list(self.tuple()[start::step])

            if stop < 0:
                return list(self.tuple()[start:stop:step])

            if start is not None and start < 0:
                return list(self.tuple()[start:stop:step])

            return list(self.iterator_range(start, stop, step))

        def __getitem__(self, i):
            r"""
            Return the item indexed by ``i``.

            .. WARNING::

                This method is only meant as a convenience shorthand for
                ``self.unrank(i)`` and
                ``self.unrank_range(start, stop, step)`` respectively, for
                casual use (e.g. in interactive sessions). Subclasses are
                hereby explicitly permitted to overload ``__getitem__``
                with a different semantic, typically for enumerated sets
                that are naturally indexed by some `I` not of the
                form `\{0, 1, \ldots\}`. In particular, generic code
                *should not* use this shorthand.

            EXAMPLES::

                sage: # needs sage.combinat
                sage: P = Partitions()
                sage: P[:5]
                [[], [1], [2], [1, 1], [3]]
                sage: P[0:5]
                [[], [1], [2], [1, 1], [3]]
                sage: P[3:5]
                [[1, 1], [3]]
                sage: P[3:10]
                [[1, 1], [3], [2, 1], [1, 1, 1], [4], [3, 1], [2, 2]]
                sage: P[3:10:2]
                [[1, 1], [2, 1], [4], [2, 2]]
                sage: P[3:]
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set
                sage: P[3]
                [1, 1]
                sage: P[-1]
                Traceback (most recent call last):
                ...
                NotImplementedError: cannot list an infinite set

            ::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.list()
                [1, 2, 3]
                sage: C[1]
                2
                sage: C[:]
                [1, 2, 3]
                sage: C[1:]
                [2, 3]
                sage: C[0:1:2]
                [1]

                sage: F = FiniteEnumeratedSet([1,2,3])
                sage: F[1:]
                [2, 3]
                sage: F[:2]
                [1, 2]
                sage: F[:2:2]
                [1]
                sage: F[1::2]
                [2]
            """
            if isinstance(i, slice):
                return self.unrank_range(i.start, i.stop, i.step)
            if i < 0:
                return self.list()[i]
            return self.unrank(i)

        def __len__(self):
            """
            Return the number of elements of ``self``.

            EXAMPLES::

                sage: len(GF(5))
                5
                sage: len(MatrixSpace(GF(2), 3, 3))                                     # needs sage.modules
                512
            """
            from sage.rings.infinity import Infinity
            try:
                c = self.cardinality()
                if c is Infinity:
                    raise NotImplementedError('infinite set')
                return int(c)
            except AttributeError:
                return len(self.tuple())

        def tuple(self):
            r"""
            Return a tuple of the elements of ``self``.

            The tuple of elements of ``x`` is created and cached on the first call
            of ``x.tuple()``. Each following call of ``x.tuple()`` returns the same tuple.

            For looping, it may be better to do ``for e in x:``, not ``for e in x.tuple():``.

            If ``x`` is not known to be finite, then an exception is raised.

            EXAMPLES::

                sage: (GF(3)^2).tuple()                                                 # needs sage.modules
                ((0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2))
                sage: R = Integers(11)
                sage: l = R.tuple(); l
                (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
                sage: l is R.tuple()
                True
            """
            try: # shortcut
                if self._list is not None:
                    return self._tuple_from_list()
            except AttributeError:
                pass

            if self.list != self._list_default:
                return tuple(self.list())

            from sage.rings.infinity import Infinity
            try:
                if self.cardinality() is Infinity:
                    raise NotImplementedError('cannot list an infinite set')
                else: # finite cardinality
                    return self._tuple_from_iterator()
            except AttributeError:
                raise NotImplementedError('unknown cardinality')
        _tuple_default = tuple

        def _tuple_from_iterator(self):
            r"""
            Return a tuple of the elements of ``self``.

            This implementation of :meth:`tuple` creates the tuple of elements and caches it for
            later uses.

            TESTS::

                sage: R = Integers(11)
                sage: R._list is None
                False
                sage: R._tuple_from_iterator()
                (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
                sage: _ is R._list
                True
            """
            # This creates one throw-away list.
            self._list_from_iterator()
            return self._tuple_from_list()

        def _tuple_from_list(self):
            r"""
            Return a tuple of the elements of ``self``.

            This implementation of :meth:`tuple` assumes that the tuple of elements is already
            cached and just returns it.

            TESTS::

                sage: R = Integers(11)
                sage: R.tuple()
                (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
                sage: R._tuple_from_list() is R.tuple()
                True
            """
            # Implementation classes may put any Sequence type in self._list.
            # Traditionally, self._list was an actual list.
            # When self._list is already a tuple, calling tuple on it is a no-op.
            return tuple(self._list)

        def list(self):
            r"""
            Return a list of the elements of ``self``.

            The elements of set ``x`` are created and cached on the first call
            of ``x.list()``. Then each call of ``x.list()`` returns a new list
            from the cached result. Thus in looping, it may be better to do
            ``for e in x:``, not ``for e in x.list():``.

            If ``x`` is not known to be finite, then an exception is raised.

            EXAMPLES::

                sage: (GF(3)^2).list()                                                  # needs sage.modules
                [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]
                sage: R = Integers(11)
                sage: R.list()
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                sage: l = R.list(); l
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                sage: l.remove(0); l
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                sage: R.list()
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

                sage: C = FiniteEnumeratedSets().example()
                sage: C.list()
                [1, 2, 3]
            """
            return list(self.tuple())
        _list_default = list # needed by the check system.

        def _list_from_iterator(self):
            r"""
            Return a list of the elements of ``self`` after cached.

            TESTS:

            Trying to list an infinite vector space raises an error
            instead of running forever (see :issue:`10470`)::

                sage: (QQ^2).list()  # indirect test                                    # needs sage.modules
                Traceback (most recent call last):
                ...
                AttributeError: 'FreeModule_ambient_field_with_category' object has no attribute 'list'...

            Here we test that for an object that does not know whether it
            is finite or not.  Calling ``x.list()`` simply tries to create
            the list (but here it fails, since the object is not
            iterable). This was fixed :issue:`11350` ::

                sage: R.<t,p> = QQ[]
                sage: Q = R.quotient(t^2-t+1)
                sage: Q.is_finite()
                Traceback (most recent call last):
                ...
                AttributeError: 'QuotientRing_generic_with_category' object has no attribute 'is_finite'...
                sage: Q.list()   # indirect test
                Traceback (most recent call last):
                ...
                AttributeError: 'QuotientRing_generic_with_category' object has no attribute 'list'...

            Here is another example. We artificially create a version of
            the ring of integers that does not know whether it is finite
            or not::

                sage: from sage.rings.integer_ring import IntegerRing_class
                sage: class MyIntegers_class(IntegerRing_class):
                ....:      def is_finite(self):
                ....:          raise NotImplementedError
                sage: MyIntegers = MyIntegers_class()
                sage: MyIntegers.is_finite()
                Traceback (most recent call last):
                ...
                NotImplementedError

            Asking for ``list(MyIntegers)`` will also raise an exception::

                sage: list(MyIntegers)  # indirect test
                Traceback (most recent call last):
                ...
                NotImplementedError: infinite set
            """
            try:
                if self._list is not None:
                    return list(self._list)
            except AttributeError:
                pass
            result = tuple(self.__iter__())
            try:
                self._list = result
            except AttributeError:
                pass
            return list(result)

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
            The ``r``-th element of ``self``.

            ``self.unrank(r)`` returns the ``r``-th element of ``self``, where
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
                ValueError: the rank must be in the range from 0 to 2
                sage: ZZ._unrank_from_iterator(-1)
                Traceback (most recent call last):
                ...
                ValueError: the rank must be greater than or equal to 0
            """
            from sage.rings.integer_ring import ZZ
            if r < 0:
                raise ValueError("the rank must be greater than or equal to 0")
            if r not in ZZ:
                raise ValueError(f"{r=} must be an integer")
            for counter, u in enumerate(self):
                if counter == r:
                    return u
            raise ValueError("the rank must be in the range from %s to %s" % (0,counter))
        unrank = _unrank_from_iterator

        def _rank_from_iterator(self, x):
            """
            The rank of an element of ``self``.

            ``self.rank(x)`` returns the rank of `x`, that is its
            position in the enumeration of ``self``. This is an
            integer between ``0`` and ``n-1`` where ``n`` is the
            cardinality of ``self``, or None if `x` is not in ``self``.

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
            yield from self.tuple()

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
            while not (f is None or f is False):
                yield f
                try:
                    f = self.next(f)
                except (TypeError, ValueError):
                    f = None

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
            Return the first element of ``self`` returned by :meth:`__iter__`.

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
                Cached version of <function ..._an_element_from_iterator at ...>
            """
            it = iter(self)
            try:
                return next(it)
            except StopIteration:
                raise EmptySetError

        # Should this be implemented from first instead?
        _an_element_ = _an_element_from_iterator

        #FIXME: use combinatorial_class_from_iterator once class_from_iterator.patch is in
        def _some_elements_from_iterator(self):
            """
            Return some elements in ``self``.

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
            Return a random element in ``self``.

            Unless otherwise stated, and for finite enumerated sets,
            the probability is uniform.

            This is a generic implementation from the category
            ``EnumeratedSets()``. It raises a :exc:`NotImplementedError`
            since one does not know whether the set is finite.

            EXAMPLES::

                sage: class broken(UniqueRepresentation, Parent):
                ....:  def __init__(self):
                ....:      Parent.__init__(self, category = EnumeratedSets())
                sage: broken().random_element()
                Traceback (most recent call last):
                ...
                NotImplementedError: unknown cardinality
            """
            raise NotImplementedError("unknown cardinality")

        def map(self, f, name=None, *, is_injective=True):
            r"""
            Return the image `\{f(x) | x \in \text{self}\}` of this
            enumerated set by `f`, as an enumerated set.

            INPUT:

            - ``is_injective`` -- boolean (default: ``True``); whether to assume
              that `f` is injective

            EXAMPLES::

                sage: R = Compositions(4).map(attrcall('partial_sums')); R
                Image of Compositions of 4 by The map *.partial_sums()
                 from Compositions of 4
                sage: R.cardinality()
                8
                sage: R.list()
                [[1, 2, 3, 4], [1, 2, 4], [1, 3, 4], [1, 4], [2, 3, 4], [2, 4], [3, 4], [4]]
                sage: [r for r in R]
                [[1, 2, 3, 4], [1, 2, 4], [1, 3, 4], [1, 4], [2, 3, 4], [2, 4], [3, 4], [4]]
                sage: R.category()
                Category of finite enumerated subobjects of sets

            .. WARNING::

                If the function is not injective, then there may be
                repeated elements::

                    sage: P = Compositions(4)
                    sage: P.list()
                    [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]
                    sage: P.map(attrcall('major_index')).list()
                    [6, 3, 4, 1, 5, 2, 3, 0]

                Pass ``is_injective=False`` to get a correct result in this case::

                    sage: P.map(attrcall('major_index'), is_injective=False).list()
                    [6, 3, 4, 1, 5, 2, 0]

            TESTS::

                sage: TestSuite(R).run(skip=['_test_an_element',
                ....:                        '_test_enumerated_set_contains',
                ....:                        '_test_some_elements'])
            """
            from sage.sets.image_set import ImageSubobject

            image = ImageSubobject(f, self, is_injective=is_injective)
            if name:
                image.rename(name)
            return image

#
#  Consistency test suite for an enumerated set:
#
        def _test_enumerated_set_contains(self, **options):
            """
            Check that the methods :meth:`.__contains__` and :meth:`.__iter__` are consistent.

            See also :class:`TestSuite`.

            TESTS::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_contains()
                sage: TestSuite(C).run()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ....:     def __contains__(self, obj):
                ....:         if obj == 3:
                ....:             return False
                ....:         else:
                ....:             return obj in C
                sage: CC = CCls()
                sage: CC._test_enumerated_set_contains()
                Traceback (most recent call last):
                ...
                AssertionError: 3 not found in An example
                of a finite enumerated set: {1,2,3}
            """
            tester = self._tester(**options)
            i = 0
            for w in self:
                tester.assertIn(w, self)
                i += 1
                if i > tester._max_runs:
                    return

        def _test_enumerated_set_iter_list(self, **options):
            """
            Check that the methods :meth:`.list` and :meth:`.__iter__` are consistent.

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
                ....:     def list(self):
                ....:         return [1,2,3,4]
                sage: CC = CCls()
                sage: CC._test_enumerated_set_iter_list()
                Traceback (most recent call last):
                ...
                AssertionError: 3 != 4

            For a large enumerated set this test does nothing:
            increase tester._max_runs if you want to actually run the
            test::

                sage: class CCls(Example):
                ....:     def list(self):
                ....:         return [1,2,3]
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
            Return the rank of ``self`` in its parent.

            See also :meth:`EnumeratedSets.ElementMethods.rank`

            EXAMPLES::

                sage: F = FiniteSemigroups().example(('a','b','c'))
                sage: L = list(F)
                sage: L[7].rank()
                7
                sage: all(x.rank() == i for i,x in enumerate(L))
                True
            """
            return self.parent().rank(self)

    Finite = LazyImport('sage.categories.finite_enumerated_sets', 'FiniteEnumeratedSets', at_startup=True)
    Infinite = LazyImport('sage.categories.infinite_enumerated_sets', 'InfiniteEnumeratedSets', at_startup=True)

    class CartesianProducts(CartesianProductsCategory):

        class ParentMethods:

            def first(self):
                r"""
                Return the first element.

                EXAMPLES::

                    sage: cartesian_product([ZZ]*10).first()
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                """
                return self._cartesian_product_of_elements(
                        tuple(c.first() for c in self.cartesian_factors()))
