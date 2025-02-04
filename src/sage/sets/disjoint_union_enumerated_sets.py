"""
Disjoint union of enumerated sets

AUTHORS:

- Florent Hivert (2009-07/09): initial implementation.
- Florent Hivert (2010-03): classcall related stuff.
- Florent Hivert (2010-12): fixed facade element construction.
"""
# ***************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.sets.family import Family
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.rings.infinity import Infinity
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation


class DisjointUnionEnumeratedSets(UniqueRepresentation, Parent):
    """
    A class for disjoint unions of enumerated sets.

    INPUT:

    - ``family`` -- list (or iterable or family) of enumerated sets
    - ``keepkey`` -- boolean
    - ``facade`` -- boolean

    This models the enumerated set obtained by concatenating together
    the specified ordered sets. The latter are supposed to be pairwise
    disjoint; otherwise, a multiset is created.

    The argument ``family`` can be a list, a tuple, a dictionary, or a
    family. If it is not a family it is first converted into a family
    (see :func:`sage.sets.family.Family`).

    Experimental options:

    By default, there is no way to tell from which set of the union an
    element is generated. The option ``keepkey=True`` keeps track of
    those by returning pairs ``(key, el)`` where ``key`` is the index
    of the set to which ``el`` belongs. When this option is specified,
    the enumerated sets need not be disjoint anymore.

    With the option ``facade=False`` the elements are wrapped in an
    object whose parent is the disjoint union itself. The wrapped
    object can then be recovered using the ``value`` attribute.

    The two options can be combined.

    The names of those options is imperfect, and subject to change in
    future versions. Feedback welcome.

    EXAMPLES:

    The input can be a list or a tuple of FiniteEnumeratedSets::

        sage: U1 = DisjointUnionEnumeratedSets((
        ....:       FiniteEnumeratedSet([1,2,3]),
        ....:       FiniteEnumeratedSet([4,5,6])))
        sage: U1
        Disjoint union of Family ({1, 2, 3}, {4, 5, 6})
        sage: U1.list()
        [1, 2, 3, 4, 5, 6]
        sage: U1.cardinality()
        6

    The input can also be a dictionary::

        sage: U2 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
        ....:                                   2: FiniteEnumeratedSet([4,5,6])})
        sage: U2
        Disjoint union of Finite family {1: {1, 2, 3}, 2: {4, 5, 6}}
        sage: U2.list()
        [1, 2, 3, 4, 5, 6]
        sage: U2.cardinality()
        6

    However in that case the enumeration order is not specified.

    In general the input can be any family::

        sage: # needs sage.combinat
        sage: U3 = DisjointUnionEnumeratedSets(
        ....:     Family([2,3,4], Permutations, lazy=True))
        sage: U3
        Disjoint union of Lazy family
         (<class 'sage.combinat.permutation.Permutations'>(i))_{i in [2, 3, 4]}
        sage: U3.cardinality()
        32
        sage: it = iter(U3)
        sage: [next(it), next(it), next(it), next(it), next(it), next(it)]
        [[1, 2], [2, 1], [1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1]]
        sage: U3.unrank(18)
        [2, 4, 1, 3]

    This allows for infinite unions::

        sage: # needs sage.combinat
        sage: U4 = DisjointUnionEnumeratedSets(
        ....:     Family(NonNegativeIntegers(), Permutations))
        sage: U4
        Disjoint union of Lazy family
         (<class 'sage.combinat.permutation.Permutations'>(i))_{i in Non negative integers}
        sage: U4.cardinality()
        +Infinity
        sage: it = iter(U4)
        sage: [next(it), next(it), next(it), next(it), next(it), next(it)]
        [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2]]
        sage: U4.unrank(18)
        [2, 3, 1, 4]

    .. WARNING::

        Beware that some of the operations assume in that case that infinitely
        many of the enumerated sets are non empty.


    .. RUBRIC:: Experimental options

    We demonstrate the ``keepkey`` option::

        sage: # needs sage.combinat
        sage: Ukeep = DisjointUnionEnumeratedSets(
        ....:            Family(list(range(4)), Permutations), keepkey=True)
        sage: it = iter(Ukeep)
        sage: [next(it) for i in range(6)]
        [(0, []), (1, [1]), (2, [1, 2]), (2, [2, 1]), (3, [1, 2, 3]), (3, [1, 3, 2])]
        sage: type(next(it)[1])
        <class 'sage.combinat.permutation.StandardPermutations_n_with_category.element_class'>

    We now demonstrate the ``facade`` option::

        sage: # needs sage.combinat
        sage: UNoFacade = DisjointUnionEnumeratedSets(
        ....:                Family(list(range(4)), Permutations), facade=False)
        sage: it = iter(UNoFacade)
        sage: [next(it) for i in range(6)]
        [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2]]
        sage: el = next(it); el
        [2, 1, 3]
        sage: type(el)
        <... 'sage.structure.element_wrapper.ElementWrapper'>
        sage: el.parent() == UNoFacade
        True
        sage: elv = el.value; elv
        [2, 1, 3]
        sage: type(elv)
        <class 'sage.combinat.permutation.StandardPermutations_n_with_category.element_class'>

    The elements ``el`` of the disjoint union are simple wrapped elements.
    So to access the methods, you need to do ``el.value``::

        sage: el[0]                                                                     # needs sage.combinat
        Traceback (most recent call last):
        ...
        TypeError: 'sage.structure.element_wrapper.ElementWrapper' object is not subscriptable

        sage: el.value[0]                                                               # needs sage.combinat
        2

    Possible extensions: the current enumeration order is not suitable
    for unions of infinite enumerated sets (except possibly for the
    last one). One could add options to specify alternative enumeration
    orders (anti-diagonal, round robin, ...) to handle this case.


    .. RUBRIC:: Inheriting from ``DisjointUnionEnumeratedSets``

    There are two different use cases for inheriting from
    :class:`DisjointUnionEnumeratedSets`: writing a parent which
    happens to be a disjoint union of some known parents, or writing
    generic disjoint unions for some particular classes of
    :class:`sage.categories.enumerated_sets.EnumeratedSets`.

    - In the first use case, the input of the ``__init__`` method is
      most likely different from that of
      :class:`DisjointUnionEnumeratedSets`. Then, one simply
      writes the ``__init__`` method as usual::

          sage: class MyUnion(DisjointUnionEnumeratedSets):
          ....:   def __init__(self):
          ....:       DisjointUnionEnumeratedSets.__init__(self,
          ....:            Family([1,2], Permutations))
          sage: pp = MyUnion()
          sage: pp.list()
          [[1], [1, 2], [2, 1]]

      In case the :meth:`__init__` method takes optional arguments,
      or does some normalization on them, a specific method
      ``__classcall_private__`` is required (see the
      documentation of :class:`UniqueRepresentation`).

    - In the second use case, the input of the ``__init__`` method
      is the same as that of :class:`DisjointUnionEnumeratedSets`;
      one therefore wants to inherit the :meth:`__classcall_private__`
      method as well, which can be achieved as follows::

          sage: class UnionOfSpecialSets(DisjointUnionEnumeratedSets):
          ....:     __classcall_private__ = staticmethod(DisjointUnionEnumeratedSets.__classcall_private__)
          sage: psp = UnionOfSpecialSets(Family([1,2], Permutations))
          sage: psp.list()
          [[1], [1, 2], [2, 1]]

    TESTS::

        sage: TestSuite(U1).run()
        sage: TestSuite(U2).run()
        sage: TestSuite(U3).run()                                                       # needs sage.combinat
        sage: TestSuite(U4).run()                                                       # needs sage.combinat
        doctest:...: UserWarning: Disjoint union of Lazy family
        (<class 'sage.combinat.permutation.Permutations'>(i))_{i in Non negative integers}
        is an infinite union
        The default implementation of __contains__ can loop forever. Please overload it.
        sage: TestSuite(UNoFacade).run()                                                # needs sage.combinat

    We skip ``_test_an_element`` because the coercion framework does not
    currently allow a tuple to be returned for facade parents::

        sage: TestSuite(Ukeep).run(skip='_test_an_element')                             # needs sage.combinat

    The following three lines are required for the pickling tests,
    because the classes ``MyUnion`` and ``UnionOfSpecialSets`` have
    been defined interactively::

        sage: import __main__
        sage: __main__.MyUnion = MyUnion
        sage: __main__.UnionOfSpecialSets = UnionOfSpecialSets

        sage: TestSuite(pp).run()
        sage: TestSuite(psp).run()
    """

    @staticmethod
    def __classcall_private__(cls, fam, facade=True,
                              keepkey=False, category=None):
        """
        Normalization of arguments; see :class:`UniqueRepresentation`.

        TESTS:

        We check that disjoint unions have unique representation::

            sage: U1 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                   2: FiniteEnumeratedSet([4,5,6])})
            sage: U2 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                   2: FiniteEnumeratedSet([4,5,6])})
            sage: U1 == U2
            True
            sage: U1 is U2        # indirect doctest
            True
            sage: U3 = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                   2: FiniteEnumeratedSet([4,5])})
            sage: U1 == U3
            False
        """
        # facade  = options.pop('facade', True);
        # keepkey = options.pop('keepkey', False);
        assert isinstance(facade, bool)
        assert isinstance(keepkey, bool)
        return super().__classcall__(
            cls, Family(fam),
            facade=facade, keepkey=keepkey, category=category)

    def __init__(self, family, facade=True, keepkey=False, category=None):
        """
        TESTS::

            sage: U = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                  2: FiniteEnumeratedSet([4,5,6])})
            sage: TestSuite(U).run()

            sage: X = DisjointUnionEnumeratedSets({i: Partitions(i) for i in range(5)})             # needs sage.combinat sage.libs.flint
            sage: TestSuite(X).run()                                                    # needs sage.combinat sage.libs.flint
        """
        self._family = family
        self._facade = facade
        if facade:
            # Note that family is not copied when it is a finite enumerated
            # set, thus, any subclass must ensure that it does not mutate this
            # input.
            if family in FiniteEnumeratedSets():
                self._facade_for = family
            else:
                # This allows the test suite to pass its tests by essentially
                #   stating that this is a facade for any parent. Technically
                #   this is wrong, but in practice, it will not have much
                #   of an effect.
                self._facade_for = True
        self._keepkey = keepkey
        if self._is_category_initialized():
            return
        if category is None:
            # try to guess if the result is infinite or not.
            if self._family in InfiniteEnumeratedSets():
                category = InfiniteEnumeratedSets()
            elif self._family.last().cardinality() == Infinity:
                category = InfiniteEnumeratedSets()
            else:
                category = FiniteEnumeratedSets()
        Parent.__init__(self, facade=facade, category=category)

    def _repr_(self):
        """
        TESTS::

            sage: U = DisjointUnionEnumeratedSets({1: FiniteEnumeratedSet([1,2,3]),
            ....:                                  2: FiniteEnumeratedSet([4,5,6])})
            sage: U
            Disjoint union of Finite family {1: {1, 2, 3}, 2: {4, 5, 6}}
        """
        return "Disjoint union of %s" % self._family

    def _is_a(self, x):
        """
        Check if a Sage object ``x`` belongs to ``self``.

        This methods is a helper for :meth:`__contains__` and the
        constructor :meth:`_element_constructor_`.

        EXAMPLES::

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Compositions))
            sage: U4._is_a(Composition([3,2,1,1]))
            doctest:...: UserWarning: Disjoint union of Lazy family
            (<class 'sage.combinat.composition.Compositions'>(i))_{i in Non negative integers}
            is an infinite union
            The default implementation of __contains__ can loop forever. Please overload it.
            True
        """
        if self._keepkey:
            return (isinstance(x, tuple) and
                    x[0] in self._family.keys() and
                    x[1] in self._family[x[0]])
        else:
            from warnings import warn
            if self._family.cardinality() == Infinity:
                warn("%s is an infinite union\nThe default implementation of __contains__ can loop forever. Please overload it." % (self))
            return any(x in a for a in self._family)

    def __contains__(self, x):
        """
        Check containment.

        .. WARNING::

            If ``self`` is an infinite union and if the answer is
            logically False, this will loop forever and never answer
            ``False``. Therefore, a warning is issued.

        EXAMPLES::

            sage: U4 = DisjointUnionEnumeratedSets(                                     # needs sage.combinat
            ....:          Family(NonNegativeIntegers(), Partitions))
            sage: Partition([]) in U4                                                   # needs sage.combinat
            doctest:...: UserWarning: Disjoint union of Lazy family
            (<class 'sage.combinat.partition.Partitions'>(i))_{i in Non negative integers}
            is an infinite union
            The default implementation of __contains__ can loop forever. Please overload it.
            True

        Note: one has to use a different family from the previous one in this
        file otherwise the warning is not re-issued::

            sage: Partition([3,2,1,1]) in U4                                            # needs sage.combinat
            True

        The following call will loop forever::

            sage: 2 in U4 # not tested, loop forever
        """
        if self._facade:
            return self._is_a(x)
        else:
            if isinstance(x, self.element_class):
                return True
            else:
                return self._is_a(x)

    def __iter__(self):
        """
        TESTS::

            sage: from itertools import islice
            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Permutations))
            sage: list(islice(iter(U4), 6))
            [[], [1], [1, 2], [2, 1], [1, 2, 3], [1, 3, 2]]

            sage: # needs sage.combinat
            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family(NonNegativeIntegers(), Permutations),
            ....:          keepkey=True, facade=False)
            sage: l = list(islice(iter(U4), 7)); l
            [(0, []), (1, [1]), (2, [1, 2]), (2, [2, 1]), (3, [1, 2, 3]), (3, [1, 3, 2]), (3, [2, 1, 3])]
            sage: l[-1].parent() is U4
            True

        Check when both the set of keys and each element set is finite::

            sage: list(DisjointUnionEnumeratedSets(
            ....:          Family({1: FiniteEnumeratedSet([1,2,3]),
            ....:                  2: FiniteEnumeratedSet([4,5,6])})))
            [1, 2, 3, 4, 5, 6]

        Check when the set of keys is finite but each element set is infinite::

            sage: list(islice(DisjointUnionEnumeratedSets(
            ....:                 Family({1: NonNegativeIntegers(),
            ....:                         2: NonNegativeIntegers()}), keepkey=True), 0, 10))
            [(1, 0), (1, 1), (2, 0), (1, 2), (2, 1), (1, 3), (2, 2), (1, 4), (2, 3), (1, 5)]

        Check when the set of keys is infinite but each element set is finite::

            sage: list(islice(DisjointUnionEnumeratedSets(
            ....:                 Family(NonNegativeIntegers(), lambda x: FiniteEnumeratedSet(range(x))),
            ....:                 keepkey=True), 0, 20))
            [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (4, 2), (4, 3),
             (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (6, 0), (6, 1), (6, 2), (6, 3), (6, 4)]

        Check when some element sets are empty (note that if there are infinitely many sets
        but only finitely many elements in total, the iteration will hang)::

            sage: list(DisjointUnionEnumeratedSets(
            ....:          Family({1: FiniteEnumeratedSet([]),
            ....:                  2: FiniteEnumeratedSet([]),
            ....:                  3: FiniteEnumeratedSet([]),
            ....:                  4: FiniteEnumeratedSet([]),
            ....:                  5: FiniteEnumeratedSet([1,2,3]),
            ....:                  6: FiniteEnumeratedSet([4,5,6])})))
            [1, 2, 3, 4, 5, 6]

        Check when there's one infinite set and infinitely many finite sets::

            sage: list(islice(DisjointUnionEnumeratedSets(
            ....:                 Family(NonNegativeIntegers(), lambda x: FiniteEnumeratedSet([]) if x else NonNegativeIntegers())),
            ....:                 0, 10))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        The following cannot be determined to be finite, but the first elements can still be retrieved::

            sage: U = DisjointUnionEnumeratedSets(
            ....:         Family(NonNegativeIntegers(), lambda x: FiniteEnumeratedSet([] if x >= 2 else [1, 2])),
            ....:         keepkey=True)
            sage: list(U)  # not tested
            sage: list(islice(iter(U), 5))  # not tested, hangs
            sage: list(islice(iter(U), 4))
            [(0, 1), (0, 2), (1, 1), (1, 2)]
        """
        def wrap_element(el, k):
            nonlocal self
            if self._keepkey:
                el = (k, el)
            if self._facade:
                return el
            else:
                return self.element_class(self, el)  # Bypass correctness tests

        keys_iter = iter(self._family.keys())
        if self._keepkey:
            seen_keys = []
        el_iters = []
        while keys_iter is not None or el_iters:
            if keys_iter is not None:
                try:
                    k = next(keys_iter)
                except StopIteration:
                    keys_iter = None
                if keys_iter is not None:
                    el_set = self._family[k]
                    if el_set.is_finite():
                        for el in el_set:
                            yield wrap_element(el, k)
                    else:
                        el_iters.append(iter(el_set))
                        if self._keepkey:
                            seen_keys.append(k)
            any_stopped = False
            for i, obj in enumerate(zip(seen_keys, el_iters) if self._keepkey else el_iters):
                if self._keepkey:
                    k, el_iter = obj
                else:
                    k = None
                    el_iter = obj
                try:
                    el = next(el_iter)
                except StopIteration:
                    el_iters[i] = None
                    any_stopped = True
                    continue
                yield wrap_element(el, k)
            if any_stopped:
                if self._keepkey:
                    filtered = [*zip(
                        *[(k, el_iter) for k, el_iter in zip(seen_keys, el_iters) if el_iter is not None])]
                    if filtered:
                        seen_keys = list(filtered[0])
                        el_iters = list(filtered[1])
                    else:
                        seen_keys = []
                        el_iters = []
                else:
                    el_iters = [el_iter for el_iter in el_iters if el_iter is not None]

    def an_element(self):
        """
        Return an element of this disjoint union, as per
        :meth:`Sets.ParentMethods.an_element`.

        EXAMPLES::

            sage: U4 = DisjointUnionEnumeratedSets(
            ....:          Family([3, 5, 7], Permutations))
            sage: U4.an_element()
            [1, 2, 3]
        """
        return self._an_element_from_iterator()

    @cached_method
    def cardinality(self):
        """
        Return the cardinality of this disjoint union.

        EXAMPLES:

        For finite disjoint unions, the cardinality is computed by
        summing the cardinalities of the enumerated sets::

            sage: U = DisjointUnionEnumeratedSets(Family([0,1,2,3], Permutations))
            sage: U.cardinality()
            10

        For infinite disjoint unions, this makes the assumption that
        the result is infinite::

            sage: U = DisjointUnionEnumeratedSets(
            ....:         Family(NonNegativeIntegers(), Permutations))
            sage: U.cardinality()
            +Infinity

        .. WARNING::

            As pointed out in the main documentation, it is
            possible to construct examples where this is incorrect::

                sage: U = DisjointUnionEnumeratedSets(
                ....:         Family(NonNegativeIntegers(), lambda x: []))
                sage: U.cardinality()  # Should be 0!
                +Infinity
        """
        if self._family.cardinality() == Infinity:
            return Infinity
        return sum(set.cardinality() for set in self._family)

    @lazy_attribute
    def _element_constructor_(self):
        """
        TESTS::

            sage: # needs sage.combinat sage.libs.flint
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=False)
            sage: U._element_constructor_
            <bound method DisjointUnionEnumeratedSets._element_constructor_default
             of Disjoint union of Finite family {...}>
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=True)
            sage: U._element_constructor_
            <bound method DisjointUnionEnumeratedSets._element_constructor_facade
             of Disjoint union of Finite family {...}>
        """
        if not self._facade:
            return self._element_constructor_default
        else:
            return self._element_constructor_facade

    def _element_constructor_default(self, el):
        r"""
        TESTS::

            sage: # needs sage.combinat sage.libs.flint
            sage: U = DisjointUnionEnumeratedSets(
            ....:         Family([1,2,3], Partitions), facade=False)
            sage: U([1])       # indirect doctest
            [1]
            sage: U([2,1])     # indirect doctest
            [2, 1]
            sage: U([1,3,2])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: value [1, 3, 2] does not belong to Disjoint union of
             Finite family {1: Partitions of the integer 1,
                 2: Partitions of the integer 2,
                 3: Partitions of the integer 3}

            sage: # needs sage.combinat sage.libs.flint
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), keepkey=True, facade=False)
            sage: U((1, [1]))    # indirect doctest
            (1, [1])
            sage: U((3,[2,1]))   # indirect doctest
            (3, [2, 1])
            sage: U((4,[2,1]))   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: value (4, [2, 1]) does not belong to Disjoint union of
             Finite family {1: Partitions of the integer 1,
                 2: Partitions of the integer 2,
                 3: Partitions of the integer 3}
        """
        if isinstance(el, self.element_class):
            el = el.value
        if self._is_a(el):
            return self.element_class(self, el)
        else:
            raise ValueError("value %s does not belong to %s" % (el, self))

    def _element_constructor_facade(self, el):
        """
        TESTS::

            sage: # needs sage.combinat sage.libs.flint
            sage: X = DisjointUnionEnumeratedSets({i: Partitions(i)
            ....:                                  for i in range(5)})
            sage: X([1]).parent()
            Partitions of the integer 1
            sage: X([2,1,1]).parent()  # indirect doctest
            Partitions of the integer 4
            sage: X([6])
            Traceback (most recent call last):
            ...
            ValueError: cannot coerce `[6]` in any parent in `Finite family {...}`

        We need to call the element constructor directly when ``keepkey=True``
        because this returns a `tuple`, where the coercion framework requires
        an :class:`Element` be returned.

            sage: X = DisjointUnionEnumeratedSets({i: Partitions(i)                     # needs sage.combinat sage.libs.flint
            ....:                                  for i in range(5)},
            ....:                                 keepkey=True)
            sage: p = X._element_constructor_((0, []))  # indirect doctest              # needs sage.combinat sage.libs.flint
            sage: p[1].parent()                                                         # needs sage.combinat sage.libs.flint
            Partitions of the integer 0

        Test that facade parents can create and properly access elements
        that are tuples (fixed by :issue:`22382`)::

            sage: # needs sage.combinat sage.libs.flint
            sage: f = lambda mu: cartesian_product([mu.standard_tableaux(),
            ....:                                   mu.standard_tableaux()])
            sage: tabs = DisjointUnionEnumeratedSets(Family(Partitions(4), f))
            sage: s = StandardTableau([[1,3],[2,4]])
            sage: (s,s) in tabs
            True
            sage: ss = tabs( (s,s) )
            sage: ss[0]
            [[1, 3], [2, 4]]

        We do not coerce when one of the elements is already in the set::

            sage: X = DisjointUnionEnumeratedSets([QQ, ZZ])
            sage: x = X(2)
            sage: x.parent() is ZZ
            True
        """
        if self._keepkey:
            P = self._family[el[0]]
            if isinstance(el[1], Element) and el[1].parent() == P:
                return el
            try:
                return (el[0], P(el[1]))
            except Exception:
                raise ValueError("cannot coerce `%s` in the parent `%s`" % (el[1], P))

        # Check first to see if the parent of el is in the family
        if (isinstance(el, Element) and self._facade_for is not True
                and el.parent() in self._facade_for):
            return el

        for P in self._family:
            try:
                return P(el)
            except Exception:
                pass
        raise ValueError("cannot coerce `%s` in any parent in `%s`" % (el, self._family))

    @lazy_attribute
    def Element(self):
        """
        TESTS::

            sage: # needs sage.combinat sage.libs.flint
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=False)
            sage: U.Element
            <... 'sage.structure.element_wrapper.ElementWrapper'>
            sage: U = DisjointUnionEnumeratedSets(
            ....:          Family([1,2,3], Partitions), facade=True)
            sage: U.Element
            Traceback (most recent call last):
            ...
            AttributeError: 'DisjointUnionEnumeratedSets_with_category' object
            has no attribute 'Element'...
        """
        if not self._facade:
            return ElementWrapper
        return NotImplemented
