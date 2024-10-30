r"""
Enumerated set of lists of integers with constraints: front-end

- :class:`IntegerLists`: class which models an enumerated set of lists
  of integers with certain constraints. This is a Python front-end
  where all user-accessible functionality should be implemented.
"""

#*****************************************************************************
#       Copyright (C) 2015 Bryan Gillespie <Brg008@gmail.com>
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#                          Anne Schilling <anne@math.ucdavis.edu>
#                          Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from inspect import ismethod
from sage.categories.enumerated_sets import EnumeratedSets
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.combinat.integer_lists.base import IntegerListsBackend


class IntegerList(ClonableArray):
    """
    Element class for :class:`IntegerLists`.
    """
    def check(self):
        """
        Check to make sure this is a valid element in its
        :class:`IntegerLists` parent.

        EXAMPLES::

            sage: C = IntegerListsLex(4)
            sage: C([4]).check()
            True
            sage: C([5]).check()
            False
        """
        return self in self.parent()


class IntegerLists(Parent):
    """
    Enumerated set of lists of integers with constraints.

    Currently, this is simply an abstract base class which should not
    be used by itself. See :class:`IntegerListsLex` for a class which
    can be used by end users.

    ``IntegerLists`` is just a Python front-end, acting as a
    :class:`Parent`, supporting element classes and so on.
    The attribute ``.backend`` which is an instance of
    :class:`sage.combinat.integer_lists.base.IntegerListsBackend` is the
    Cython back-end which implements all operations such as iteration.

    The front-end (i.e. this class) and the back-end are supposed to be
    orthogonal: there is no imposed correspondence between front-ends
    and back-ends.

    For example, the set of partitions of 5 and the set of weakly
    decreasing sequences which sum to 5 might be implemented by the
    same back-end, but they will be presented to the user by a
    different front-end.

    EXAMPLES::

        sage: from sage.combinat.integer_lists import IntegerLists
        sage: L = IntegerLists(5)
        sage: L
        Integer lists of sum 5 satisfying certain constraints

        sage: IntegerListsLex(2, length=3, name="A given name")
        A given name
    """
    backend = None
    backend_class = IntegerListsBackend

    Element = IntegerList

    def __init__(self, *args, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.integer_lists import IntegerLists
            sage: C = IntegerLists(2, length=3)
            sage: C == loads(dumps(C))
            True
        """
        if "name" in kwds:
            self.rename(kwds.pop("name"))

        if "element_class" in kwds:
            self.Element = kwds.pop("element_class")

        if "element_constructor" in kwds:
            element_constructor = kwds.pop("element_constructor")
        elif issubclass(self.Element, ClonableArray):
            # Not all element classes support check=False, but
            # ClonableArray certainly does.
            element_constructor = self._element_constructor_nocheck
        else:
            element_constructor = self._element_constructor_default
        self._element_constructor_ = element_constructor

        category = kwds.pop("category", None)
        if category is None:
            category = EnumeratedSets().Finite()

        # Let self.backend be some IntegerListsBackend
        self.backend = self.backend_class(*args, **kwds)

        Parent.__init__(self, category=category)

    def __eq__(self, other):
        r"""
        Return whether ``self == other``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: D = IntegerListsLex(2, length=3); L = D.list()
            sage: E = IntegerListsLex(2, min_length=3)
            sage: F = IntegerListsLex(2, length=3, element_constructor=list)
            sage: G = IntegerListsLex(4, length=3)
            sage: C == C
            True
            sage: C == D
            True
            sage: C == E
            False
            sage: C == F
            False
            sage: C == None
            False
            sage: C == G
            False

        This is a minimal implementation enabling pickling tests. It
        is safe, but one would want the two following objects to be
        detected as equal::

            sage: C = IntegerListsLex(2, ceiling=[1,1,1])
            sage: D = IntegerListsLex(2, ceiling=[1,1,1])
            sage: C == D
            False

        TESTS:

        This used to fail due to poor equality testing. See
        :issue:`17979`, comment 433::

            sage: DisjointUnionEnumeratedSets(Family([2,2],
            ....:     lambda n: IntegerListsLex(n, length=2))).list()
            [[2, 0], [1, 1], [0, 2], [2, 0], [1, 1], [0, 2]]
            sage: DisjointUnionEnumeratedSets(Family([2,2],
            ....:     lambda n: IntegerListsLex(n, length=1))).list()
            [[2], [2]]
        """
        if self.__class__ != other.__class__:
            return False
        if self.backend != other.backend:
            return False
        a = self._element_constructor_
        b = other._element_constructor_
        if ismethod(a):
            a = a.__func__
        if ismethod(b):
            b = b.__func__
        return a == b

    def __ne__(self, other):
        r"""
        Return whether ``self != other``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: D = IntegerListsLex(2, length=3); L = D.list()
            sage: E = IntegerListsLex(2, max_length=3)
            sage: C != D
            False
            sage: C != E
            True
        """
        return not self == other

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: D = IntegerListsLex(2, max_length=3)
            sage: hash(C) == hash(C)
            True
        """
        a = self._element_constructor_
        if ismethod(a):
            a = a.__func__
        return hash((self.__class__, a))

    def __iter__(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C)     # indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        return self._element_iter(self.backend._iter(), self._element_constructor_)

    @staticmethod
    def _element_iter(itr, constructor):
        """
        Given an iterator ``itr`` and an element constructor
        ``constructor``, iterate over ``constructor(v)`` where `v`
        are the values of ``itr``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C._element_iter(C._iter(), tuple))
            [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
        """
        for v in itr:
            yield constructor(v)

    def __getattr__(self, name):
        """
        Get an attribute of the implementation backend.

        Ideally, this would be done using multiple inheritance, but
        Python doesn't support that for built-in types.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: C.min_length
            3

        TESTS:

        Check that uninitialized instances do not lead to infinite
        recursion because there is no ``backend`` attribute::

            sage: from sage.combinat.integer_lists import IntegerLists
            sage: L = IntegerLists.__new__(IntegerLists)
            sage: L.foo
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute 'foo'...
        """
        return getattr(self.backend, name)

    def __contains__(self, item):
        """
        Return ``True`` if ``item`` meets the constraints imposed by
        the arguments.

        EXAMPLES::

            sage: C = IntegerListsLex(n=2, max_length=3, min_slope=0)
            sage: all(l in C for l in C)
            True
        """
        return self.backend._contains(item)

    def _element_constructor_default(self, l):
        """
        Default element constructor.

        EXAMPLES::

            sage: L = IntegerListsLex(4)
            sage: L._element_constructor_default([1,2,3])
            [1, 2, 3]

        We construct a variant of :class:`IntegerLists` with a custom
        element class::

            sage: class MyElt(list):
            ....:     def __init__(self, parent, x):
            ....:         list.__init__(self, x)
            sage: from sage.combinat.integer_lists import IntegerLists
            sage: class MyIntegersLists(IntegerLists):
            ....:     Element = MyElt
            sage: L = MyIntegersLists(5)
            sage: L._element_constructor_
            <bound method IntegerLists._element_constructor_default of Integer lists of sum 5 satisfying certain constraints>
        """
        return self.element_class(self, l)

    def _element_constructor_nocheck(self, l):
        r"""
        A variant of the standard element constructor that passes
        ``check=False`` to the element class.

        EXAMPLES::

            sage: L = IntegerListsLex(4, max_slope=0)
            sage: L._element_constructor_nocheck([1,2,3])
            [1, 2, 3]

        When relevant, this is assigned to
        ``self._element_constructor_`` by :meth:`__init__`, to avoid
        overhead when constructing elements from trusted data in the
        iterator::

            sage: L._element_constructor_
            <bound method IntegerLists._element_constructor_nocheck of ...>
            sage: L._element_constructor_([1,2,3])
            [1, 2, 3]
        """
        return self.element_class(self, l, check=False)
