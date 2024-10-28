# sage_setup: distribution = sagemath-objects
r"""
Elements, Array and Lists With Clone Protocol, demonstration classes

This module demonstrate the usage of the various classes defined in
:mod:`~sage.structure.list_clone`
"""

#*****************************************************************************
#  Copyright (C) 2011 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.sets_cat import Sets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone cimport (
    ClonableArray, ClonableList, NormalizedClonableList, ClonableIntArray )
from sage.structure.parent import Parent

cdef class IncreasingArray(ClonableArray):
    """
    A small extension class for testing
    :class:`~sage.structure.list_clone.ClonableArray`.

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingArrays
        sage: TestSuite(IncreasingArrays()([1,2,3])).run()
        sage: TestSuite(IncreasingArrays()([])).run()
    """

    cpdef check(self):
        """
        Check that ``self`` is increasing.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()([1,2,3]) # indirect doctest
            [1, 2, 3]
            sage: IncreasingArrays()([3,2,1]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        cdef int i
        for i in range(len(self)-1):
            if self._getitem(i) > self._getitem(i+1):
                raise ValueError("array is not increasing")


class IncreasingArrays(UniqueRepresentation, Parent):
    """
    A small (incomplete) parent for testing
    :class:`~sage.structure.list_clone.ClonableArray`

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingArrays
        sage: IncreasingArrays().element_class
        <... 'sage.structure.list_clone_demo.IncreasingArray'>
    """

    def __init__(self):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()
            <sage.structure.list_clone_demo.IncreasingArrays_with_category object at ...>
            sage: IncreasingArrays() == IncreasingArrays()
            True
        """
        Parent.__init__(self, category = Sets())

    def _element_constructor_(self, *args, **keywords):
        """
        TESTS::

            sage: from sage.structure.list_clone_demo import IncreasingArrays
            sage: IncreasingArrays()([1])     # indirect doctest
            [1]
        """
        return self.element_class(self, *args, **keywords)

    Element = IncreasingArray


class IncreasingLists(IncreasingArrays):
    """
    A small (incomplete) parent for testing
    :class:`~sage.structure.list_clone.ClonableList`

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingLists
        sage: IncreasingLists().element_class
        <... 'sage.structure.list_clone_demo.IncreasingList'>
    """
    Element = IncreasingList


cdef class IncreasingList(ClonableList):
    """
    A small extension class for testing
    :class:`~sage.structure.list_clone.ClonableList`

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingLists
        sage: TestSuite(IncreasingLists()([1,2,3])).run()
        sage: TestSuite(IncreasingLists()([])).run()
    """

    cpdef check(self):
        """
        Check that ``self`` is increasing.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingLists
            sage: IncreasingLists()([1,2,3]) # indirect doctest
            [1, 2, 3]
            sage: IncreasingLists()([3,2,1]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        cdef int i
        for i in range(len(self)-1):
            if self._getitem(i) >= self._getitem(i+1):
                raise ValueError("array is not increasing")


cdef class IncreasingIntArray(ClonableIntArray):
    """
    A small extension class for testing
    :class:`~sage.structure.list_clone.ClonableIntArray`.

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: TestSuite(IncreasingIntArrays()([1,2,3])).run()
        sage: TestSuite(IncreasingIntArrays()([])).run()
    """

    cpdef check(self):
        """
        Check that ``self`` is increasing.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import IncreasingIntArrays
            sage: IncreasingIntArrays()([1,2,3]) # indirect doctest
            [1, 2, 3]
            sage: IncreasingIntArrays()([3,2,1]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: array is not increasing
        """
        cdef int i
        if not self:
            return
        for i in range(len(self)-1):
            if self._getitem(i) >= self._getitem(i+1):
                raise ValueError("array is not increasing")


class IncreasingIntArrays(IncreasingArrays):
    """
    A small (incomplete) parent for testing
    :class:`~sage.structure.list_clone.ClonableIntArray`

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: IncreasingIntArrays().element_class
        <... 'sage.structure.list_clone_demo.IncreasingIntArray'>
    """
    Element = IncreasingIntArray


cdef class SortedList(NormalizedClonableList):
    """
    A small extension class for testing
    :class:`~sage.structure.list_clone.NormalizedClonableList`.

    TESTS::

        sage: from sage.structure.list_clone_demo import IncreasingIntArrays
        sage: TestSuite(IncreasingIntArrays()([1,2,3])).run()
        sage: TestSuite(IncreasingIntArrays()([])).run()
    """
    cpdef normalize(self):
        """
        Normalize ``self``.

        Sort the list stored in ``self``.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import SortedList, SortedLists
            sage: l = SortedList(SortedLists(), [3,1,2], False, False)
            sage: l         # indirect doctest
            [1, 2, 3]
            sage: l[1] = 5; l
            [1, 5, 3]
            sage: l.normalize(); l
            [1, 3, 5]
        """
        self._require_mutable()
        self._get_list().sort()

    cpdef check(self):
        """
        Check that ``self`` is strictly increasing.

        EXAMPLES::

            sage: from sage.structure.list_clone_demo import SortedList, SortedLists
            sage: SortedLists()([1,2,3]) # indirect doctest
            [1, 2, 3]
            sage: SortedLists()([3,2,2]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: list is not strictly increasing
        """
        for i in range(len(self)-1):
            if self._getitem(i) >= self._getitem(i+1):
                raise ValueError("list is not strictly increasing")


class SortedLists(IncreasingLists):
    """
    A small (incomplete) parent for testing
    :class:`~sage.structure.list_clone.NormalizedClonableList`

    TESTS::

        sage: from sage.structure.list_clone_demo import SortedList, SortedLists
        sage: SL = SortedLists()
        sage: SL([3,1,2])
        [1, 2, 3]
    """
    Element = SortedList
