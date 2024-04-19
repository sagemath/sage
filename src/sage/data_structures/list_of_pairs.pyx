# sage_setup: distribution = sagemath-categories
r"""
A data structure to store lists of integer pairs of large size.
"""

# ****************************************************************************
#       Copyright (C) 2022 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_reallocarray, check_allocarray, sig_free
from sage.rings.integer cimport smallInteger

# Should be a power of two.
# Should be neither exposed nor modified.
cdef size_t length_per_list = 16348

cdef class ListOfPairs:
    def __dealloc__(self):
        cdef size_t n_lists = self.length // length_per_list
        cdef size_t i
        for i in range(n_lists):
            sig_free(self._lists[i])
        sig_free(self._lists)

    cdef inline int enlarge(self) except -1:
        """
        Increase size of list by one.
        """
        if self.length % length_per_list:
            self.length += 1
            return 0

        cdef size_t n_lists = self.length // length_per_list
        self._lists = <pair_s**> check_reallocarray(self._lists, n_lists + 1, sizeof(pair_s*))
        self._lists[n_lists] = <pair_s*> check_allocarray(length_per_list, sizeof(pair_s))
        self.length += 1

    cdef inline pair_s* get(self, size_t index) except NULL:
        """
        Return a pointer to a pair of the list corresponding to the ``index``.
        """
        if not (0 <= index < self.length):
            raise IndexError

        cdef size_t list_index = index // length_per_list
        cdef size_t index_in_list = index - list_index * length_per_list

        return &self._lists[list_index][index_in_list]

    def __getitem__(self, size_t index):
        r"""
        Get item of specified index.

        EXAMPLES::

            sage: from sage.data_structures.list_of_pairs import ListOfPairs
            sage: l = ListOfPairs()
            sage: l[0] = [1, 5]
            sage: l[0]
            (1, 5)
            sage: l[1]
            Traceback (most recent call last):
            ...
            IndexError
            sage: l[-1]
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to size_t
        """
        cdef pair_s* pair = self.get(index)
        return (smallInteger(pair.first), smallInteger(pair.second))

    def __setitem__(self, size_t index, value):
        r"""
        Set item of specified index.

        Allows increasing the size of the list by at most 1.

        EXAMPLES::

            sage: from sage.data_structures.list_of_pairs import ListOfPairs
            sage: l = ListOfPairs()
            sage: l[0] = (2, 1)
            sage: l[1] = (1, 2)
            sage: l[0]
            (2, 1)
            sage: l[1]
            (1, 2)
            sage: l[10] = (5, 3)
            Traceback (most recent call last):
            ...
            IndexError
            sage: l[2] = 2
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        cdef size_t first, second
        (first, second) = value

        if index == self.length:
            self.add(first, second)
            return

        cdef pair_s* pair_pt = self.get(index)
        pair_pt.first = first
        pair_pt.second = second

    cdef inline int add(self, size_t first, size_t second) except -1:
        """
        Add a pair to the list.
        """
        self.enlarge()
        cdef pair_s* last = self.get(self.length - 1)
        last.first = first
        last.second = second
