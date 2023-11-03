#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element cimport Element


# Cython-0.17.2 disallows inline cpdef in non-final classes
# This restriction will be lifted at one point, then we can set
# some of the methods to be inline again, that is,
# revert the patch form https://github.com/sagemath/sage/issues/13740

cdef class ClonableElement(Element):
    cdef bint _is_immutable
    cdef bint _needs_check
    cdef long int  _hash

    cpdef bint _require_mutable(self) except -2
    cpdef bint is_mutable(self) noexcept
    cpdef bint is_immutable(self) noexcept
    cpdef set_immutable(self) noexcept

    cpdef _set_mutable(self) noexcept

    cpdef ClonableElement clone(self, bint check=?) noexcept

cdef class ClonableArray(ClonableElement):
    cdef list _list

    cpdef list _get_list(self) noexcept
    cpdef _set_list(self, list lst) noexcept
    cpdef ClonableArray __copy__(self) noexcept
    cpdef check(self) noexcept
    cpdef object _getitem(self, int key) noexcept
    cpdef _setitem(self, int key, value) noexcept
    cpdef int index(self, key, start=*, stop=*) except -1
    cpdef int count(self, key) except -1
    cpdef long int _hash_(self) except? -1

cdef class ClonableList(ClonableArray):
    cpdef append(self, el) noexcept
    cpdef extend(self, it) noexcept
    cpdef insert(self, int index, el) noexcept
    cpdef pop(self, int index=*) noexcept
    cpdef remove(self, el) noexcept

cdef class NormalizedClonableList(ClonableList):
    cpdef normalize(self) noexcept

cdef class ClonableIntArray(ClonableElement):
    cdef int _len
    cdef int* _list

    cpdef _alloc_(self, int size) noexcept
    cpdef ClonableIntArray __copy__(self) noexcept
    cpdef check(self) noexcept
    cpdef object _getitem(self, int key) noexcept
    cpdef _setitem(self, int item, value) noexcept
    cpdef int index(self, int item) except -1
    cpdef long int _hash_(self) except? -1
    cpdef list list(self) noexcept
