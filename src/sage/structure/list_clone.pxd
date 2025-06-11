# sage_setup: distribution = sagemath-objects
# ***************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************

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
    cpdef set_immutable(self)

    cpdef _set_mutable(self)

    cpdef ClonableElement clone(self, bint check=?)

cdef class ClonableArray(ClonableElement):
    cdef list _list

    cpdef list _get_list(self)
    cpdef _set_list(self, list lst)
    cpdef ClonableArray __copy__(self)
    cpdef check(self)
    cpdef object _getitem(self, int key)
    cpdef _setitem(self, int key, value)
    cpdef int index(self, key, start=*, stop=*) except -1
    cpdef int count(self, key) except -1
    cpdef long int _hash_(self) except? -1

cdef class ClonableList(ClonableArray):
    cpdef append(self, el)
    cpdef extend(self, it)
    cpdef insert(self, int index, el)
    cpdef pop(self, int index=*)
    cpdef remove(self, el)

cdef class NormalizedClonableList(ClonableList):
    cpdef normalize(self)

cdef class ClonableIntArray(ClonableElement):
    cdef int _len
    cdef int* _list

    cpdef _alloc_(self, int size)
    cpdef ClonableIntArray __copy__(self)
    cpdef check(self)
    cpdef object _getitem(self, int key)
    cpdef _setitem(self, int item, value)
    cpdef int index(self, int item) except -1
    cpdef long int _hash_(self) except? -1
    cpdef list list(self)
