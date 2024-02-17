#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cpdef fibers(f, domain) noexcept

from sage.structure.parent cimport Parent
from sage.structure.list_clone cimport ClonableIntArray

cdef class FiniteSetMap_MN(ClonableIntArray):
    cpdef _setimage(self, int i, int j) noexcept
    cpdef _getimage(self, int i) noexcept
    cpdef setimage(self, i, j) noexcept
    cpdef getimage(self, i) noexcept
    cpdef domain(self) noexcept
    cpdef codomain(self) noexcept
    cpdef image_set(self) noexcept
    cpdef fibers(self) noexcept
    cpdef items(self) noexcept
    cpdef FiniteSetMap_MN _compose_internal_(self, FiniteSetMap_MN other,
                                             Parent resParent) noexcept
    cpdef check(self) noexcept

cdef class FiniteSetMap_Set(FiniteSetMap_MN): pass

cpdef FiniteSetMap_Set FiniteSetMap_Set_from_list(cls, parent, lst) noexcept
cpdef FiniteSetMap_Set FiniteSetMap_Set_from_dict(cls, parent, d) noexcept

cdef class FiniteSetEndoMap_N(FiniteSetMap_MN): pass
cdef class FiniteSetEndoMap_Set(FiniteSetMap_Set): pass
