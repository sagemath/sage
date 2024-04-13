# sage_setup: distribution = sagemath-categories
#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cpdef fibers(f, domain)

from sage.structure.parent cimport Parent
from sage.structure.list_clone cimport ClonableIntArray

cdef class FiniteSetMap_MN(ClonableIntArray):
    cpdef _setimage(self, int i, int j)
    cpdef _getimage(self, int i)
    cpdef setimage(self, i, j)
    cpdef getimage(self, i)
    cpdef domain(self)
    cpdef codomain(self)
    cpdef image_set(self)
    cpdef fibers(self)
    cpdef items(self)
    cpdef FiniteSetMap_MN _compose_internal_(self, FiniteSetMap_MN other,
                                             Parent resParent)
    cpdef check(self)

cdef class FiniteSetMap_Set(FiniteSetMap_MN): pass

cpdef FiniteSetMap_Set FiniteSetMap_Set_from_list(cls, parent, lst)
cpdef FiniteSetMap_Set FiniteSetMap_Set_from_dict(cls, parent, d)

cdef class FiniteSetEndoMap_N(FiniteSetMap_MN): pass
cdef class FiniteSetEndoMap_Set(FiniteSetMap_Set): pass
