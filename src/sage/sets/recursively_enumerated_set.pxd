#*****************************************************************************
#       Copyright (C) 2014 Sage
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#
###############################################################################

cimport sage.structure.parent

cdef class RecursivelyEnumeratedSet_generic(sage.structure.parent.Parent):
    cdef readonly _seeds
    cdef public successors
    cdef readonly str _enumeration
    cdef readonly _max_depth
    cdef readonly _graded_component

    cpdef seeds(self)
    cpdef graded_component(self, depth)

cdef class RecursivelyEnumeratedSet_symmetric(RecursivelyEnumeratedSet_generic):
    cdef set _get_next_graded_component(self, set A, set B)

    cpdef graded_component(self, depth)

cdef class RecursivelyEnumeratedSet_graded(RecursivelyEnumeratedSet_generic):
    cdef set _get_next_graded_component(self, set B)

    cpdef graded_component(self, depth)
