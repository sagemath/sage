# sage_setup: distribution = sagemath-gap
#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gap.gap_includes cimport Obj, UInt
from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element, ModuleElement, RingElement

cdef Obj make_gap_list(sage_list) except NULL
cdef Obj make_gap_matrix(sage_list, gap_ring) except NULL
cdef Obj make_gap_record(sage_dict) except NULL
cdef Obj make_gap_integer(sage_int) except NULL
cdef Obj make_gap_string(sage_string) except NULL

cdef GapElement make_any_gap_element(parent, Obj obj)
cdef GapElement make_GapElement(parent, Obj obj)
cdef GapElement_List make_GapElement_List(parent, Obj obj)
cdef GapElement_Record make_GapElement_Record(parent, Obj obj)
cdef GapElement_Integer make_GapElement_Integer(parent, Obj obj)
cdef GapElement_Rational make_GapElement_Rational(parent, Obj obj)
cdef GapElement_String make_GapElement_String(parent, Obj obj)
cdef GapElement_Boolean make_GapElement_Boolean(parent, Obj obj)
cdef GapElement_Function make_GapElement_Function(parent, Obj obj)
cdef GapElement_Permutation make_GapElement_Permutation(parent, Obj obj)

cdef char *capture_stdout(Obj, Obj) noexcept
cdef char *gap_element_str(Obj) noexcept
cdef char *gap_element_repr(Obj) noexcept


cdef class GapElement(RingElement):

    # the pointer to the GAP object (memory managed by GASMAN)
    cdef Obj value

    # comparison
    cdef bint _compare_by_id
    cdef bint _compare_equal(self, Element other) except -2
    cdef bint _compare_less(self, Element other) except -2
    cpdef _set_compare_by_id(self)
    cpdef _assert_compare_by_id(self)

    cdef _initialize(self, parent, Obj obj)
    cpdef _type_number(self)
    cpdef is_bool(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _mod_(self, other)
    cpdef _pow_(self, other)

    cpdef GapElement deepcopy(self, bint mut)

cdef class GapElement_Integer(GapElement):
    pass

cdef class GapElement_Rational(GapElement):
    pass

cdef class GapElement_IntegerMod(GapElement):
    cpdef GapElement_Integer lift(self)

cdef class GapElement_FiniteField(GapElement):
    cpdef GapElement_Integer lift(self)

cdef class GapElement_Cyclotomic(GapElement):
    pass

cdef class GapElement_Ring(GapElement):
    pass

cdef class GapElement_String(GapElement):
    pass

cdef class GapElement_Boolean(GapElement):
    pass

cdef class GapElement_Function(GapElement):
    pass

cdef class GapElement_MethodProxy(GapElement_Function):
    cdef GapElement first_argument

cdef class GapElement_Record(GapElement):
    cpdef UInt record_name_to_index(self, name) noexcept

cdef class GapElement_RecordIterator():
    cdef GapElement_Record rec
    cdef UInt i

cdef class GapElement_List(GapElement):
    pass

cdef class GapElement_Permutation(GapElement):
    pass
