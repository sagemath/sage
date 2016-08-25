#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .util cimport *
from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element, ModuleElement, RingElement

cdef libGAP_Obj make_gap_list(sage_list) except NULL
cdef libGAP_Obj make_gap_record(sage_dict) except NULL
cdef libGAP_Obj make_gap_integer(sage_dict) except NULL
cdef libGAP_Obj make_gap_string(sage_string) except NULL

cdef GapElement make_any_gap_element(parent, libGAP_Obj obj)
cdef GapElement make_GapElement(parent, libGAP_Obj obj)
cdef GapElement_List make_GapElement_List(parent, libGAP_Obj obj)
cdef GapElement_Record make_GapElement_Record(parent, libGAP_Obj obj)
cdef GapElement_Integer make_GapElement_Integer(parent, libGAP_Obj obj)
cdef GapElement_Rational make_GapElement_Rational(parent, libGAP_Obj obj)
cdef GapElement_String make_GapElement_String(parent, libGAP_Obj obj)
cdef GapElement_Boolean make_GapElement_Boolean(parent, libGAP_Obj obj)
cdef GapElement_Function make_GapElement_Function(parent, libGAP_Obj obj)



cdef class GapElement(RingElement):

    # the pointer to the GAP object (memory managed by GASMAN)
    cdef libGAP_Obj value

    # comparison
    cdef bint _compare_by_id
    cdef bint _compare_equal(self, Element other) except -2
    cdef bint _compare_less(self, Element other) except -2
    cpdef _set_compare_by_id(self)
    cpdef _assert_compare_by_id(self)

    cdef _initialize(self, parent, libGAP_Obj obj)
    cpdef _type_number(self)
    cpdef is_bool(self)


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
    cpdef libGAP_UInt record_name_to_index(self, bytes py_name)

cdef class GapElement_RecordIterator(object):
    cdef GapElement_Record rec
    cdef libGAP_UInt i

cdef class GapElement_List(GapElement):
    pass

cdef class GapElement_Permutation(GapElement):
    pass
