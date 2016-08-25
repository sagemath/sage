#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport sage.structure.parent as parent
from sage.structure.coerce_dict cimport MonoDict


cdef class Parent(parent.Parent):

    # returns a Morphism from S to self, or None
    cpdef coerce_map_from_c(self, S)
    cdef coerce_map_from_c_impl(self, S)

    # returns the Action by/on self on/by S
    # corresponding to op and self_on_left
    cpdef get_action_c(self, S, op, bint self_on_left)
    cdef get_action_c_impl(self, S, op, bint self_on_left)



    cdef public MonoDict _has_coerce_map_from

    #########################################
    # Canonical Coercion Methods
    cpdef has_coerce_map_from_c(self, S)
    cdef has_coerce_map_from_c_impl(self, S)
    cpdef _coerce_c(self, x)
    cdef _coerce_c_impl(self, x)

    cdef _an_element_c_impl(self)
    cpdef _an_element_c(self)
