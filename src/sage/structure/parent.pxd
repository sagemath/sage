# sage_setup: distribution = sagemath-objects
#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport sage.structure.category_object
from sage.structure.coerce_dict cimport MonoDict, TripleDict

cdef class Parent(sage.structure.category_object.CategoryObject):
    cdef _element_constructor
    cdef bint _element_init_pass_parent
    cdef public _convert_method_name
    cdef public _initial_coerce_list
    cdef public _initial_action_list
    cdef public _initial_convert_list
    cdef readonly bint _coercions_used

    # Flags, see below
    cdef int flags
    cdef inline bint get_flag(self, int flag) noexcept:
        return self.flags & flag

    cpdef register_coercion(self, mor)
    cpdef register_action(self, action)
    cpdef register_conversion(self, mor)
    cpdef register_embedding(self, embedding)

    cpdef bint is_exact(self) except -2

    # Called from the __init__ method to set up coercion.
    cdef int init_coerce(self, bint warn=*) except -1

    # returns whether or not there is a Morphism from S to self
    cpdef bint has_coerce_map_from(self, S) except -2

    # returns a Morphism from S to self, or None
    cpdef coerce_map_from(self, S)
    cpdef _internal_coerce_map_from(self, S)
    cpdef _coerce_map_from_(self, S)

    # returns a Map from S to self, or None
    cpdef convert_map_from(self, S)
    cpdef _internal_convert_map_from(self, S)
    cpdef _convert_map_from_(self, S)
    cdef convert_method_map(self, S, method_name)

    # returns the Action by/on self on/by S
    # corresponding to op and self_on_left
    cpdef get_action(self, S, op=*, bint self_on_left=*, self_el=*, S_el=*)
    cpdef _get_action_(self, S, op, bint self_on_left)

    # coerce x into self
    cpdef coerce(self, x)

    cpdef an_element(self)
    cdef public object _cache_an_element

    # For internal use
    cpdef _generic_convert_map(self, S, category=*)
    cpdef _generic_coerce_map(self, S)
    cdef discover_coerce_map_from(self, S)
    cdef discover_convert_map_from(self, S)
    cdef discover_action(self, S, op, bint self_on_left, self_el=*, S_el=*)

    # List consisting of Morphisms (from anything to self)
    # and Parents for which the __call__ method of self
    # results in natural coercion.
    # Initialized at ring creation.
    cdef list _coerce_from_list
    # List of the domains of the registered coercions, to make
    # sure that the maps in _coerce_from_list remain valid.
    # This is important, since they are fundamental for discovering
    # new coercions by backtracking.
    cdef list _registered_domains
    # Hashtable of everything we've (possibly recursively) discovered so far.
    cdef MonoDict _coerce_from_hash

    # List consisting of Actions (either by or on self)
    # and Parents for which self._rmul_ and/or self._lmul_
    # do the correct thing.
    # Initialized at ring creation.
    cdef list _action_list

    # List consisting of Morphisms (from anything to self)
    # and Parents for which the __call__ method of self
    # does not result in type errors
    # Initialized at ring creation.
    cdef list _convert_from_list
    # Hashtable of everything we've (possibly recursively) discovered so far.
    cdef MonoDict _convert_from_hash
    # An optional single Morphism that describes a canonical coercion out of self
    cdef _embedding

    # Write-only hashtable of all actions discovered using this parent.
    # This is only needed to keep a strong reference to actions, to
    # prevent them being garbage collected prematurely.
    cdef TripleDict _action_hash


cdef class Set_generic(Parent):
    pass


# Flags for Parent.flags
cdef enum:
    # If this flag is set, call __richcmp__ on elements without
    # coercion. This allows a completely custom comparison function.
    Parent_richcmp_element_without_coercion = 1
