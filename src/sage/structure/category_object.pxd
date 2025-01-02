# sage_setup: distribution = sagemath-objects
#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject

cpdef check_default_category(default_category, category)

cdef class CategoryObject(SageObject):
    cdef public dict _cached_methods
    cdef _category
    cdef public _base
    cdef public _names
    cdef public _factory_data
    cdef object __weakref__
    cdef long _hash_value

    cdef getattr_from_category(self, name)

cpdef normalize_names(Py_ssize_t ngens, names)
cpdef bint certify_names(names) except -1
