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

from sage.libs.gap.gap_includes cimport Obj

############################################################################
### Hooking into the GAP memory management #################################
############################################################################

cdef class ObjWrapper():
    cdef Obj value

cdef ObjWrapper wrap_obj(Obj obj)

# returns the refcount dictionary for debugging purposes
cpdef get_owned_objects()

# Reference count GAP objects that you want to prevent from being
# garbage collected
cdef void reference_obj(Obj obj) noexcept
cdef void dereference_obj(Obj obj) noexcept

# callback from the GAP memory manager so we can mark all_gap_elements.values()
cdef void gasman_callback() noexcept with gil


############################################################################
### Initialization of GAP ##################################################
############################################################################

cdef initialize()


############################################################################
### Evaluate string in GAP #################################################
############################################################################

# Evaluate a string
cdef Obj gap_eval(str gap_string) except? NULL
