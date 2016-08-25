#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .gap_includes cimport *

############################################################################
### Hooking into the GAP memory management #################################
############################################################################

cdef class ObjWrapper(object):
    cdef libGAP_Obj value

cdef ObjWrapper wrap_obj(libGAP_Obj obj)

# a dictionary to keep all GAP elements
cdef dict owned_objects_refcount

# returns the refcount dictionary for debugging purposes
cpdef get_owned_objects()

# Reference count GAP objects that you want to prevent from being
# garbage collected
cdef void reference_obj(libGAP_Obj obj)
cdef void dereference_obj(libGAP_Obj obj)

# callback from the GAP memory manager so we can mark all_gap_elements.values()
cdef void gasman_callback()


############################################################################
### Initialization of libGAP ###############################################
############################################################################

# To ensure that we call initialize_libgap only once.
cdef bint _gap_is_initialized = False
cdef initialize()


############################################################################
### Helper to protect temporary objects from deletion ######################
############################################################################

# Hold a reference (inside the GAP kernel) to obj so that it doesn't
# get deleted this works by assigning it to a global variable. This is
# very simple, but you can't use it to keep two objects alive. Be
# careful.
cdef libGAP_UInt reference_holder
cdef void hold_reference(libGAP_Obj obj)


############################################################################
### Evaluate string in GAP #################################################
############################################################################

# Evaluate a string
cdef libGAP_Obj gap_eval(str gap_string) except? NULL


############################################################################
### Debug functions ########################################################
############################################################################

# Return details of the GAP memory pool
cpdef memory_usage()
