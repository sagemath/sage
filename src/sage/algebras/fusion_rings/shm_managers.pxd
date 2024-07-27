# sage_setup: distribution = sagemath-singular
cimport numpy as np
from sage.rings.number_field.number_field_base cimport NumberField
from sage.rings.number_field.number_field_element cimport NumberFieldElement_absolute

cdef class KSHandler:
    cdef list obj_cache
    cdef np.ndarray ks_dat
    cdef NumberField field
    cdef public object shm

    cdef bint contains(self, int idx) noexcept
    cdef NumberFieldElement_absolute get(self, int idx)
    cdef setitem(self, int idx, rhs)
    cpdef update(self, list eqns)

cdef class FvarsHandler:
    cdef dict sext_to_idx, obj_cache
    cdef int bytes, ngens
    cdef np.ndarray fvars
    cdef NumberField field
    cdef object fvars_t, pid_list
    cdef Py_ssize_t child_id
    cdef public object shm
