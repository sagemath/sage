# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport mat_GF2E_c
from sage.libs.ntl.ntl_GF2EContext cimport ntl_GF2EContext_class
from sage.libs.ntl.ntl_GF2E cimport ntl_GF2E

cdef class ntl_mat_GF2E():
    cdef mat_GF2E_c x
    cdef ntl_GF2EContext_class c
    cdef ntl_GF2E _new_element(self)
    cdef ntl_mat_GF2E _new(self)
