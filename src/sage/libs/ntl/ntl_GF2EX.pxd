from sage.libs.ntl.types cimport GF2EX_c
from sage.libs.ntl.ntl_GF2EContext cimport ntl_GF2EContext_class
from sage.libs.ntl.ntl_GF2E cimport ntl_GF2E

cdef class ntl_GF2EX():
    cdef GF2EX_c x
    cdef ntl_GF2EContext_class c
    cdef ntl_GF2E _new_element(self)
    cdef ntl_GF2EX _new(self)
