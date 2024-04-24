from sage.libs.ntl.ZZ_pEX cimport ZZ_pEX_c
from sage.libs.ntl.ntl_ZZ_pEContext cimport ZZ_pEContext_ptrs

ctypedef ZZ_pEX_c celement
ctypedef ZZ_pEContext_ptrs *cparent

include "polynomial_template_header.pxi"

cdef class Polynomial_ZZ_pEX(Polynomial_template):
    pass

