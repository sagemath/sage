# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.ZZ_pEX cimport ZZ_pEX_c
from sage.libs.ntl.ntl_ZZ_pEContext cimport ZZ_pEContext_ptrs
from sage.rings.integer cimport Integer

ctypedef ZZ_pEX_c celement
ctypedef ZZ_pEContext_ptrs *cparent

include "polynomial_template_header.pxi"

cdef class Polynomial_ZZ_pEX(Polynomial_template):
    cdef _powmod_bigexp(Polynomial_ZZ_pEX self, Integer exp, Polynomial_ZZ_pEX modulus)
