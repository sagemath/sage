from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial

from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

from sage.libs.ntl.ntl_lzz_p cimport ntl_zz_p
from sage.libs.ntl.ntl_lzz_pX cimport ntl_zz_pX
from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

from sage.rings.integer cimport Integer

from sage.libs.ntl.ZZ_pX cimport *
from sage.libs.ntl.lzz_pX cimport *

# Word-sized modulus variant

cdef class Polynomial_dense_mod_n(Polynomial):
    cdef object __poly
    cdef object __singular

cdef class Polynomial_dense_modn_ntl_zz(Polynomial_dense_mod_n):
    cdef zz_pX_c x
    cdef ntl_zz_pContext_class c
    cdef Polynomial_dense_modn_ntl_zz _new(self)
    cpdef _mod_(self, right)

# Large modulus variant

ctypedef ZZ_pX_c celement
ctypedef ZZ_pContext_c *cparent

include "polynomial_template_header.pxi"

cdef class Polynomial_ZZ_pX(Polynomial_template):
    cdef inline Polynomial_ZZ_pX _new(self)

cdef class Polynomial_dense_modn_ntl_ZZ(Polynomial_ZZ_pX):
    pass

cdef class Polynomial_dense_mod_p(Polynomial_ZZ_pX):
    pass
