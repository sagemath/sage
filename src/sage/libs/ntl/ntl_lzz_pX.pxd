# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.lzz_p cimport *
from sage.libs.ntl.lzz_pX cimport *

from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

cdef class ntl_zz_pX():
    cdef zz_pX_c x
    cdef ntl_zz_pContext_class c
    cdef ntl_zz_pX _new(self)
