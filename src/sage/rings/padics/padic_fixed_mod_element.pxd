from sage.libs.gmp.types cimport mpz_t
from sage.libs.pari.gen cimport gen as pari_gen

ctypedef mpz_t celement
include "FM_template_header.pxi"

cdef class pAdicFixedModElement(FMElement):
    cdef lift_c(self)
    cdef pari_gen _to_gen(self)
