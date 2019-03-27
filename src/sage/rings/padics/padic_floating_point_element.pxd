from sage.libs.gmp.types cimport mpz_t
from cypari2.gen cimport Gen as pari_gen

ctypedef mpz_t celement
include "FP_template_header.pxi"

cdef class pAdicFloatingPointElement(FPElement):
    cdef lift_c(self)
    cdef pari_gen _to_gen(self)

from sage.rings.padics.pow_computer cimport PowComputer_base
cdef class PowComputer_(PowComputer_base):
    pass
