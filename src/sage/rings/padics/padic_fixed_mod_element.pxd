from sage.libs.gmp.types cimport mpz_t
from cypari2.gen cimport Gen as pari_gen
from sage.rings.padics.padic_floating_point_element cimport FPElement

ctypedef mpz_t celement
include "FM_template_header.pxi"

cdef class pAdicFixedModElement(FMElement):
    cdef lift_c(self)
    cdef pari_gen _to_gen(self)

from sage.rings.padics.pow_computer cimport PowComputer_base
cdef class PowComputer_(PowComputer_base):
    pass
