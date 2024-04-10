# sage_setup: distribution = sagemath-pari
from sage.libs.gmp.types cimport mpz_t
from cypari2.gen cimport Gen as pari_gen
from sage.rings.padics.padic_capped_relative_element cimport CRElement

ctypedef mpz_t celement
include "CA_template_header.pxi"

cdef class pAdicCappedAbsoluteElement(CAElement):
    cdef lift_c(self)
    cdef pari_gen _to_gen(self)

from sage.rings.padics.pow_computer cimport PowComputer_base
cdef class PowComputer_(PowComputer_base):
    pass
