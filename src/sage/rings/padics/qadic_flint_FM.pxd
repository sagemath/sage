# sage_setup: distribution = sagemath-flint

from cypari2.gen cimport Gen as pari_gen
from sage.libs.flint.types cimport fmpz_poly_t
from sage.rings.padics.qadic_flint_FP cimport FPElement
from sage.rings.padics.pow_computer_flint cimport PowComputer_flint_unram

cdef class PowComputer_(PowComputer_flint_unram):
    pass
ctypedef fmpz_poly_t celement

include "FM_template_header.pxi"

cdef class qAdicFixedModElement(FMElement):
    pass
