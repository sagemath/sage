# sage_setup: distribution = sagemath-ntl
from sage.rings.padics.padic_ext_element cimport pAdicExtElement
from sage.rings.padics.pow_computer_ext cimport PowComputer_ZZ_pX

cdef class pAdicZZpXElement(pAdicExtElement):
    cdef PowComputer_ZZ_pX prime_pow
