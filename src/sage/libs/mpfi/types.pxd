# sage_setup: distribution = sagemath-ntl
from sage.libs.mpfr.types cimport __mpfr_struct

cdef extern from "mpfi.h":
    ctypedef struct __mpfi_struct:
        __mpfr_struct left
        __mpfr_struct right
    ctypedef __mpfi_struct mpfi_t[1]
    ctypedef __mpfi_struct* mpfi_ptr
    ctypedef __mpfi_struct* mpfi_srcptr
