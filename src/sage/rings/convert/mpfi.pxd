# sage_setup: distribution = sagemath-flint

from sage.libs.mpfi.types cimport mpfi_ptr

cdef int mpfi_set_sage(mpfi_ptr re, mpfi_ptr im, x, field, int base) except -1
cdef int mpfi_interv_sage(mpfi_ptr re, mpfi_ptr im, x, y, field, int base) except -1
cdef int mpfi_set_via_RR(mpfi_ptr re, x, field) except -1
