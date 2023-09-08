# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport mat_ZZ_c

cdef class ntl_mat_ZZ():
    cdef mat_ZZ_c x
    cdef long __nrows, __ncols
