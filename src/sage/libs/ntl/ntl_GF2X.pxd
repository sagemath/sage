# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport GF2X_c

cdef class ntl_GF2X():
    cdef GF2X_c x
