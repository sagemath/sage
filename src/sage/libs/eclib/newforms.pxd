# sage_setup: distribution = sagemath-eclib
from sage.libs.eclib cimport newforms

cdef class ECModularSymbol:
    cdef newforms* nfs
    cdef int n
    cdef object _E
    cdef int sign
