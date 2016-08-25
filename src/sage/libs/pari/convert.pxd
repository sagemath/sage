from .types cimport GEN
from .gen cimport gen

cpdef integer_to_gen(x)
cpdef gen_to_integer(gen x)

cdef GEN gtoi(GEN g0) except NULL
cdef GEN PyLong_AsGEN(x)
cdef PyLong_FromGEN(GEN g)
