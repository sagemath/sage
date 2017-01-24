from .types cimport GEN
from .gen cimport gen

cpdef integer_to_gen(x)
cpdef gen_to_integer(gen x)

cdef GEN gtoi(GEN g0) except NULL
cdef GEN PyLong_AsGEN(x)
cdef PyLong_FromGEN(GEN g)

cdef gen new_t_POL_from_int_star(int* vals, unsigned long length, long varnum)
cdef gen new_gen_from_double(double)
cdef gen new_t_COMPLEX_from_double(double re, double im)
