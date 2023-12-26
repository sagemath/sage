# sage_setup: distribution = sagemath-eclib
from sage.libs.eclib cimport homspace

cdef class ModularSymbols:
    cdef homspace* H
