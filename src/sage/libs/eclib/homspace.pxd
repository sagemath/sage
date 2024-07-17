from sage.libs.eclib cimport homspace

cdef class ModularSymbols:
    cdef homspace* H
