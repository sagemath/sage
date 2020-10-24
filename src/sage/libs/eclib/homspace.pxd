from ..eclib cimport homspace

cdef class ModularSymbols:
    cdef homspace* H
