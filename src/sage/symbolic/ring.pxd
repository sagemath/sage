# sage_setup: distribution = sagemath-symbolics
cimport sage.rings.abc

cdef class SymbolicRing(sage.rings.abc.SymbolicRing):
    cdef public dict symbols
