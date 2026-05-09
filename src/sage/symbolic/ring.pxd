from sage.structure.parent cimport Parent

cdef class SymbolicRing(Parent):
    cdef public dict symbols
