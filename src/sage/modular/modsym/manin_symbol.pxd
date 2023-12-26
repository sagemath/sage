from sage.structure.element cimport Element

cdef class ManinSymbol(Element):
    cdef public i, u, v
