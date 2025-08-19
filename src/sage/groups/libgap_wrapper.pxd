from sage.structure.element cimport MultiplicativeGroupElement
from sage.libs.gap.element cimport GapElement

cdef class ElementLibGAP(MultiplicativeGroupElement):
    cdef GapElement _libgap
    cpdef GapElement gap(self)
    cpdef _mul_(self, other)
