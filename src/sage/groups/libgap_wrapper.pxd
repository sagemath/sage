from sage.structure.element cimport MultiplicativeGroupElement, MonoidElement
from sage.libs.gap.element cimport GapElement


cdef class ElementLibGAP(MultiplicativeGroupElement):
    cdef GapElement _libgap
    cpdef GapElement gap(self)
