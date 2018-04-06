from sage.structure.element cimport MultiplicativeGroupElement, Element, MonoidElement, Matrix
from sage.groups.libgap_wrapper cimport ElementLibGAP

cpdef is_MatrixGroupElement(x)

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    cdef public Matrix _matrix

    cpdef _act_on_(self, x, bint self_on_left)
    cpdef _mul_(self, other)
    cpdef list list(self)

cdef class MatrixGroupElement_gap(ElementLibGAP):
    cpdef _act_on_(self, x, bint self_on_left)
    cpdef list list(self)

