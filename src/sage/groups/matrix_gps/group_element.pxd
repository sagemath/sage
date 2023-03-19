from sage.structure.element cimport MultiplicativeGroupElement, Element, MonoidElement, Matrix

cpdef is_MatrixGroupElement(x)

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    cdef public Matrix _matrix

    cpdef _act_on_(self, x, bint self_on_left)
    cpdef _mul_(self, other)
    cpdef list list(self)
