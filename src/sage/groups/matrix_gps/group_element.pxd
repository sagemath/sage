from sage.structure.element cimport MultiplicativeGroupElement, Element, MonoidElement, Matrix

cpdef is_MatrixGroupElement(x) noexcept

cdef class MatrixGroupElement_generic(MultiplicativeGroupElement):
    cdef public Matrix _matrix

    cpdef _act_on_(self, x, bint self_on_left) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef list list(self) noexcept
