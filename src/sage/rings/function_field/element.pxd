from sage.structure.element cimport FieldElement


cdef class FunctionFieldElement(FieldElement):
    cdef readonly object _x
    cdef readonly object _matrix

    cdef FunctionFieldElement _new_c(self)
    cpdef bint is_nth_power(self, n) noexcept
    cpdef FunctionFieldElement nth_root(self, n)
