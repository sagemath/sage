# sage_setup: distribution = sagemath-combinat
from sage.structure.element cimport Element

cdef class Spin(Element):
    cdef bint* _value
    cdef int _n
    cdef long _hash

    cdef Spin _new_c(self, bint* value)

cdef class Spin_crystal_type_B_element(Spin):
    cpdef Spin e(self, int i)
    cpdef Spin f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Spin_crystal_type_D_element(Spin):
    cpdef Spin e(self, int i)
    cpdef Spin f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept
