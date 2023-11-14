from sage.structure.element cimport CommutativeAlgebraElement, ModuleElement, RingElement, Element


cdef class LaurentPolynomial(CommutativeAlgebraElement):
    cdef LaurentPolynomial _new_c(self) noexcept
    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef _floordiv_(self, other) noexcept
    cpdef long number_of_terms(self) except -1
    cpdef dict dict(self) noexcept

cdef class LaurentPolynomial_univariate(LaurentPolynomial):
    cdef ModuleElement __u
    cdef long __n
    cpdef _normalize(self) noexcept
    cpdef _unsafe_mutate(self, i, value) noexcept

