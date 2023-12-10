from sage.structure.element import Element
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from sage.rings.integer cimport Integer
from sage.rings.polynomial.commutative_polynomial cimport CommutativePolynomial
from sage.rings.polynomial.polynomial_compiled cimport CompiledPolynomialFunction


cdef class Polynomial(CommutativePolynomial):
    cdef Polynomial _new_generic(self, list coeffs) noexcept
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled
    cpdef Polynomial truncate(self, long n) noexcept
    cpdef Polynomial inverse_series_trunc(self, long prec) noexcept
    cdef long _hash_c(self) except -1
    cpdef constant_coefficient(self) noexcept
    cpdef Polynomial _new_constant_poly(self, a, Parent P) noexcept
    cpdef list list(self, bint copy=*) noexcept
    cpdef _mul_generic(self, right) noexcept
    cdef _square_generic(self) noexcept

    cpdef bint is_zero(self) except -1
    cpdef bint is_one(self) except -1
    cpdef bint is_term(self) except -1

    cpdef dict _mpoly_dict_recursive(self, tuple variables=*, base_ring=*) noexcept

    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef _floordiv_(self, right) noexcept
    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n) noexcept
    cpdef Polynomial _power_trunc(self, unsigned long n, long prec) noexcept
    cdef Polynomial _mul_term(self, Polynomial term, bint term_on_right) noexcept

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n) noexcept

    cdef get_coeff_c(self, Py_ssize_t i) noexcept
    cdef get_unsafe(self, Py_ssize_t i) noexcept
    cpdef long number_of_terms(self) noexcept

    # See 23227
    cpdef _add_(self, right) noexcept
    cpdef _mul_(self, right) noexcept
    cpdef _floordiv_(self, right) noexcept

    cdef public dict _cached_methods

cdef class Polynomial_generic_dense(Polynomial):
    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P) noexcept
    cdef list _coeffs
    cdef int _normalize(self) except -1
    cpdef list list(self, bint copy=*) noexcept

cdef class Polynomial_generic_dense_inexact(Polynomial_generic_dense):
    pass

cpdef is_Polynomial(f) noexcept
cpdef Polynomial generic_power_trunc(Polynomial p, Integer n, long prec) noexcept
cpdef list _dict_to_list(dict x, zero) noexcept

cpdef bint polynomial_is_variable(x) noexcept

