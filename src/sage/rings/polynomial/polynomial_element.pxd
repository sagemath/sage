from sage.structure.element import Element
from sage.structure.element cimport Element, CommutativeAlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from sage.rings.integer cimport Integer
from sage.rings.polynomial.commutative_polynomial cimport CommutativePolynomial
from sage.rings.polynomial.polynomial_compiled cimport CompiledPolynomialFunction


cdef class Polynomial(CommutativePolynomial):
    cdef Polynomial _new_generic(self, list coeffs)
    cdef char _is_gen
    cdef CompiledPolynomialFunction _compiled
    cpdef Polynomial truncate(self, long n)
    cpdef Polynomial inverse_series_trunc(self, long prec)
    cdef long _hash_c(self) except -1
    cpdef constant_coefficient(self)
    cpdef Polynomial _new_constant_poly(self, a, Parent P)
    cpdef list list(self, bint copy=*)
    cpdef _mul_generic(self, right)
    cdef _square_generic(self)

    cpdef bint is_zero(self) except -1
    cpdef bint is_one(self) except -1
    cpdef bint is_term(self) except -1

    cpdef dict _mpoly_dict_recursive(self, tuple variables=*, base_ring=*)

    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _floordiv_(self, right)
    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n)
    cpdef Polynomial _power_trunc(self, unsigned long n, long prec)
    cdef Polynomial _mul_term(self, Polynomial term, bint term_on_right)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long n)

    cdef get_coeff_c(self, Py_ssize_t i)
    cdef get_unsafe(self, Py_ssize_t i)
    cpdef long number_of_terms(self) noexcept

    # See 23227
    cpdef _add_(self, right)
    cpdef _mul_(self, right)
    cpdef _floordiv_(self, right)

    cdef public dict _cached_methods

cdef class Polynomial_generic_dense(Polynomial):
    cdef Polynomial_generic_dense _new_c(self, list coeffs, Parent P)
    cdef list _coeffs
    cdef int _normalize(self) except -1
    cpdef list list(self, bint copy=*)

cdef class Polynomial_generic_dense_inexact(Polynomial_generic_dense):
    pass

cpdef is_Polynomial(f)
cpdef Polynomial generic_power_trunc(Polynomial p, Integer n, long prec)
cpdef list _dict_to_list(dict x, zero)

cpdef bint polynomial_is_variable(x) noexcept
