from sage.structure.element cimport Element
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.polynomial.polydict cimport PolyDict
from sage.rings.polynomial.polydict cimport ETuple

from sage.rings.padics.padic_generic_element cimport pAdicGenericElement


cdef class TateAlgebraTerm(MonoidElement):
    cdef _field
    cdef pAdicGenericElement _coeff
    cdef ETuple _exponent

    cpdef _mul_(self, other) noexcept
    cdef TateAlgebraTerm _floordiv_c(self, TateAlgebraTerm other) noexcept
    cpdef _floordiv_(self, other) noexcept

    cdef TateAlgebraTerm _new_c(self) noexcept
    cdef long _valuation_c(self) noexcept
    cdef long _cmp_c(self, TateAlgebraTerm other) except? 300
    cdef Element _call_c(self, list arg) noexcept
    cpdef TateAlgebraTerm monomial(self) noexcept
    cpdef TateAlgebraTerm monic(self) noexcept
    cdef TateAlgebraTerm _gcd_c(self, TateAlgebraTerm other) noexcept
    cdef TateAlgebraTerm _lcm_c(self, TateAlgebraTerm other) noexcept
    cdef bint _divides_c(self, TateAlgebraTerm other, bint integral) noexcept


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    cdef _prec
    cdef PolyDict _poly
    cdef list _terms
    cdef list _terms_nonzero
    cdef bint _is_normalized

    cdef _normalize(self) noexcept
    cdef TateAlgebraElement _new_c(self) noexcept
    cdef list _terms_c(self, bint include_zero=*) noexcept
    cpdef valuation(self) noexcept
    cdef TateAlgebraElement _term_mul_c(self, TateAlgebraTerm term) noexcept
    cdef TateAlgebraElement _positive_lshift_c(self, n) noexcept
    cdef TateAlgebraElement _lshift_c(self, n) noexcept
    cpdef TateAlgebraElement monic(self) noexcept
    cdef _quo_rem_c(self, list divisors, bint quo, bint rem, bint integral) noexcept
    cdef _quo_rem_check(self, divisors, bint quo, bint rem) noexcept
    cdef TateAlgebraElement _Spoly_c(self, TateAlgebraElement other) noexcept

