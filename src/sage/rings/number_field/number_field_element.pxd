# sage_setup: distribution = sagemath-flint
cimport sage.structure.element
from sage.libs.gmp.types cimport mpz_t
from sage.rings.integer cimport Integer
from sage.rings.number_field.number_field_element_base cimport NumberFieldElement_base
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.parent cimport Parent
from sage.libs.ntl.types cimport ZZ_c, ZZX_c
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ


cdef class NumberFieldElement(NumberFieldElement_base):
    cdef ZZX_c _numerator
    cdef ZZ_c _denominator
    # Pointers to the defining polynomial (with numerator) for the field.
    # I keep these as pointers for arithmetic speed.
    cdef ntl_ZZX _fld_numerator
    cdef ntl_ZZ _fld_denominator
    cdef object __multiplicative_order
    cdef object __pari
    cdef object __matrix

    cdef _new(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

    cpdef _add_(self, other)
    cpdef _mul_(self, other)

    cpdef _copy_for_parent(self, Parent parent)

    cdef number_field(self)

    cdef void _ntl_coeff_as_mpz(self, mpz_t z, long i) noexcept
    cdef void _ntl_denom_as_mpz(self, mpz_t z) noexcept

    cdef void _reduce_c_(self) noexcept

    cpdef list _coefficients(self)

    cpdef bint is_rational(self) noexcept
    cpdef bint is_one(self) noexcept
    cdef int _randomize(self, num_bound, den_bound, distribution) except -1


cdef class NumberFieldElement_absolute(NumberFieldElement):
    pass

cdef class NumberFieldElement_relative(NumberFieldElement):
    pass

# TODO: cyclotomic and/or quadratic classes? (Both for differing implementations and speed).

cdef class OrderElement_absolute(NumberFieldElement_absolute):
    cdef object _number_field

cdef class OrderElement_relative(NumberFieldElement_relative):
    cdef object _number_field
