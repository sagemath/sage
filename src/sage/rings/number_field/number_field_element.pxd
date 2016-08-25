cimport sage.structure.element
from sage.libs.gmp.types cimport mpz_t
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport FieldElement, RingElement, ModuleElement
from sage.structure.parent_base cimport ParentWithBase
from sage.libs.ntl.types cimport ZZ_c, ZZX_c
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ

cdef class NumberFieldElement(FieldElement):
    cdef ZZX_c __numerator
    cdef ZZ_c __denominator
    # Pointers to the defining polynomial (with numerator) for the field.
    # I keep these as pointers for arithmetic speed.
    cdef ntl_ZZX __fld_numerator
    cdef ntl_ZZ __fld_denominator
    cdef object __multiplicative_order
    cdef object __pari
    cdef object __matrix

    cdef _new(self)

    cdef number_field(self)

    cdef void _ntl_coeff_as_mpz(self, mpz_t z, long i)
    cdef void _ntl_denom_as_mpz(self, mpz_t z)

    cdef void _reduce_c_(self)

    cpdef list _coefficients(self)

    cpdef bint is_rational(self)
    cpdef bint is_one(self)
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
