from sage.libs.gmp.types cimport mpq_t

cimport sage.structure.element
cimport sage.rings.integer as integer

cpdef rational_power_parts(a, Rational b, factor_limit=?) noexcept

cdef class Rational(sage.structure.element.FieldElement):
    cdef mpq_t value

    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef _pow_(self, other) noexcept
    cdef __set_value(self, x, unsigned int base) noexcept
    cdef void set_from_mpq(Rational self, mpq_t value) noexcept
    cdef _lshift(self, long int exp) noexcept
    cdef _rshift(self, long int exp) noexcept

    cdef _val_unit(self, integer.Integer p) noexcept
