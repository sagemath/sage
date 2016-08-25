from sage.libs.gmp.types cimport mpq_t

cimport sage.structure.element
cimport sage.rings.integer as integer

cpdef rational_power_parts(a, b, factor_limit=?)

cdef class Rational(sage.structure.element.FieldElement):
    cdef mpq_t value

    cdef __set_value(self, x, unsigned int base)
    cdef void set_from_mpq(Rational self, mpq_t value)
    cdef _lshift(self, long int exp)
    cdef _rshift(self, long int exp)

    cdef _val_unit(self, integer.Integer p)
