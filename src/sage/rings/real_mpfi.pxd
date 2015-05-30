from sage.libs.mpfi cimport *

cimport sage.rings.ring

cimport sage.structure.element
from sage.structure.element cimport RingElement

from rational import Rational
from rational cimport Rational

cimport real_mpfr

cdef class RealIntervalFieldElement(sage.structure.element.RingElement)  # forward decl

cdef class RealIntervalField_class(sage.rings.ring.Field):
    cdef int __prec
    cdef bint sci_not
    # Cache RealField instances for the lower, upper, and middle bounds.
    # These have the same precision as the interval field;
    # __lower_field rounds down, __upper_field rounds up.
    # These fields with their roundings are not used for computation
    # in this module, but they do affect the printing and the return
    # values of lower() and upper().  Consider a 3-bit
    # interval containing exactly the floating-point number 1.25.
    # In round-to-nearest or round-down, this prints as 1.2; in round-up,
    # this prints as 1.3.  The straightforward options, then, are to
    # print this interval as [1.2 ... 1.2] (which does not even contain
    # the true value, 1.25), or to print it as [1.2 ... 1.3] (which
    # gives the impression that the upper and lower bounds are not
    # equal, even though they really are).  Neither of these is very
    # satisfying, but I have chosen the latter for now.
    cdef real_mpfr.RealField_class __lower_field
    cdef real_mpfr.RealField_class __middle_field
    cdef real_mpfr.RealField_class __upper_field
    cdef RealIntervalFieldElement _new(self)


cdef class RealIntervalFieldElement(sage.structure.element.RingElement):
    cdef mpfi_t value
    cdef char init
    cdef RealIntervalFieldElement _new(self)

    cdef RealIntervalFieldElement abs(RealIntervalFieldElement self)
    cdef Rational _simplest_rational_helper(self)
    cpdef _str_question_style(self, int base, int error_digits, e, bint prefer_sci)
