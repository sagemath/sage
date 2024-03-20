# sage_setup: distribution = sagemath-flint

from sage.rings.rational cimport Rational
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.modules.vector_real_double_dense cimport Vector_real_double_dense
from sage.rings.real_mpfi cimport RealIntervalFieldElement


cdef class interval_bernstein_polynomial:
    cdef int min_variations
    cdef int max_variations

    cdef Rational lower
    cdef Rational upper

    # 1: positive; -1: negative; 0: unknown (note that 0 does NOT mean zero)
    cdef int lsign
    cdef int usign
    cdef int level
    cdef lft
    cdef RealIntervalFieldElement slope_err
    cdef int bitsize
    cdef int scale_log2

    cdef void update_variations(self, interval_bernstein_polynomial bp1, interval_bernstein_polynomial bp2) noexcept
    cdef int degree(self) noexcept

cdef class interval_bernstein_polynomial_integer(interval_bernstein_polynomial):
    cdef Vector_integer_dense coeffs

    cdef int error

    cdef void _set_bitsize(self) noexcept
    cdef void _count_variations(self) noexcept

cdef class interval_bernstein_polynomial_float(interval_bernstein_polynomial):
    cdef Vector_real_double_dense coeffs

    cdef double neg_err
    cdef double pos_err

    cdef void _count_variations(self) noexcept

# forward declaration
cdef class rr_gap

cdef class island:
    cdef interval_bernstein_polynomial bp
    cdef ancestors
    cdef target_width
    cdef rr_gap lgap
    cdef rr_gap rgap
    cdef known_done

cdef class rr_gap:
    cdef Rational lower
    cdef Rational upper
    cdef int sign
    cdef island lisland
    cdef island risland


cdef class context:
    cdef random
    cdef int do_logging
    cdef int wordsize
    cdef seed # for debug printing
    cdef dc_log
    cdef be_log

    cdef void dc_log_append(self, x) noexcept
    cdef void be_log_append(self, x) noexcept

cdef class ocean:
    cdef context ctx
    cdef bpf
    cdef mapping
    cdef island endpoint
    cdef rr_gap lgap
    cdef rr_gap rgap
    cdef int msb
    cdef int prec
