# distutils: libraries = flint
# distutils: depends = flint/mpfr_vec.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    mpfr_ptr _mpfr_vec_init(slong len, flint_bitcnt_t prec)
    # Returns a vector of the given length of initialised ``mpfr``'s
    # with the given exact precision.

    void _mpfr_vec_clear(mpfr_ptr vec, slong len)
    # Clears the given vector.

    void _mpfr_vec_zero(mpfr_ptr vec, slong len)
    # Zeros the vector ``(vec, len)``.

    void _mpfr_vec_set(mpfr_ptr vec1, mpfr_srcptr vec2, slong len)
    # Copies the vector ``vec2`` of the given length into ``vec1``.
    # No check is made to ensure ``vec1`` and ``vec2`` are different.

    void _mpfr_vec_add(mpfr_ptr res, mpfr_srcptr vec1, mpfr_srcptr vec2, slong len)
    # Adds the given vectors of the given length together and stores the
    # result in ``res``.

    void _mpfr_vec_scalar_mul_mpfr(mpfr_ptr res, mpfr_srcptr vec, slong len, mpfr_t c)
    # Multiplies the vector with given length by the scalar `c` and
    # sets ``res`` to the result.

    void _mpfr_vec_scalar_mul_2exp(mpfr_ptr res, mpfr_srcptr vec, slong len, flint_bitcnt_t exp)
    # Multiplies the given vector of the given length by ``2^exp``.

    void _mpfr_vec_scalar_product(mpfr_t res, mpfr_srcptr vec1, mpfr_srcptr vec2, slong len)
