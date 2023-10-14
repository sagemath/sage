# distutils: libraries = flint
# distutils: depends = flint/fmpz_extras.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    slong fmpz_allocated_bytes(const fmpz_t x)
    # Returns the total number of bytes heap-allocated internally by this object.
    # The count excludes the size of the structure itself. Add
    # ``sizeof(fmpz)`` to get the size of the object as a whole.

    void fmpz_adiv_q_2exp(fmpz_t z, const fmpz_t x, flint_bitcnt_t exp)
    # Sets *z* to `x / 2^{exp}`, rounded away from zero.

    void fmpz_ui_mul_ui(fmpz_t x, ulong a, ulong b)
    # Sets *x* to *a* times *b*.

    void fmpz_max(fmpz_t z, const fmpz_t x, const fmpz_t y)

    void fmpz_min(fmpz_t z, const fmpz_t x, const fmpz_t y)
    # Sets *z* to the maximum (respectively minimum) of *x* and *y*.

    void fmpz_add_inline(fmpz_t z, const fmpz_t x, const fmpz_t y)

    void fmpz_add_si_inline(fmpz_t z, const fmpz_t x, slong y)

    void fmpz_add_ui_inline(fmpz_t z, const fmpz_t x, ulong y)
    # Sets *z* to the sum of *x* and *y*.

    void fmpz_sub_si_inline(fmpz_t z, const fmpz_t x, slong y)
    # Sets *z* to the difference of *x* and *y*.

    void fmpz_add2_fmpz_si_inline(fmpz_t z, const fmpz_t x, const fmpz_t y, slong c)
    # Sets *z* to the sum of *x*, *y*, and *c*.

    mp_size_t _fmpz_size(const fmpz_t x)
    # Returns the number of limbs required to represent *x*.

    slong _fmpz_sub_small(const fmpz_t x, const fmpz_t y)
    # Computes the difference of *x* and *y* and returns the result as
    # an *slong*. The result is clamped between -*WORD_MAX* and *WORD_MAX*,
    # i.e. between `\pm (2^{63}-1)` inclusive on a 64-bit machine.

    void _fmpz_set_si_small(fmpz_t x, slong v)
    # Sets *x* to the integer *v* which is required to be a value
    # between *COEFF_MIN* and *COEFF_MAX* so that promotion to
    # a bignum cannot occur.

    void fmpz_set_mpn_large(fmpz_t z, mp_srcptr src, mp_size_t n, int negative)
    # Sets *z* to the integer represented by the *n* limbs in the array *src*,
    # or minus this value if *negative* is 1.
    # Requires `n \ge 2` and that the top limb of *src* is nonzero.
    # Note that *fmpz_set_ui*, *fmpz_neg_ui* can be used for single-limb integers.

    void fmpz_lshift_mpn(fmpz_t z, mp_srcptr src, mp_size_t n, int negative, flint_bitcnt_t shift)
    # Sets *z* to the integer represented by the *n* limbs in the array *src*,
    # or minus this value if *negative* is 1, shifted left by *shift* bits.
    # Requires `n \ge 1` and that the top limb of *src* is nonzero.
