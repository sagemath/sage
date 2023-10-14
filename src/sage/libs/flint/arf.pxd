# distutils: libraries = flint
# distutils: depends = flint/arf.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void arf_init(arf_t x)
    # Initializes the variable *x* for use. Its value is set to zero.

    void arf_clear(arf_t x)
    # Clears the variable *x*, freeing or recycling its allocated memory.

    slong arf_allocated_bytes(const arf_t x)
    # Returns the total number of bytes heap-allocated internally by this object.
    # The count excludes the size of the structure itself. Add
    # ``sizeof(arf_struct)`` to get the size of the object as a whole.

    void arf_zero(arf_t res)

    void arf_one(arf_t res)

    void arf_pos_inf(arf_t res)

    void arf_neg_inf(arf_t res)

    void arf_nan(arf_t res)
    # Sets *res* respectively to 0, 1, `+\infty`, `-\infty`, NaN.

    bint arf_is_zero(const arf_t x)

    bint arf_is_one(const arf_t x)

    bint arf_is_pos_inf(const arf_t x)

    bint arf_is_neg_inf(const arf_t x)

    bint arf_is_nan(const arf_t x)
    # Returns nonzero iff *x* respectively equals 0, 1, `+\infty`, `-\infty`, NaN.

    bint arf_is_inf(const arf_t x)
    # Returns nonzero iff *x* equals either `+\infty` or `-\infty`.

    bint arf_is_normal(const arf_t x)
    # Returns nonzero iff *x* is a finite, nonzero floating-point value, i.e.
    # not one of the special values 0, `+\infty`, `-\infty`, NaN.

    bint arf_is_special(const arf_t x)
    # Returns nonzero iff *x* is one of the special values
    # 0, `+\infty`, `-\infty`, NaN, i.e. not a finite, nonzero
    # floating-point value.

    bint arf_is_finite(const arf_t x)
    # Returns nonzero iff *x* is a finite floating-point value,
    # i.e. not one of the values `+\infty`, `-\infty`, NaN.
    # (Note that this is not equivalent to the negation of
    # :func:`arf_is_inf`.)

    void arf_set(arf_t res, const arf_t x)

    void arf_set_mpz(arf_t res, const mpz_t x)

    void arf_set_fmpz(arf_t res, const fmpz_t x)

    void arf_set_ui(arf_t res, ulong x)

    void arf_set_si(arf_t res, slong x)

    void arf_set_mpfr(arf_t res, const mpfr_t x)

    void arf_set_d(arf_t res, double x)
    # Sets *res* to the exact value of *x*.

    void arf_swap(arf_t x, arf_t y)
    # Swaps *x* and *y* efficiently.

    void arf_init_set_ui(arf_t res, ulong x)

    void arf_init_set_si(arf_t res, slong x)
    # Initializes *res* and sets it to *x* in a single operation.

    int arf_set_round(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

    int arf_set_round_si(arf_t res, slong x, slong prec, arf_rnd_t rnd)

    int arf_set_round_ui(arf_t res, ulong x, slong prec, arf_rnd_t rnd)

    int arf_set_round_mpz(arf_t res, const mpz_t x, slong prec, arf_rnd_t rnd)

    int arf_set_round_fmpz(arf_t res, const fmpz_t x, slong prec, arf_rnd_t rnd)
    # Sets *res* to *x*, rounded to *prec* bits in the direction
    # specified by *rnd*.

    void arf_set_si_2exp_si(arf_t res, slong m, slong e)

    void arf_set_ui_2exp_si(arf_t res, ulong m, slong e)

    void arf_set_fmpz_2exp(arf_t res, const fmpz_t m, const fmpz_t e)
    # Sets *res* to `m \cdot 2^e`.

    int arf_set_round_fmpz_2exp(arf_t res, const fmpz_t x, const fmpz_t e, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x \cdot 2^e`, rounded to *prec* bits in the direction
    # specified by *rnd*.

    void arf_get_fmpz_2exp(fmpz_t m, fmpz_t e, const arf_t x)
    # Sets *m* and *e* to the unique integers such that
    # `x = m \cdot 2^e` and *m* is odd,
    # provided that *x* is a nonzero finite fraction.
    # If *x* is zero, both *m* and *e* are set to zero. If *x* is
    # infinite or NaN, the result is undefined.

    void arf_frexp(arf_t m, fmpz_t e, const arf_t x)
    # Writes *x* as `m \cdot 2^e`, where `0.5 \le |m| < 1` if *x* is a normal
    # value. If *x* is a special value, copies this to *m* and sets *e* to zero.
    # Note: for the inverse operation (*ldexp*), use :func:`arf_mul_2exp_fmpz`.

    double arf_get_d(const arf_t x, arf_rnd_t rnd)
    # Returns *x* rounded to a double in the direction specified by *rnd*.
    # This method rounds correctly when overflowing or underflowing
    # the double exponent range (this was not the case in an earlier version).

    int arf_get_mpfr(mpfr_t res, const arf_t x, mpfr_rnd_t rnd)
    # Sets the MPFR variable *res* to the value of *x*. If the precision of *x*
    # is too small to allow *res* to be represented exactly, it is rounded in
    # the specified MPFR rounding mode. The return value (-1, 0 or 1)
    # indicates the direction of rounding, following the convention
    # of the MPFR library.
    # If *x* has an exponent too large or small to fit in the MPFR type, the
    # result overflows to an infinity or underflows to a (signed) zero,
    # and the corresponding MPFR exception flags are set.

    int arf_get_fmpz(fmpz_t res, const arf_t x, arf_rnd_t rnd)
    # Sets *res* to *x* rounded to the nearest integer in the direction
    # specified by *rnd*. If rnd is *ARF_RND_NEAR*, rounds to the nearest
    # even integer in case of a tie. Returns inexact (beware: accordingly
    # returns whether *x* is *not* an integer).
    # This method aborts if *x* is infinite or NaN, or if the exponent of *x*
    # is so large that allocating memory for the result fails.
    # Warning: this method will allocate a huge amount of memory to store
    # the result if the exponent of *x* is huge. Memory allocation could
    # succeed even if the required space is far larger than the physical
    # memory available on the machine, resulting in swapping. It is recommended
    # to check that *x* is within a reasonable range before calling this method.

    slong arf_get_si(const arf_t x, arf_rnd_t rnd)
    # Returns *x* rounded to the nearest integer in the direction specified by
    # *rnd*. If *rnd* is *ARF_RND_NEAR*, rounds to the nearest even integer
    # in case of a tie. Aborts if *x* is infinite, NaN, or the value is
    # too large to fit in a slong.

    int arf_get_fmpz_fixed_fmpz(fmpz_t res, const arf_t x, const fmpz_t e)

    int arf_get_fmpz_fixed_si(fmpz_t res, const arf_t x, slong e)
    # Converts *x* to a mantissa with predetermined exponent, i.e. sets *res* to
    # an integer *y* such that `y \times 2^e \approx x`, truncating if necessary.
    # Returns 0 if exact and 1 if truncation occurred.
    # The warnings for :func:`arf_get_fmpz` apply.

    void arf_floor(arf_t res, const arf_t x)

    void arf_ceil(arf_t res, const arf_t x)
    # Sets *res* to `\lfloor x \rfloor` and `\lceil x \rceil` respectively.
    # The result is always represented exactly, requiring no more bits to
    # store than the input. To round the result to a floating-point number
    # with a lower precision, call :func:`arf_set_round` afterwards.

    void arf_get_fmpq(fmpq_t res, const arf_t x)
    # Set *res* to the exact rational value of *x*.
    # This method aborts if *x* is infinite or NaN, or if the exponent of *x*
    # is so large that allocating memory for the result fails.

    bint arf_equal(const arf_t x, const arf_t y)
    bint arf_equal_si(const arf_t x, slong y)
    bint arf_equal_ui(const arf_t x, ulong y)
    bint arf_equal_d(const arf_t x, double y)
    # Returns nonzero iff *x* and *y* are exactly equal. NaN is not
    # treated specially, i.e. NaN compares as equal to itself.
    # For comparison with a *double*, the values -0 and +0 are
    # both treated as zero, and all NaN values are treated as identical.

    int arf_cmp(const arf_t x, const arf_t y)

    int arf_cmp_si(const arf_t x, slong y)

    int arf_cmp_ui(const arf_t x, ulong y)

    int arf_cmp_d(const arf_t x, double y)
    # Returns negative, zero, or positive, depending on whether *x* is
    # respectively smaller, equal, or greater compared to *y*.
    # Comparison with NaN is undefined.

    int arf_cmpabs(const arf_t x, const arf_t y)

    int arf_cmpabs_ui(const arf_t x, ulong y)

    int arf_cmpabs_d(const arf_t x, double y)

    int arf_cmpabs_mag(const arf_t x, const mag_t y)
    # Compares the absolute values of *x* and *y*.

    int arf_cmp_2exp_si(const arf_t x, slong e)

    int arf_cmpabs_2exp_si(const arf_t x, slong e)
    # Compares *x* (respectively its absolute value) with `2^e`.

    int arf_sgn(const arf_t x)
    # Returns `-1`, `0` or `+1` according to the sign of *x*. The sign
    # of NaN is undefined.

    void arf_min(arf_t res, const arf_t a, const arf_t b)

    void arf_max(arf_t res, const arf_t a, const arf_t b)
    # Sets *res* respectively to the minimum and the maximum of *a* and *b*.

    slong arf_bits(const arf_t x)
    # Returns the number of bits needed to represent the absolute value
    # of the mantissa of *x*, i.e. the minimum precision sufficient to represent
    # *x* exactly. Returns 0 if *x* is a special value.

    bint arf_is_int(const arf_t x)
    # Returns nonzero iff *x* is integer-valued.

    bint arf_is_int_2exp_si(const arf_t x, slong e)
    # Returns nonzero iff *x* equals `n 2^e` for some integer *n*.

    void arf_abs_bound_lt_2exp_fmpz(fmpz_t res, const arf_t x)
    # Sets *res* to the smallest integer *b* such that `|x| < 2^b`.
    # If *x* is zero, infinity or NaN, the result is undefined.

    void arf_abs_bound_le_2exp_fmpz(fmpz_t res, const arf_t x)
    # Sets *res* to the smallest integer *b* such that `|x| \le 2^b`.
    # If *x* is zero, infinity or NaN, the result is undefined.

    slong arf_abs_bound_lt_2exp_si(const arf_t x)
    # Returns the smallest integer *b* such that `|x| < 2^b`, clamping
    # the result to lie between -*ARF_PREC_EXACT* and *ARF_PREC_EXACT*
    # inclusive. If *x* is zero, -*ARF_PREC_EXACT* is returned,
    # and if *x* is infinity or NaN, *ARF_PREC_EXACT* is returned.

    void arf_get_mag(mag_t res, const arf_t x)
    # Sets *res* to an upper bound for the absolute value of *x*.

    void arf_get_mag_lower(mag_t res, const arf_t x)
    # Sets *res* to a lower bound for the absolute value of *x*.

    void arf_set_mag(arf_t res, const mag_t x)
    # Sets *res* to *x*. This operation is exact.

    void mag_init_set_arf(mag_t res, const arf_t x)
    # Initializes *res* and sets it to an upper bound for *x*.

    void mag_fast_init_set_arf(mag_t res, const arf_t x)
    # Initializes *res* and sets it to an upper bound for *x*.
    # Assumes that the exponent of *res* is small (this function is unsafe).

    void arf_mag_set_ulp(mag_t res, const arf_t x, slong prec)
    # Sets *res* to the magnitude of the unit in the last place (ulp) of *x*
    # at precision *prec*.

    void arf_mag_add_ulp(mag_t res, const mag_t x, const arf_t y, slong prec)
    # Sets *res* to an upper bound for the sum of *x* and the
    # magnitude of the unit in the last place (ulp) of *y*
    # at precision *prec*.

    void arf_mag_fast_add_ulp(mag_t res, const mag_t x, const arf_t y, slong prec)
    # Sets *res* to an upper bound for the sum of *x* and the
    # magnitude of the unit in the last place (ulp) of *y*
    # at precision *prec*. Assumes that all exponents are small.

    void arf_init_set_shallow(arf_t z, const arf_t x)

    void arf_init_set_mag_shallow(arf_t z, const mag_t x)
    # Initializes *z* to a shallow copy of *x*. A shallow copy just involves
    # copying struct data (no heap allocation is performed).
    # The target variable *z* may not be cleared or modified in any way (it can
    # only be used as constant input to functions), and may not be used after
    # *x* has been cleared. Moreover, after *x* has been assigned shallowly
    # to *z*, no modification of *x* is permitted as slong as *z* is in use.

    void arf_init_neg_shallow(arf_t z, const arf_t x)

    void arf_init_neg_mag_shallow(arf_t z, const mag_t x)
    # Initializes *z* shallowly to the negation of *x*.

    void arf_randtest(arf_t res, flint_rand_t state, slong bits, slong mag_bits)
    # Generates a finite random number whose mantissa has precision at most
    # *bits* and whose exponent has at most *mag_bits* bits. The
    # values are distributed non-uniformly: special bit patterns are generated
    # with high probability in order to allow the test code to exercise corner
    # cases.

    void arf_randtest_not_zero(arf_t res, flint_rand_t state, slong bits, slong mag_bits)
    # Identical to :func:`arf_randtest`, except that zero is never produced
    # as an output.

    void arf_randtest_special(arf_t res, flint_rand_t state, slong bits, slong mag_bits)
    # Identical to :func:`arf_randtest`, except that the output occasionally
    # is set to an infinity or NaN.

    void arf_urandom(arf_t res, flint_rand_t state, slong bits, arf_rnd_t rnd)
    # Sets *res* to a uniformly distributed random number in the interval
    # `[0, 1]`. The method uses rounding from integers to floats based on the
    # rounding mode *rnd*.

    void arf_debug(const arf_t x)
    # Prints information about the internal representation of *x*.

    void arf_print(const arf_t x)
    # Prints *x* as an integer mantissa and exponent.

    void arf_printd(const arf_t x, slong d)
    # Prints *x* as a decimal floating-point number, rounding to *d* digits.
    # Rounding is faithful (at most 1 ulp error).

    char * arf_get_str(const arf_t x, slong d)
    # Returns *x* as a decimal floating-point number, rounding to *d* digits.
    # Rounding is faithful (at most 1 ulp error).

    void arf_fprint(FILE * file, const arf_t x)
    # Prints *x* as an integer mantissa and exponent to the stream *file*.

    void arf_fprintd(FILE * file, const arf_t y, slong d)
    # Prints *x* as a decimal floating-point number to the stream *file*,
    # rounding to *d* digits.
    # Rounding is faithful (at most 1 ulp error).

    char * arf_dump_str(const arf_t x)
    # Allocates a string and writes a binary representation of *x* to it that can
    # be read by :func:`arf_load_str`. The returned string needs to be
    # deallocated with *flint_free*.

    int arf_load_str(arf_t x, const char * str)
    # Parses *str* into *x*. Returns a nonzero value if *str* is not formatted
    # correctly.

    int arf_dump_file(FILE * stream, const arf_t x)
    # Writes a binary representation of *x* to *stream* that can be read by
    # :func:`arf_load_file`. Returns a nonzero value if the data could not be
    # written.

    int arf_load_file(arf_t x, FILE * stream)
    # Reads *x* from *stream*. Returns a nonzero value if the data is not
    # formatted correctly or the read failed. Note that the data is assumed to be
    # delimited by a whitespace or end-of-file, i.e., when writing multiple
    # values with :func:`arf_dump_file` make sure to insert a whitespace to
    # separate consecutive values.

    void arf_abs(arf_t res, const arf_t x)
    # Sets *res* to the absolute value of *x* exactly.

    void arf_neg(arf_t res, const arf_t x)
    # Sets *res* to `-x` exactly.

    int arf_neg_round(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
    # Sets *res* to `-x`.

    int arf_add(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_add_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

    int arf_add_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

    int arf_add_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x + y`.

    int arf_add_fmpz_2exp(arf_t res, const arf_t x, const fmpz_t y, const fmpz_t e, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x + y 2^e`.

    int arf_sub(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_sub_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

    int arf_sub_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

    int arf_sub_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x - y`.

    void arf_mul_2exp_si(arf_t res, const arf_t x, slong e)

    void arf_mul_2exp_fmpz(arf_t res, const arf_t x, const fmpz_t e)
    # Sets *res* to `x 2^e` exactly.

    int arf_mul(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_mul_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

    int arf_mul_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

    int arf_mul_mpz(arf_t res, const arf_t x, const mpz_t y, slong prec, arf_rnd_t rnd)

    int arf_mul_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x \cdot y`.

    int arf_addmul(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_addmul_ui(arf_t z, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

    int arf_addmul_si(arf_t z, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

    int arf_addmul_mpz(arf_t z, const arf_t x, const mpz_t y, slong prec, arf_rnd_t rnd)

    int arf_addmul_fmpz(arf_t z, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
    # Performs a fused multiply-add `z = z + x \cdot y`, updating *z* in-place.

    int arf_submul(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_submul_ui(arf_t z, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

    int arf_submul_si(arf_t z, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

    int arf_submul_mpz(arf_t z, const arf_t x, const mpz_t y, slong prec, arf_rnd_t rnd)

    int arf_submul_fmpz(arf_t z, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
    # Performs a fused multiply-subtract `z = z - x \cdot y`, updating *z* in-place.

    int arf_fma(arf_t res, const arf_t x, const arf_t y, const arf_t z, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x \cdot y + z`. This is equivalent to an *addmul* except
    # that *res* and *z* can be separate variables.

    int arf_sosq(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x^2 + y^2`, rounded to *prec* bits in the direction specified by *rnd*.

    int arf_sum(arf_t res, arf_srcptr terms, slong len, slong prec, arf_rnd_t rnd)
    # Sets *res* to the sum of the array *terms* of length *len*, rounded to
    # *prec* bits in the direction specified by *rnd*. The sum is computed as if
    # done without any intermediate rounding error, with only a single rounding
    # applied to the final result. Unlike repeated calls to :func:`arf_add` with
    # infinite precision, this function does not overflow if the magnitudes of
    # the terms are far apart. Warning: this function is implemented naively,
    # and the running time is quadratic with respect to *len* in the worst case.

    void arf_approx_dot(arf_t res, const arf_t initial, int subtract, arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)
    # Computes an approximate dot product, with the same meaning of
    # the parameters as :func:`arb_dot`.
    # This operation is not correctly rounded: the final rounding is done
    # in the direction ``rnd`` but intermediate roundings are
    # implementation-defined.

    int arf_div(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_div_ui(arf_t res, const arf_t x, ulong y, slong prec, arf_rnd_t rnd)

    int arf_ui_div(arf_t res, ulong x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_div_si(arf_t res, const arf_t x, slong y, slong prec, arf_rnd_t rnd)

    int arf_si_div(arf_t res, slong x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_div_fmpz(arf_t res, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)

    int arf_fmpz_div(arf_t res, const fmpz_t x, const arf_t y, slong prec, arf_rnd_t rnd)

    int arf_fmpz_div_fmpz(arf_t res, const fmpz_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x / y`, rounded to *prec* bits in the direction specified by *rnd*,
    # returning nonzero iff the operation is inexact. The result is NaN if *y* is zero.

    int arf_sqrt(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)

    int arf_sqrt_ui(arf_t res, ulong x, slong prec, arf_rnd_t rnd)

    int arf_sqrt_fmpz(arf_t res, const fmpz_t x, slong prec, arf_rnd_t rnd)
    # Sets *res* to `\sqrt{x}`. The result is NaN if *x* is negative.

    int arf_rsqrt(arf_t res, const arf_t x, slong prec, arf_rnd_t rnd)
    # Sets *res* to `1/\sqrt{x}`. The result is NaN if *x* is
    # negative, and `+\infty` if *x* is zero.

    int arf_root(arf_t res, const arf_t x, ulong k, slong prec, arf_rnd_t rnd)
    # Sets *res* to `x^{1/k}`. The result is NaN if *x* is negative.
    # Warning: this function is a wrapper around the MPFR root function.
    # It gets slow and uses much memory for large *k*.
    # Consider working with :func:`arb_root_ui` for large *k* instead of using this
    # function directly.

    int arf_complex_mul(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)

    int arf_complex_mul_fallback(arf_t e, arf_t f, const arf_t a, const arf_t b, const arf_t c, const arf_t d, slong prec, arf_rnd_t rnd)
    # Computes the complex product `e + fi = (a + bi)(c + di)`, rounding both
    # `e` and `f` correctly to *prec* bits in the direction specified by *rnd*.
    # The first bit in the return code indicates inexactness of `e`, and the
    # second bit indicates inexactness of `f`.
    # If any of the components *a*, *b*, *c*, *d* is zero, two real
    # multiplications and no additions are done. This convention is used even
    # if any other part contains an infinity or NaN, and the behavior
    # with infinite/NaN input is defined accordingly.
    # The *fallback* version is implemented naively, for testing purposes.
    # No squaring optimization is implemented.

    int arf_complex_sqr(arf_t e, arf_t f, const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)
    # Computes the complex square `e + fi = (a + bi)^2`. This function has
    # identical semantics to :func:`arf_complex_mul` (with `c = a, b = d`),
    # but is faster.

    int _arf_get_integer_mpn(mp_ptr y, mp_srcptr xp, mp_size_t xn, slong exp)
    # Given a floating-point number *x* represented by *xn* limbs at *xp*
    # and an exponent *exp*, writes the integer part of *x* to
    # *y*, returning whether the result is inexact.
    # The correct number of limbs is written (no limbs are written
    # if the integer part of *x* is zero).
    # Assumes that ``xp[0]`` is nonzero and that the
    # top bit of ``xp[xn-1]`` is set.

    int _arf_set_mpn_fixed(arf_t z, mp_srcptr xp, mp_size_t xn, mp_size_t fixn, int negative, slong prec, arf_rnd_t rnd)
    # Sets *z* to the fixed-point number having *xn* total limbs and *fixn*
    # fractional limbs, negated if *negative* is set, rounding *z* to *prec*
    # bits in the direction *rnd* and returning whether the result is inexact.
    # Both *xn* and *fixn* must be nonnegative and not so large
    # that the bit shift would overflow an *slong*, but otherwise no
    # assumptions are made about the input.

    int _arf_set_round_ui(arf_t z, ulong x, int sgnbit, slong prec, arf_rnd_t rnd)
    # Sets *z* to the integer *x*, negated if *sgnbit* is 1, rounded to *prec*
    # bits in the direction specified by *rnd*. There are no assumptions on *x*.

    int _arf_set_round_uiui(arf_t z, slong * fix, mp_limb_t hi, mp_limb_t lo, int sgnbit, slong prec, arf_rnd_t rnd)
    # Sets the mantissa of *z* to the two-limb mantissa given by *hi* and *lo*,
    # negated if *sgnbit* is 1, rounded to *prec* bits in the direction specified
    # by *rnd*. Requires that not both *hi* and *lo* are zero.
    # Writes the exponent shift to *fix* without writing the exponent of *z*
    # directly.

    int _arf_set_round_mpn(arf_t z, slong * exp_shift, mp_srcptr x, mp_size_t xn, int sgnbit, slong prec, arf_rnd_t rnd)
    # Sets the mantissa of *z* to the mantissa given by the *xn* limbs in *x*,
    # negated if *sgnbit* is 1, rounded to *prec* bits in the direction
    # specified by *rnd*. Returns the inexact flag. Requires that *xn* is positive
    # and that the top limb of *x* is nonzero. If *x* has leading zero bits,
    # writes the shift to *exp_shift*. This method does not write the exponent of
    # *z* directly. Requires that *x* does not point to the limbs of *z*.
