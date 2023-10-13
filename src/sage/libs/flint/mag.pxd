# distutils: libraries = flint
# distutils: depends = flint/mag.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void mag_init(mag_t x)
    # Initializes the variable *x* for use. Its value is set to zero.

    void mag_clear(mag_t x)
    # Clears the variable *x*, freeing or recycling its allocated memory.

    void mag_swap(mag_t x, mag_t y)
    # Swaps *x* and *y* efficiently.

    mag_ptr _mag_vec_init(long n)
    # Allocates a vector of length *n*. All entries are set to zero.

    void _mag_vec_clear(mag_ptr v, long n)
    # Clears a vector of length *n*.

    long mag_allocated_bytes(const mag_t x)
    # Returns the total number of bytes heap-allocated internally by this object.
    # The count excludes the size of the structure itself. Add
    # ``sizeof(mag_struct)`` to get the size of the object as a whole.

    void mag_zero(mag_t res)
    # Sets *res* to zero.

    void mag_one(mag_t res)
    # Sets *res* to one.

    void mag_inf(mag_t res)
    # Sets *res* to positive infinity.

    int mag_is_special(const mag_t x)
    # Returns nonzero iff *x* is zero or positive infinity.

    int mag_is_zero(const mag_t x)
    # Returns nonzero iff *x* is zero.

    int mag_is_inf(const mag_t x)
    # Returns nonzero iff *x* is positive infinity.

    int mag_is_finite(const mag_t x)
    # Returns nonzero iff *x* is not positive infinity (since there is no
    # NaN value, this function is exactly the logical negation of :func:`mag_is_inf`).

    void mag_init_set(mag_t res, const mag_t x)
    # Initializes *res* and sets it to the value of *x*. This operation is always exact.

    void mag_set(mag_t res, const mag_t x)
    # Sets *res* to the value of *x*. This operation is always exact.

    void mag_set_d(mag_t res, double x)

    void mag_set_ui(mag_t res, unsigned long x)

    void mag_set_fmpz(mag_t res, const fmpz_t x)
    # Sets *res* to an upper bound for `|x|`. The operation may be inexact
    # even if *x* is exactly representable.

    void mag_set_d_lower(mag_t res, double x)

    void mag_set_ui_lower(mag_t res, unsigned long x)

    void mag_set_fmpz_lower(mag_t res, const fmpz_t x)
    # Sets *res* to a lower bound for `|x|`.
    # The operation may be inexact even if *x* is exactly representable.

    void mag_set_d_2exp_fmpz(mag_t res, double x, const fmpz_t y)

    void mag_set_fmpz_2exp_fmpz(mag_t res, const fmpz_t x, const fmpz_t y)

    void mag_set_ui_2exp_si(mag_t res, unsigned long x, long y)
    # Sets *res* to an upper bound for `|x| \cdot 2^y`.

    void mag_set_d_2exp_fmpz_lower(mag_t res, double x, const fmpz_t y)

    void mag_set_fmpz_2exp_fmpz_lower(mag_t res, const fmpz_t x, const fmpz_t y)
    # Sets *res* to a lower bound for `|x| \cdot 2^y`.

    double mag_get_d(const mag_t x)
    # Returns a *double* giving an upper bound for *x*.

    double mag_get_d_log2_approx(const mag_t x)
    # Returns a *double* approximating `\log_2(x)`, suitable for estimating
    # magnitudes (warning: not a rigorous bound).
    # The value is clamped between *COEFF_MIN* and *COEFF_MAX*.

    void mag_get_fmpq(fmpq_t res, const mag_t x)

    void mag_get_fmpz(fmpz_t res, const mag_t x)

    void mag_get_fmpz_lower(fmpz_t res, const mag_t x)
    # Sets *res*, respectively, to the exact rational number represented by *x*,
    # the integer exactly representing the ceiling function of *x*, or the
    # integer exactly representing the floor function of *x*.
    # These functions are unsafe: the user must check in advance that *x* is of
    # reasonable magnitude. If *x* is infinite or has a bignum exponent, an
    # abort will be raised. If the exponent otherwise is too large or too small,
    # the available memory could be exhausted resulting in undefined behavior.

    int mag_equal(const mag_t x, const mag_t y)
    # Returns nonzero iff *x* and *y* have the same value.

    int mag_cmp(const mag_t x, const mag_t y)
    # Returns negative, zero, or positive, depending on whether *x*
    # is smaller, equal, or larger than *y*.

    int mag_cmp_2exp_si(const mag_t x, long y)
    # Returns negative, zero, or positive, depending on whether *x*
    # is smaller, equal, or larger than `2^y`.

    void mag_min(mag_t res, const mag_t x, const mag_t y)

    void mag_max(mag_t res, const mag_t x, const mag_t y)
    # Sets *res* respectively to the smaller or the larger of *x* and *y*.

    void mag_print(const mag_t x)
    # Prints *x* to standard output.

    void mag_fprint(FILE * file, const mag_t x)
    # Prints *x* to the stream *file*.

    char * mag_dump_str(const mag_t x)
    # Allocates a string and writes a binary representation of *x* to it that can
    # be read by :func:`mag_load_str`. The returned string needs to be
    # deallocated with *flint_free*.

    int mag_load_str(mag_t x, const char * str)
    # Parses *str* into *x*. Returns a nonzero value if *str* is not formatted
    # correctly.

    int mag_dump_file(FILE * stream, const mag_t x)
    # Writes a binary representation of *x* to *stream* that can be read by
    # :func:`mag_load_file`. Returns a nonzero value if the data could not be
    # written.

    int mag_load_file(mag_t x, FILE * stream)
    # Reads *x* from *stream*. Returns a nonzero value if the data is not
    # formatted correctly or the read failed. Note that the data is assumed to be
    # delimited by a whitespace or end-of-file, i.e., when writing multiple
    # values with :func:`mag_dump_file` make sure to insert a whitespace to
    # separate consecutive values.

    void mag_randtest(mag_t res, flint_rand_t state, long expbits)
    # Sets *res* to a random finite value, with an exponent up to *expbits* bits large.

    void mag_randtest_special(mag_t res, flint_rand_t state, long expbits)
    # Like :func:`mag_randtest`, but also sometimes sets *res* to infinity.

    void mag_add(mag_t res, const mag_t x, const mag_t y)

    void mag_add_ui(mag_t res, const mag_t x, unsigned long y)
    # Sets *res* to an upper bound for `x + y`.

    void mag_add_lower(mag_t res, const mag_t x, const mag_t y)

    void mag_add_ui_lower(mag_t res, const mag_t x, unsigned long y)
    # Sets *res* to a lower bound for `x + y`.

    void mag_add_2exp_fmpz(mag_t res, const mag_t x, const fmpz_t e)
    # Sets *res* to an upper bound for `x + 2^e`.

    void mag_add_ui_2exp_si(mag_t res, const mag_t x, unsigned long y, long e)
    # Sets *res* to an upper bound for `x + y 2^e`.

    void mag_sub(mag_t res, const mag_t x, const mag_t y)
    # Sets *res* to an upper bound for `\max(x-y, 0)`.

    void mag_sub_lower(mag_t res, const mag_t x, const mag_t y)
    # Sets *res* to a lower bound for `\max(x-y, 0)`.

    void mag_mul_2exp_si(mag_t res, const mag_t x, long y)

    void mag_mul_2exp_fmpz(mag_t res, const mag_t x, const fmpz_t y)
    # Sets *res* to `x \cdot 2^y`. This operation is exact.

    void mag_mul(mag_t res, const mag_t x, const mag_t y)

    void mag_mul_ui(mag_t res, const mag_t x, unsigned long y)

    void mag_mul_fmpz(mag_t res, const mag_t x, const fmpz_t y)
    # Sets *res* to an upper bound for `xy`.

    void mag_mul_lower(mag_t res, const mag_t x, const mag_t y)

    void mag_mul_ui_lower(mag_t res, const mag_t x, unsigned long y)

    void mag_mul_fmpz_lower(mag_t res, const mag_t x, const fmpz_t y)
    # Sets *res* to a lower bound for `xy`.

    void mag_addmul(mag_t z, const mag_t x, const mag_t y)
    # Sets *z* to an upper bound for `z + xy`.

    void mag_div(mag_t res, const mag_t x, const mag_t y)

    void mag_div_ui(mag_t res, const mag_t x, unsigned long y)

    void mag_div_fmpz(mag_t res, const mag_t x, const fmpz_t y)
    # Sets *res* to an upper bound for `x / y`.

    void mag_div_lower(mag_t res, const mag_t x, const mag_t y)
    # Sets *res* to a lower bound for `x / y`.

    void mag_inv(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `1 / x`.

    void mag_inv_lower(mag_t res, const mag_t x)
    # Sets *res* to a lower bound for `1 / x`.

    void mag_fast_init_set(mag_t x, const mag_t y)
    # Initialises *x* and sets it to the value of *y*.

    void mag_fast_zero(mag_t res)
    # Sets *res* to zero.

    int mag_fast_is_zero(const mag_t x)
    # Returns nonzero iff *x* to zero.

    void mag_fast_mul(mag_t res, const mag_t x, const mag_t y)
    # Sets *res* to an upper bound for `xy`.

    void mag_fast_addmul(mag_t z, const mag_t x, const mag_t y)
    # Sets *z* to an upper bound for `z + xy`.

    void mag_fast_add_2exp_si(mag_t res, const mag_t x, long e)
    # Sets *res* to an upper bound for `x + 2^e`.

    void mag_fast_mul_2exp_si(mag_t res, const mag_t x, long e)
    # Sets *res* to an upper bound for `x 2^e`.

    void mag_pow_ui(mag_t res, const mag_t x, unsigned long e)

    void mag_pow_fmpz(mag_t res, const mag_t x, const fmpz_t e)
    # Sets *res* to an upper bound for `x^e`.

    void mag_pow_ui_lower(mag_t res, const mag_t x, unsigned long e)

    void mag_pow_fmpz_lower(mag_t res, const mag_t x, const fmpz_t e)
    # Sets *res* to a lower bound for `x^e`.

    void mag_sqrt(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `\sqrt{x}`.

    void mag_sqrt_lower(mag_t res, const mag_t x)
    # Sets *res* to a lower bound for `\sqrt{x}`.

    void mag_rsqrt(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `1/\sqrt{x}`.

    void mag_rsqrt_lower(mag_t res, const mag_t x)
    # Sets *res* to an lower bound for `1/\sqrt{x}`.

    void mag_hypot(mag_t res, const mag_t x, const mag_t y)
    # Sets *res* to an upper bound for `\sqrt{x^2 + y^2}`.

    void mag_root(mag_t res, const mag_t x, unsigned long n)
    # Sets *res* to an upper bound for `x^{1/n}`.

    void mag_log(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `\log(\max(1,x))`.

    void mag_log_lower(mag_t res, const mag_t x)
    # Sets *res* to a lower bound for `\log(\max(1,x))`.

    void mag_neg_log(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `-\log(\min(1,x))`, i.e. an upper
    # bound for `|\log(x)|` for `x \le 1`.

    void mag_neg_log_lower(mag_t res, const mag_t x)
    # Sets *res* to a lower bound for `-\log(\min(1,x))`, i.e. a lower
    # bound for `|\log(x)|` for `x \le 1`.

    void mag_log_ui(mag_t res, unsigned long n)
    # Sets *res* to an upper bound for `\log(n)`.

    void mag_log1p(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `\log(1+x)`. The bound is computed
    # accurately for small *x*.

    void mag_exp(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `\exp(x)`.

    void mag_exp_lower(mag_t res, const mag_t x)
    # Sets *res* to a lower bound for `\exp(x)`.

    void mag_expinv(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `\exp(-x)`.

    void mag_expinv_lower(mag_t res, const mag_t x)
    # Sets *res* to a lower bound for `\exp(-x)`.

    void mag_expm1(mag_t res, const mag_t x)
    # Sets *res* to an upper bound for `\exp(x) - 1`. The bound is computed
    # accurately for small *x*.

    void mag_exp_tail(mag_t res, const mag_t x, unsigned long N)
    # Sets *res* to an upper bound for `\sum_{k=N}^{\infty} x^k / k!`.

    void mag_binpow_uiui(mag_t res, unsigned long m, unsigned long n)
    # Sets *res* to an upper bound for `(1 + 1/m)^n`.

    void mag_geom_series(mag_t res, const mag_t x, unsigned long N)
    # Sets *res* to an upper bound for `\sum_{k=N}^{\infty} x^k`.

    void mag_const_pi(mag_t res)

    void mag_const_pi_lower(mag_t res)
    # Sets *res* to an upper (respectively lower) bound for `\pi`.

    void mag_atan(mag_t res, const mag_t x)

    void mag_atan_lower(mag_t res, const mag_t x)
    # Sets *res* to an upper (respectively lower) bound for `\operatorname{atan}(x)`.

    void mag_cosh(mag_t res, const mag_t x)

    void mag_cosh_lower(mag_t res, const mag_t x)

    void mag_sinh(mag_t res, const mag_t x)

    void mag_sinh_lower(mag_t res, const mag_t x)
    # Sets *res* to an upper or lower bound for `\cosh(x)` or `\sinh(x)`.

    void mag_fac_ui(mag_t res, unsigned long n)
    # Sets *res* to an upper bound for `n!`.

    void mag_rfac_ui(mag_t res, unsigned long n)
    # Sets *res* to an upper bound for `1/n!`.

    void mag_bin_uiui(mag_t res, unsigned long n, unsigned long k)
    # Sets *res* to an upper bound for the binomial coefficient `{n \choose k}`.

    void mag_bernoulli_div_fac_ui(mag_t res, unsigned long n)
    # Sets *res* to an upper bound for `|B_n| / n!` where `B_n` denotes
    # a Bernoulli number.

    void mag_polylog_tail(mag_t res, const mag_t z, long s, unsigned long d, unsigned long N)
    # Sets *res* to an upper bound for
    # .. math ::
    # \sum_{k=N}^{\infty} \frac{z^k \log^d(k)}{k^s}.
    # The bounding strategy is described in :ref:`algorithms_polylogarithms`.
    # Note: in applications where `s` in this formula may be
    # real or complex, the user can simply
    # substitute any convenient integer `s'` such that `s' \le \operatorname{Re}(s)`.

    void mag_hurwitz_zeta_uiui(mag_t res, unsigned long s, unsigned long a)
    # Sets *res* to an upper bound for `\zeta(s,a) = \sum_{k=0}^{\infty} (k+a)^{-s}`.
    # We use the formula
    # .. math ::
    # \zeta(s,a) \le \frac{1}{a^s} + \frac{1}{(s-1) a^{s-1}}
    # which is obtained by estimating the sum by an integral.
    # If `s \le 1` or `a = 0`, the bound is infinite.

from .mag_macros cimport *
