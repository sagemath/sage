# distutils: libraries = flint
# distutils: depends = flint/fmpz.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    fmpz PTR_TO_COEFF(__mpz_struct * ptr)

    __mpz_struct * COEFF_TO_PTR(fmpz f)

    int COEFF_IS_MPZ(fmpz f)

    __mpz_struct * _fmpz_new_mpz()

    void _fmpz_clear_mpz(fmpz f)

    void _fmpz_cleanup_mpz_content()

    void _fmpz_cleanup()

    __mpz_struct * _fmpz_promote(fmpz_t f)

    __mpz_struct * _fmpz_promote_val(fmpz_t f)

    void _fmpz_demote(fmpz_t f)

    void _fmpz_demote_val(fmpz_t f)

    void fmpz_init(fmpz_t f)
    # A small ``fmpz_t`` is initialised, i.e. just a ``slong``.
    # The value is set to zero.

    void fmpz_init2(fmpz_t f, unsigned long limbs)
    # Initialises the given ``fmpz_t`` to have space for the given
    # number of limbs.
    # If ``limbs`` is zero then a small ``fmpz_t`` is allocated,
    # i.e. just a ``slong``.  The value is also set to zero.  It is
    # not necessary to call this function except to save time.  A call
    # to ``fmpz_init`` will do just fine.

    void fmpz_clear(fmpz_t f)
    # Clears the given ``fmpz_t``, releasing any memory associated
    # with it, either back to the stack or the OS, depending on
    # whether the reentrant or non-reentrant version of FLINT is built.

    void fmpz_init_set(fmpz_t f, const fmpz_t g)

    void fmpz_init_set_ui(fmpz_t f, unsigned long g)

    void fmpz_init_set_si(fmpz_t f, long g)
    # Initialises `f` and sets it to the value of `g`.

    void fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    # Generates a random signed integer whose absolute value has precisely
    # the given number of bits.

    void fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    # Generates a random signed integer whose absolute value has a number
    # of bits which is random from `0` up to ``bits`` inclusive.

    void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    # Generates a random unsigned integer whose value has a number
    # of bits which is random from `0` up to ``bits`` inclusive.

    void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    # As per ``fmpz_randtest``, but the result will not be `0`.
    # If ``bits`` is set to `0`, an exception will result.

    void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m)
    # Generates a random integer in the range `0` to `m - 1` inclusive.

    void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m)
    # Generates a random integer in the range `0` to `m - 1` inclusive,
    # with an increased probability of generating values close to
    # the endpoints.

    void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m)
    # Generates a random integer in the range `(-m/2, m/2]`, with an
    # increased probability of generating values close to the
    # endpoints or close to zero.

    void fmpz_randprime(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits, int proved)
    # Generates a random prime number with the given number of bits.
    # The generation is performed by choosing a random number and then
    # finding the next largest prime, and therefore does not quite
    # give a uniform distribution over the set of primes with that
    # many bits.
    # Random number generation is performed using the standard Flint
    # random number generator, which is not suitable for cryptographic use.
    # If ``proved`` is nonzero, then the integer returned is
    # guaranteed to actually be prime.

    long fmpz_get_si(const fmpz_t f)
    # Returns `f` as a ``slong``.  The result is undefined
    # if `f` does not fit into a ``slong``.

    unsigned long fmpz_get_ui(const fmpz_t f)
    # Returns `f` as an ``ulong``.  The result is undefined
    # if `f` does not fit into an ``ulong`` or is negative.

    void fmpz_get_uiui(mp_limb_t * hi, mp_limb_t * low, const fmpz_t f)
    # If `f` consists of two limbs, then ``*hi`` and ``*low`` are set to the high
    # and low limbs, otherwise ``*low`` is set to the low limb and ``*hi`` is set
    # to `0`.

    mp_limb_t fmpz_get_nmod(const fmpz_t f, nmod_t mod)
    # Returns `f \mod n`.

    double fmpz_get_d(const fmpz_t f)
    # Returns `f` as a ``double``, rounding down towards zero if
    # `f` cannot be represented exactly. The outcome is undefined
    # if `f` is too large to fit in the normal range of a double.

    void fmpz_set_mpf(fmpz_t f, const mpf_t x)
    # Sets `f` to the ``mpf_t`` `x`, rounding down towards zero if
    # the value of `x` is fractional.

    void fmpz_get_mpf(mpf_t x, const fmpz_t f)
    # Sets the value of the ``mpf_t`` `x` to the value of `f`.

    void fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd)
    # Sets the value of `x` from `f`, rounded toward the given
    # direction ``rnd``.
    # **Note:** Requires that ``mpfr.h`` has been included before any FLINT
    # header is included.

    double fmpz_get_d_2exp(long * exp, const fmpz_t f)
    # Returns `f` as a normalized ``double`` along with a `2`-exponent
    # ``exp``, i.e. if `r` is the return value then `f = r 2^{exp}`,
    # to within 1 ULP.

    void fmpz_get_mpz(mpz_t x, const fmpz_t f)
    # Sets the ``mpz_t`` `x` to the same value as `f`.

    int fmpz_get_mpn(mp_ptr *n, fmpz_t n_in)
    # Sets the ``mp_ptr`` `n` to the same value as `n_{in}`. Returned
    # integer is number of limbs allocated to `n`, minimum number of limbs
    # required to hold the value stored in `n_{in}`.

    char * fmpz_get_str(char * str, int b, const fmpz_t f)
    # Returns the representation of `f` in base `b`, which can vary
    # between `2` and `62`, inclusive.
    # If ``str`` is ``NULL``, the result string is allocated by
    # the function.  Otherwise, it is up to the caller to ensure that
    # the allocated block of memory is sufficiently large.

    void fmpz_set_si(fmpz_t f, long val)
    # Sets `f` to the given ``slong`` value.

    void fmpz_set_ui(fmpz_t f, unsigned long val)
    # Sets `f` to the given ``ulong`` value.

    void fmpz_set_d(fmpz_t f, double c)
    # Sets `f` to the ``double`` `c`, rounding down towards zero if
    # the value of `c` is fractional. The outcome is undefined if `c` is
    # infinite, not-a-number, or subnormal.

    void fmpz_set_d_2exp(fmpz_t f, double d, long exp)
    # Sets `f` to the nearest integer to `d 2^{exp}`.

    void fmpz_neg_ui(fmpz_t f, unsigned long val)
    # Sets `f` to the given ``ulong`` value, and then negates `f`.

    void fmpz_set_uiui(fmpz_t f, mp_limb_t hi, mp_limb_t lo)
    # Sets `f` to ``lo``, plus ``hi`` shifted to the left by
    # ``FLINT_BITS``.

    void fmpz_neg_uiui(fmpz_t f, mp_limb_t hi, mp_limb_t lo)
    # Sets `f` to ``lo``, plus ``hi`` shifted to the left by
    # ``FLINT_BITS``, and then negates `f`.

    void fmpz_set_signed_uiui(fmpz_t f, unsigned long hi, unsigned long lo)
    # Sets `f` to ``lo``, plus ``hi`` shifted to the left by
    # ``FLINT_BITS``, interpreted as a signed two's complement
    # integer with ``2 * FLINT_BITS`` bits.

    void fmpz_set_signed_uiuiui(fmpz_t f, unsigned long hi, unsigned long mid, unsigned long lo)
    # Sets `f` to ``lo``, plus ``mid`` shifted to the left by
    # ``FLINT_BITS``, plus ``hi`` shifted to the left by
    # ``2*FLINT_BITS`` bits, interpreted as a signed two's complement
    # integer with ``3 * FLINT_BITS`` bits.

    void fmpz_set_ui_array(fmpz_t out, const unsigned long * input, long n)
    # Sets ``out`` to the nonnegative integer
    # ``in[0] + in[1]*X  + ... + in[n - 1]*X^(n - 1)``
    # where ``X = 2^FLINT_BITS``. It is assumed that ``n > 0``.

    void fmpz_set_signed_ui_array(fmpz_t out, const unsigned long * input, long n)
    # Sets ``out`` to the integer represented in ``in[0], ..., in[n - 1]``
    # as a signed two's complement integer with ``n * FLINT_BITS`` bits.
    # It is assumed that ``n > 0``. The function operates as a call to
    # :func:`fmpz_set_ui_array` followed by a symmetric remainder modulo
    # `2^{n\cdot FLINT\_BITS}`.

    void fmpz_get_ui_array(unsigned long * out, long n, const fmpz_t input)
    # Assuming that the nonnegative integer ``in`` can be represented in the
    # form ``out[0] + out[1]*X + ... + out[n - 1]*X^(n - 1)``,
    # where `X = 2^{FLINT\_BITS}`, sets the corresponding elements of ``out``
    # so that this is true. It is assumed that ``n > 0``.

    void fmpz_get_signed_ui_array(unsigned long * out, long n, const fmpz_t input)
    # Retrieves the value of `in` modulo `2^{n * FLINT\_BITS}` and puts the `n`
    # words of the result in ``out[0], ..., out[n-1]``. This will give a signed
    # two's complement representation of `in` (assuming `in` doesn't overflow the array).

    void fmpz_get_signed_uiui(unsigned long * hi, unsigned long * lo, const fmpz_t input)
    # Retrieves the value of `in` modulo `2^{2 * FLINT\_BITS}` and puts the high
    # and low words into ``*hi`` and ``*lo`` respectively.

    void fmpz_set_mpz(fmpz_t f, const mpz_t x)
    # Sets `f` to the given ``mpz_t`` value.

    int fmpz_set_str(fmpz_t f, const char * str, int b)
    # Sets `f` to the value given in the null-terminated string ``str``,
    # in base `b`. The base `b` can vary between `2` and `62`, inclusive.
    # Returns `0` if the string contains a valid input and `-1` otherwise.

    void fmpz_set_ui_smod(fmpz_t f, mp_limb_t x, mp_limb_t m)
    # Sets `f` to the signed remainder `y \equiv x \bmod m` satisfying
    # `-m/2 < y \leq m/2`, given `x` which is assumed to satisfy
    # `0 \leq x < m`.

    void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f)
    # Sets the uninitialised ``mpz_t`` `z` to the value of the
    # readonly ``fmpz_t`` `f`.
    # Note that it is assumed that `f` does not change during
    # the lifetime of `z`.
    # The integer `z` has to be cleared by a call to
    # :func:`flint_mpz_clear_readonly`.
    # The suggested use of the two functions is as follows::
    # fmpz_t f;
    # ...
    # {
    # mpz_t z;
    # flint_mpz_init_set_readonly(z, f);
    # foo(..., z);
    # flint_mpz_clear_readonly(z);
    # }
    # This provides a convenient function for user code, only
    # requiring to work with the types ``fmpz_t`` and ``mpz_t``.
    # In critical code, the following approach may be favourable::
    # fmpz_t f;
    # ...
    # {
    # __mpz_struct *z;
    # z = _fmpz_promote_val(f);
    # foo(..., z);
    # _fmpz_demote_val(f);
    # }

    void flint_mpz_clear_readonly(mpz_t z)
    # Clears the readonly ``mpz_t`` `z`.

    void fmpz_init_set_readonly(fmpz_t f, const mpz_t z)
    # Sets the uninitialised ``fmpz_t`` `f` to a readonly
    # version of the integer `z`.
    # Note that the value of `z` is assumed to remain constant
    # throughout the lifetime of `f`.
    # The ``fmpz_t`` `f` has to be cleared by calling the
    # function :func:`fmpz_clear_readonly`.
    # The suggested use of the two functions is as follows::
    # mpz_t z;
    # ...
    # {
    # fmpz_t f;
    # fmpz_init_set_readonly(f, z);
    # foo(..., f);
    # fmpz_clear_readonly(f);
    # }

    void fmpz_clear_readonly(fmpz_t f)
    # Clears the readonly ``fmpz_t`` `f`.

    int fmpz_read(fmpz_t f)
    # Reads a multiprecision integer from ``stdin``.  The format is
    # an optional minus sign, followed by one or more digits.  The
    # first digit should be non-zero unless it is the only digit.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive number.
    # This convention is adopted in light of the return values of
    # ``scanf`` from the standard library and ``mpz_inp_str``
    # from MPIR.

    int fmpz_fread(FILE * file, fmpz_t f)
    # Reads a multiprecision integer from the stream ``file``.  The
    # format is an optional minus sign, followed by one or more digits.
    # The first digit should be non-zero unless it is the only digit.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive number.
    # This convention is adopted in light of the return values of
    # ``scanf`` from the standard library and ``mpz_inp_str``
    # from MPIR.

    size_t fmpz_inp_raw(fmpz_t x, FILE *fin )
    # Reads a multiprecision integer from the stream ``file``.  The
    # format is raw binary format write by :func:`fmpz_out_raw`.
    # In case of success, return a positive number, indicating number of bytes read.
    # In case of failure 0.
    # This function calls the ``mpz_inp_raw`` function in library gmp. So that it
    # can read the raw data written by ``mpz_inp_raw`` directly.

    int fmpz_print(const fmpz_t x)
    # Prints the value `x` to ``stdout``, without a carriage return (CR).
    # The value is printed as either `0`, the decimal digits of a
    # positive integer, or a minus sign followed by the digits of
    # a negative integer.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive number.
    # This convention is adopted in light of the return values of
    # ``flint_printf`` from the standard library and ``mpz_out_str``
    # from MPIR.

    int fmpz_fprint(FILE * file, const fmpz_t x)
    # Prints the value `x` to ``file``, without a carriage return (CR).
    # The value is printed as either `0`, the decimal digits of a
    # positive integer, or a minus sign followed by the digits of
    # a negative integer.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive number.
    # This convention is adopted in light of the return values of
    # ``flint_printf`` from the standard library and ``mpz_out_str``
    # from MPIR.

    size_t fmpz_out_raw(FILE *fout, const fmpz_t x )
    # Writes the value `x` to ``file``.
    # The value is written in raw binary format. The integer is written in
    # portable format, with 4 bytes of size information, and that many bytes
    # of limbs. Both the size and the limbs are written in decreasing
    # significance order (i.e., in big-endian).
    # The output can be read with ``fmpz_inp_raw``.
    # In case of success, return a positive number, indicating number of bytes written.
    # In case of failure, return 0.
    # The output of this can also be read by ``mpz_inp_raw`` from GMP >= 2,
    # since this function calls the ``mpz_inp_raw`` function in library gmp.

    size_t fmpz_sizeinbase(const fmpz_t f, int b)
    # Returns the size of the absolute value of `f` in base `b`, measured in
    # numbers of digits. The base `b` can be between `2` and `62`, inclusive.

    flint_bitcnt_t fmpz_bits(const fmpz_t f)
    # Returns the number of bits required to store the absolute
    # value of `f`.  If `f` is `0` then `0` is returned.

    mp_size_t fmpz_size(const fmpz_t f)
    # Returns the number of limbs required to store the absolute
    # value of `f`.  If `f` is zero then `0` is returned.

    int fmpz_sgn(const fmpz_t f)
    # Returns `-1` if the sign of `f` is negative, `+1` if it is positive,
    # otherwise returns `0`.

    flint_bitcnt_t fmpz_val2(const fmpz_t f)
    # Returns the exponent of the largest power of two dividing `f`, or
    # equivalently the number of trailing zeros in the binary expansion of `f`.
    # If `f` is zero then `0` is returned.

    void fmpz_swap(fmpz_t f, fmpz_t g)
    # Efficiently swaps `f` and `g`.  No data is copied.

    void fmpz_set(fmpz_t f, const fmpz_t g)
    # Sets `f` to the same value as `g`.

    void fmpz_zero(fmpz_t f)
    # Sets `f` to zero.

    void fmpz_one(fmpz_t f)
    # Sets `f` to one.

    int fmpz_abs_fits_ui(const fmpz_t f)
    # Returns whether the absolute value of `f`
    # fits into an ``ulong``.

    int fmpz_fits_si(const fmpz_t f)
    # Returns whether the value of `f` fits into a ``slong``.

    void fmpz_setbit(fmpz_t f, unsigned long i)
    # Sets bit index `i` of `f`.

    int fmpz_tstbit(const fmpz_t f, unsigned long i)
    # Test bit index `i` of `f` and return `0` or `1`, accordingly.

    mp_limb_t fmpz_abs_lbound_ui_2exp(long * exp, const fmpz_t x, int bits)
    # For nonzero `x`, returns a mantissa `m` with exactly ``bits`` bits and
    # sets ``exp`` to an exponent `e`, such that `|x| \ge m 2^e`. The number
    # of bits must be between 1 and ``FLINT_BITS`` inclusive.
    # The mantissa is guaranteed to be correctly rounded.

    mp_limb_t fmpz_abs_ubound_ui_2exp(long * exp, const fmpz_t x, int bits)
    # For nonzero `x`, returns a mantissa `m` with exactly ``bits`` bits
    # and sets ``exp`` to an exponent `e`, such that `|x| \le m 2^e`.
    # The number of bits must be between 1 and ``FLINT_BITS`` inclusive.
    # The mantissa is either correctly rounded or one unit too large
    # (possibly meaning that the exponent is one too large,
    # if the mantissa is a power of two).

    int fmpz_cmp(const fmpz_t f, const fmpz_t g)

    int fmpz_cmp_ui(const fmpz_t f, unsigned long g)

    int fmpz_cmp_si(const fmpz_t f, long g)
    # Returns a negative value if `f < g`, positive value if `g < f`,
    # otherwise returns `0`.

    int fmpz_cmpabs(const fmpz_t f, const fmpz_t g)
    # Returns a negative value if `\lvert f\rvert < \lvert g\rvert`, positive value if
    # `\lvert g\rvert < \lvert f \rvert`, otherwise returns `0`.

    int fmpz_cmp2abs(const fmpz_t f, const fmpz_t g)
    # Returns a negative value if `\lvert f\rvert < \lvert 2g\rvert`, positive value if
    # `\lvert 2g\rvert < \lvert f \rvert`, otherwise returns `0`.

    int fmpz_equal(const fmpz_t f, const fmpz_t g)

    int fmpz_equal_ui(const fmpz_t f, unsigned long g)

    int fmpz_equal_si(const fmpz_t f, long g)
    # Returns `1` if `f` is equal to `g`, otherwise returns `0`.

    int fmpz_is_zero(const fmpz_t f)
    # Returns `1` if `f` is `0`, otherwise returns `0`.

    int fmpz_is_one(const fmpz_t f)
    # Returns `1` if `f` is equal to one, otherwise returns `0`.

    int fmpz_is_pm1(const fmpz_t f)
    # Returns `1` if `f` is equal to one or minus one, otherwise returns `0`.

    int fmpz_is_even(const fmpz_t f)
    # Returns whether the integer `f` is even.

    int fmpz_is_odd(const fmpz_t f)
    # Returns whether the integer `f` is odd.

    void fmpz_neg(fmpz_t f1, const fmpz_t f2)
    # Sets `f_1` to `-f_2`.

    void fmpz_abs(fmpz_t f1, const fmpz_t f2)
    # Sets `f_1` to the absolute value of `f_2`.

    void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_add_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    void fmpz_add_si(fmpz_t f, const fmpz_t g, long h)
    # Sets `f` to `g + h`.

    void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_sub_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    void fmpz_sub_si(fmpz_t f, const fmpz_t g, long h)
    # Sets `f` to `g - h`.

    void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_mul_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    void fmpz_mul_si(fmpz_t f, const fmpz_t g, long h)
    # Sets `f` to `g \times h`.

    void fmpz_mul2_uiui(fmpz_t f, const fmpz_t g, unsigned long x, unsigned long y)
    # Sets `f` to `g \times x \times y` where `x` and `y` are of type ``ulong``.

    void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, unsigned long e)
    # Sets `f` to `g \times 2^e`.
    # Note: Assumes that ``e + FLINT_BITS`` does not overflow.

    void fmpz_one_2exp(fmpz_t f, unsigned long e)
    # Sets `f` to `2^e`.

    void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    void fmpz_addmul_si(fmpz_t f, const fmpz_t g, long h)
    # Sets `f` to `f + g \times h`.

    void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_submul_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    void fmpz_submul_si(fmpz_t f, const fmpz_t g, long h)
    # Sets `f` to `f - g \times h`.

    void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d)
    # Sets `f` to `a \times b + c \times d`.

    void fmpz_fmms(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d)
    # Sets `f` to `a \times b - c \times d`.

    void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)

    void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)

    void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)

    void fmpz_ndiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)

    void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)

    void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)

    void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)

    void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, long h)

    void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, long h)

    void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, long h)

    void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, unsigned long h)

    void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, unsigned long h)

    void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, unsigned long h)

    void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, unsigned long exp)

    void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, unsigned long exp)

    void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, unsigned long exp)

    void fmpz_fdiv_r(fmpz_t s, const fmpz_t g, const fmpz_t h)

    void fmpz_cdiv_r_2exp(fmpz_t s, const fmpz_t g, unsigned long exp)

    void fmpz_fdiv_r_2exp(fmpz_t s, const fmpz_t g, unsigned long exp)

    void fmpz_tdiv_r_2exp(fmpz_t s, const fmpz_t g, unsigned long exp)
    # Sets `f` to the quotient of `g` by `h` and/or `s` to the remainder. For the
    # ``2exp`` functions, ``g = 2^exp``. `If `h` is `0` an exception is raised.
    # Rounding is made in the following way:
    # * ``fdiv`` rounds the quotient via floor rounding.
    # * ``cdiv`` rounds the quotient via ceil rounding.
    # * ``tdiv`` rounds the quotient via truncation, i.e. rounding towards zero.
    # * ``ndiv`` rounds the quotient such that the remainder has the smallest
    # absolute value. In case of ties, it rounds the quotient towards zero.

    unsigned long fmpz_cdiv_ui(const fmpz_t g, unsigned long h)

    unsigned long fmpz_fdiv_ui(const fmpz_t g, unsigned long h)

    unsigned long fmpz_tdiv_ui(const fmpz_t g, unsigned long h)

    void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h)

    void fmpz_divexact_si(fmpz_t f, const fmpz_t g, long h)

    void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    # Sets `f` to the quotient of `g` and `h`, assuming that the
    # division is exact, i.e. `g` is a multiple of `h`.  If `h`
    # is `0` an exception is raised.

    void fmpz_divexact2_uiui(fmpz_t f, const fmpz_t g, unsigned long x, unsigned long y)
    # Sets `f` to the quotient of `g` and `h = x \times y`, assuming that
    # the division is exact, i.e. `g` is a multiple of `h`.
    # If `x` or `y` is `0` an exception is raised.

    int fmpz_divisible(const fmpz_t f, const fmpz_t g)

    int fmpz_divisible_si(const fmpz_t f, long g)
    # Returns `1` if there is an integer `q` with `f = q g` and `0` if there is
    # none.

    int fmpz_divides(fmpz_t q, const fmpz_t g, const fmpz_t h)
    # Returns `1` if there is an integer `q` with `f = q g` and sets `q` to the
    # quotient. Otherwise returns `0` and sets `q` to `0`.

    void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    # Sets `f` to the remainder of `g` divided by `h` such that the remainder is
    # positive. Assumes that `h` is not zero.

    unsigned long fmpz_mod_ui(fmpz_t f, const fmpz_t g, unsigned long h)
    # Sets `f` to the remainder of `g` divided by `h` such that the remainder is
    # positive and also returns this value. Raises an exception if `h` is zero.

    void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    # Sets `f` to the signed remainder `y \equiv g \bmod h` satisfying
    # `-\lvert h \rvert/2 < y \leq \lvert h\rvert/2`.

    void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f)
    # Compute a precomputed inverse ``inv`` of ``f`` for use in the
    # ``preinvn`` functions listed below.

    void fmpz_preinvn_clear(fmpz_preinvn_t inv)
    # Clean up the resources used by a precomputed inverse created with the
    # :func:`fmpz_preinvn_init` function.

    void fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h, const fmpz_preinvn_t hinv)
    # As per :func:`fmpz_fdiv_qr`, but takes a precomputed inverse ``hinv``
    # of `h` constructed using :func:`fmpz_preinvn`.
    # This function will be faster than :func:`fmpz_fdiv_qr_preinvn` when the
    # number of limbs of `h` is at least ``PREINVN_CUTOFF``.

    void fmpz_pow_ui(fmpz_t f, const fmpz_t g, unsigned long x)
    void fmpz_ui_pow_ui(fmpz_t f, unsigned long g, unsigned long x)
    # Sets `f` to `g^x`.  Defines `0^0 = 1`.

    int fmpz_pow_fmpz(fmpz_t f, const fmpz_t g, const fmpz_t x)
    # Sets `f` to `g^x`. Defines `0^0 = 1`. Return `1` for success and `0` for
    # failure. The function throws only if `x` is negative.

    void fmpz_powm_ui(fmpz_t f, const fmpz_t g, unsigned long e, const fmpz_t m)

    void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)
    # Sets `f` to `g^e \bmod{m}`.  If `e = 0`, sets `f` to `1`.
    # Assumes that `m \neq 0`, raises an ``abort`` signal otherwise.

    long fmpz_clog(const fmpz_t x, const fmpz_t b)
    long fmpz_clog_ui(const fmpz_t x, unsigned long b)
    # Returns `\lceil\log_b x\rceil`.
    # Assumes that `x \geq 1` and `b \geq 2` and that
    # the return value fits into a signed ``slong``.

    long fmpz_flog(const fmpz_t x, const fmpz_t b)
    long fmpz_flog_ui(const fmpz_t x, unsigned long b)
    # Returns `\lfloor\log_b x\rfloor`.
    # Assumes that `x \geq 1` and `b \geq 2` and that
    # the return value fits into a signed ``slong``.

    double fmpz_dlog(const fmpz_t x)
    # Returns a double precision approximation of the
    # natural logarithm of `x`.
    # The accuracy depends on the implementation of the floating-point
    # logarithm provided by the C standard library. The result can
    # typically be expected to have a relative error no greater than 1-2 bits.

    int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p)
    # If `p` is prime, set `b` to a square root of `a` modulo `p` if `a` is a
    # quadratic residue modulo `p` and return `1`, otherwise return `0`.
    # If `p` is not prime the return value is with high probability `0`,
    # indicating that `p` is not prime, or `a` is not a square modulo `p`.
    # If `p` is not prime and the return value is `1`, the value of `b` is
    # meaningless.

    void fmpz_sqrt(fmpz_t f, const fmpz_t g)
    # Sets `f` to the integer part of the square root of `g`, where
    # `g` is assumed to be non-negative.  If `g` is negative, an exception
    # is raised.

    void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g)
    # Sets `f` to the integer part of the square root of `g`, where `g` is
    # assumed to be non-negative, and sets `r` to the remainder, that is,
    # the difference `g - f^2`.  If `g` is negative, an exception is raised.
    # The behaviour is undefined if `f` and `r` are aliases.

    int fmpz_is_square(const fmpz_t f)
    # Returns nonzero if `f` is a perfect square and zero otherwise.

    int fmpz_root(fmpz_t r, const fmpz_t f, long n)
    # Set `r` to the integer part of the `n`-th root of `f`. Requires that
    # `n > 0` and that if `n` is even then `f` be non-negative, otherwise an
    # exception is raised. The function returns `1` if the root was exact,
    # otherwise `0`.

    int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f)
    # If `f` is a perfect power `r^k` set ``root`` to `r` and return `k`,
    # otherwise return `0`. Note that `-1, 0, 1` are all considered perfect
    # powers. No guarantee is made about `r` or `k` being the smallest
    # possible value. Negative values of `f` are permitted.

    void fmpz_fac_ui(fmpz_t f, unsigned long n)
    # Sets `f` to the factorial `n!` where `n` is an ``ulong``.

    void fmpz_fib_ui(fmpz_t f, unsigned long n)
    # Sets `f` to the Fibonacci number `F_n` where `n` is an
    # ``ulong``.

    void fmpz_bin_uiui(fmpz_t f, unsigned long n, unsigned long k)
    # Sets `f` to the binomial coefficient `{n \choose k}`.

    void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, unsigned long a, unsigned long b)
    # Sets `r` to the rising factorial `(x+a) (x+a+1) (x+a+2) \cdots (x+b-1)`.
    # Assumes `b > a`.

    void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, unsigned long k)
    # Sets `r` to the rising factorial `x (x+1) (x+2) \cdots (x+k-1)`.

    void fmpz_rfac_uiui(fmpz_t r, unsigned long x, unsigned long k)
    # Sets `r` to the rising factorial `x (x+1) (x+2) \cdots (x+k-1)`.

    void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, unsigned long exp)
    # Sets `f` to the product of `g` and `h` divided by ``2^exp``, rounding
    # down towards zero.

    void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, long x, unsigned long exp)
    # Sets `f` to the product of `g` and `x` divided by ``2^exp``, rounding
    # down towards zero.

    void fmpz_gcd_ui(fmpz_t f, const fmpz_t g, unsigned long h)

    void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
    # Sets `f` to the greatest common divisor of `g` and `h`.  The
    # result is always positive, even if one of `g` and `h` is
    # negative.

    void fmpz_gcd3(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c)
    # Sets `f` to the greatest common divisor of `a`, `b` and `c`.
    # This is equivalent to calling ``fmpz_gcd`` twice, but may be faster.

    void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h)
    # Sets `f` to the least common multiple of `g` and `h`.  The
    # result is always nonnegative, even if one of `g` and `h` is
    # negative.

    void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)
    # Given integers `f, g` with `0 \leq f < g`, computes the
    # greatest common divisor `d = \gcd(f, g)` and the modular
    # inverse `a = f^{-1} \pmod{g}`, whenever `f \neq 0`.
    # Assumes that `d` and `a` are not aliased.

    void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
    # Computes the extended GCD of `f` and `g`, i.e. the values `a` and `b` such
    # that `af + bg = d`, where `d = \gcd(f, g)`. Here `a` will be the same as
    # calling ``fmpz_gcdinv`` when `f < g` (or vice versa for `b` when `g < f`).
    # To obtain the canonical solution to Bézout's identity, call
    # ``fmpz_xgcd_canonical_bezout`` instead. This is also faster.
    # Assumes that there is no aliasing among the outputs.

    void fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
    # Computes the extended GCD `\operatorname{xgcd}(f, g) = (d, a, b)` such that
    # the solution is the canonical solution to Bézout's identity. We define the
    # canonical solution to satisfy one of the following if one of the given
    # conditions apply:
    # .. math ::
    # \operatorname{xgcd}(\pm g, g) &= \bigl(|g|, 0, \operatorname{sgn}(g)\bigr)
    # \operatorname{xgcd}(f, 0) &= \bigl(|f|, \operatorname{sgn}(f), 0\bigr)
    # \operatorname{xgcd}(0, g) &= \bigl(|g|, 0, \operatorname{sgn}(g)\bigr)
    # \operatorname{xgcd}(f, \mp 1) &= (1, 0, \mp 1)
    # \operatorname{xgcd}(\mp 1, g) &= (1, \mp 1, 0)\quad g \neq 0, \pm 1
    # \operatorname{xgcd}(\mp 2 d, g) &=
    # \bigl(d, {\textstyle\frac{d - |g|}{\mp 2 d}}, \operatorname{sgn}(g)\bigr)
    # \operatorname{xgcd}(f, \mp 2 d) &=
    # \bigl(d, \operatorname{sgn}(f), {\textstyle\frac{d - |g|}{\mp 2 d}}\bigr).
    # If the pair `(f, g)` does not satisfy any of these conditions, the solution
    # `(d, a, b)` will satisfy the following:
    # .. math ::
    # |a| < \Bigl| \frac{g}{2 d} \Bigr|,
    # \qquad |b| < \Bigl| \frac{f}{2 d} \Bigr|.
    # Assumes that there is no aliasing among the outputs.

    void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, fmpz_t r2, fmpz_t r1, const fmpz_t L)
    # This function is an implementation of Lehmer extended GCD with early
    # termination, as used in the ``qfb`` module. It terminates early when
    # remainders fall below the specified bound. The initial values ``r1``
    # and ``r2`` are treated as successive remainders in the Euclidean
    # algorithm and are replaced with the last two remainders computed. The
    # values ``co1`` and ``co2`` are the last two cofactors and satisfy
    # the identity ``co2*r1 - co1*r2 == +/- r2_orig`` upon termination, where
    # ``r2_orig`` is the starting value of ``r2`` supplied, and ``r1``
    # and ``r2`` are the final values.
    # Aliasing of inputs is not allowed. Similarly aliasing of inputs and outputs
    # is not allowed.

    long _fmpz_remove(fmpz_t x, const fmpz_t f, double finv)
    # Removes all factors `f` from `x` and returns the number of such.
    # Assumes that `x` is non-zero, that `f > 1` and that ``finv``
    # is the precomputed ``double`` inverse of `f` whenever `f` is
    # a small integer and `0` otherwise.
    # Does not support aliasing.

    long fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f)
    # Remove all occurrences of the factor `f > 1` from the
    # integer ``op`` and sets ``rop`` to the resulting
    # integer.
    # If ``op`` is zero, sets ``rop`` to ``op`` and
    # returns `0`.
    # Returns an ``abort`` signal if any of the assumptions
    # are violated.

    int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    # Sets `f` to the inverse of `g` modulo `h`.  The value of `h` may
    # not be `0` otherwise an exception results.  If the inverse exists
    # the return value will be non-zero, otherwise the return value will
    # be `0` and the value of `f` undefined. As a special case, we
    # consider any number invertible modulo `h = \pm 1`, with inverse 0.

    void fmpz_negmod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    # Sets `f` to `-g \pmod{h}`, assuming `g` is reduced modulo `h`.

    int fmpz_jacobi(const fmpz_t a, const fmpz_t n)
    # Computes the Jacobi symbol `\left(\frac{a}{n}\right)` for any `a` and odd positive `n`.

    int fmpz_kronecker(const fmpz_t a, const fmpz_t n)
    # Computes the Kronecker symbol `\left(\frac{a}{n}\right)` for any `a` and any `n`.

    void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride, fmpz_t xlength, const fmpz_t a, const fmpz_t b, const fmpz_t n)
    # Set `xstart`, `xstride`, and `xlength` so that the solution set for `x` modulo `n` in `a x = b \bmod n` is exactly `\{xstart + xstride\,i \mid 0 \le i < xlength\}`.
    # This function essentially gives a list of possibilities for the fraction `a/b` modulo `n`.
    # The outputs may not be aliased, and `n` should be positive.

    int fmpz_bit_pack(mp_limb_t * arr, flint_bitcnt_t shift, flint_bitcnt_t bits, const fmpz_t coeff, int negate, int borrow)
    # Shifts the given coefficient to the left by ``shift`` bits and adds
    # it to the integer in ``arr`` in a field of the given number of bits::
    # shift  bits  --------------
    # X X X C C C C 0 0 0 0 0 0 0
    # An optional borrow of `1` can be subtracted from ``coeff`` before
    # it is packed.  If ``coeff`` is negative after the borrow, then a
    # borrow will be returned by the function.
    # The value of ``shift`` is assumed to be less than ``FLINT_BITS``.
    # All but the first ``shift`` bits of ``arr`` are assumed to be zero
    # on entry to the function.
    # The value of ``coeff`` may also be optionally (and notionally) negated
    # before it is used, by setting the ``negate`` parameter to `-1`.

    int fmpz_bit_unpack(fmpz_t coeff, mp_limb_t * arr, flint_bitcnt_t shift, flint_bitcnt_t bits, int negate, int borrow)
    # A bit field of the given number of bits is extracted from ``arr``,
    # starting after ``shift`` bits, and placed into ``coeff``.  An
    # optional borrow of `1` may be added to the coefficient.  If the result
    # is negative, a borrow of `1` is returned.  Finally, the resulting
    # ``coeff`` may be negated by setting the ``negate`` parameter to `-1`.
    # The value of ``shift`` is expected to be less than ``FLINT_BITS``.

    void fmpz_bit_unpack_unsigned(fmpz_t coeff, const mp_limb_t * arr, flint_bitcnt_t shift, flint_bitcnt_t bits)
    # A bit field of the given number of bits is extracted from ``arr``,
    # starting after ``shift`` bits, and placed into ``coeff``.
    # The value of ``shift`` is expected to be less than ``FLINT_BITS``.

    void fmpz_complement(fmpz_t r, const fmpz_t f)
    # The variable ``r`` is set to the ones-complement of ``f``.

    void fmpz_clrbit(fmpz_t f, unsigned long i)
    # Sets the ``i``\th bit in ``f`` to zero.

    void fmpz_combit(fmpz_t f, unsigned long i)
    # Complements the ``i``\th bit in ``f``.

    void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b)
    # Sets ``r`` to the bit-wise logical ``and`` of ``a`` and ``b``.

    void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b)
    # Sets ``r`` to the bit-wise logical (inclusive) ``or`` of
    # ``a`` and ``b``.

    void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b)
    # Sets ``r`` to the bit-wise logical exclusive ``or`` of
    # ``a`` and ``b``.

    unsigned long fmpz_popcnt(const fmpz_t a)
    # Returns the number of '1' bits in the given Z (aka Hamming weight or
    # population count).
    # The return value is undefined if the input is negative.

    void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1, unsigned long r2, unsigned long m2, int sign)
    # Uses the Chinese Remainder Theorem to compute the unique integer
    # `0 \le x < M` (if sign = 0) or `-M/2 < x \le M/2` (if sign = 1)
    # congruent to `r_1` modulo `m_1` and `r_2` modulo `m_2`,
    # where `M = m_1 \times m_2`. The result `x` is stored in ``out``.
    # It is assumed that `m_1` and `m_2` are positive integers greater
    # than `1` and coprime.
    # If sign = 0, it is assumed that `0 \le r_1 < m_1` and `0 \le r_2 < m_2`.
    # Otherwise, it is assumed that `-m_1 \le r_1 < m_1` and `0 \le r_2 < m_2`.

    void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, fmpz_t r2, fmpz_t m2, int sign)
    # Use the Chinese Remainder Theorem to set ``out`` to the unique value
    # `0 \le x < M` (if sign = 0) or `-M/2 < x \le M/2` (if sign = 1)
    # congruent to `r_1` modulo `m_1` and `r_2` modulo `m_2`,
    # where `M = m_1 \times m_2`.
    # It is assumed that `m_1` and `m_2` are positive integers greater
    # than `1` and coprime.
    # If sign = 0, it is assumed that `0 \le r_1 < m_1` and `0 \le r_2 < m_2`.
    # Otherwise, it is assumed that `-m_1 \le r_1 < m_1` and `0 \le r_2 < m_2`.

    void fmpz_multi_mod_ui(mp_limb_t * out, const fmpz_t input, const fmpz_comb_t comb, fmpz_comb_temp_t temp)
    # Reduces the multiprecision integer ``in`` modulo each of the primes
    # stored in the ``comb`` structure. The array ``out`` will be filled
    # with the residues modulo these primes. The structure ``temp`` is
    # temporary space which must be provided by :func:`fmpz_comb_temp_init` and
    # cleared by :func:`fmpz_comb_temp_clear`.

    void fmpz_multi_CRT_ui(fmpz_t output, mp_srcptr residues, const fmpz_comb_t comb, fmpz_comb_temp_t ctemp, int sign)
    # This function takes a set of residues modulo the list of primes
    # contained in the ``comb`` structure and reconstructs a multiprecision
    # integer modulo the product of the primes which has
    # these residues modulo the corresponding primes.
    # If `N` is the product of all the primes then ``out`` is normalised to
    # be in the range `[0, N)` if sign = 0 and the range `[-(N-1)/2, N/2]`
    # if sign = 1. The array ``temp`` is temporary
    # space which must be provided by :func:`fmpz_comb_temp_init` and
    # cleared by :func:`fmpz_comb_temp_clear`.

    void fmpz_comb_init(fmpz_comb_t comb, mp_srcptr primes, long num_primes)
    # Initialises a ``comb`` structure for multimodular reduction and
    # recombination.  The array ``primes`` is assumed to contain
    # ``num_primes`` primes each of ``FLINT_BITS - 1`` bits. Modular
    # reductions and recombinations will be done modulo this list of primes.
    # The ``primes`` array must not be ``free``'d until the ``comb``
    # structure is no longer required and must be cleared by the user.

    void fmpz_comb_temp_init(fmpz_comb_temp_t temp, const fmpz_comb_t comb)
    # Creates temporary space to be used by multimodular and CRT functions
    # based on an initialised ``comb`` structure.

    void fmpz_comb_clear(fmpz_comb_t comb)
    # Clears the given ``comb`` structure, releasing any memory it uses.

    void fmpz_comb_temp_clear(fmpz_comb_temp_t temp)
    # Clears temporary space ``temp`` used by multimodular and CRT functions
    # using the given ``comb`` structure.

    void fmpz_multi_CRT_init(fmpz_multi_CRT_t CRT)
    # Initialize ``CRT`` for Chinese remaindering.

    int fmpz_multi_CRT_precompute(fmpz_multi_CRT_t CRT, const fmpz * moduli, long len)
    # Configure ``CRT`` for repeated Chinese remaindering of ``moduli``. The number of moduli, ``len``, should be positive.
    # A return of ``0`` indicates that the compilation failed and future
    # calls to :func:`fmpz_multi_CRT_precomp` will leave the output undefined.
    # A return of ``1`` indicates that the compilation was successful, which occurs if and only
    # if either (1) ``len == 1`` and ``modulus + 0`` is nonzero, or (2) no modulus is `0,1,-1` and all moduli are pairwise relatively prime.

    void fmpz_multi_CRT_precomp(fmpz_t output, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign)
    # Set ``output`` to an integer of smallest absolute value that is congruent to ``values + i`` modulo the ``moduli + i``
    # in ``P``.

    int fmpz_multi_CRT(fmpz_t output, const fmpz * moduli, const fmpz * values, long len, int sign)
    # Perform the same operation as :func:`fmpz_multi_CRT_precomp` while internally constructing and destroying the precomputed data.
    # All of the remarks in :func:`fmpz_multi_CRT_precompute` apply.

    void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P)
    # Free all space used by ``CRT``.

    int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t a)
    # Returns `1` if `n` is a strong probable prime to base `a`, otherwise it
    # returns `0`.

    int fmpz_is_probabprime_lucas(const fmpz_t n)
    # Performs a Lucas probable prime test with parameters chosen by Selfridge's
    # method `A` as per [BaiWag1980]_.
    # Return `1` if `n` is a Lucas probable prime, otherwise return `0`. This
    # function declares some composites probably prime, but no primes composite.

    int fmpz_is_probabprime_BPSW(const fmpz_t n)
    # Perform a Baillie-PSW probable prime test with parameters chosen by
    # Selfridge's method `A` as per [BaiWag1980]_.
    # Return `1` if `n` is a Lucas probable prime, otherwise return `0`.
    # There are no known composites passed as prime by this test, though
    # infinitely many probably exist. The test will declare no primes
    # composite.

    int fmpz_is_probabprime(const fmpz_t p)
    # Performs some trial division and then some probabilistic primality tests.
    # If `p` is definitely composite, the function returns `0`, otherwise it
    # is declared probably prime, i.e. prime for most practical purposes, and
    # the function returns `1`. The chance of declaring a composite prime is
    # very small.
    # Subsequent calls to the same function do not increase the probability of
    # the number being prime.

    int fmpz_is_prime_pseudosquare(const fmpz_t n)
    # Return `0` is `n` is composite. If `n` is too large (greater than about
    # `94` bits) the function fails silently and returns `-1`, otherwise, if
    # `n` is proven prime by the pseudosquares method, return `1`.
    # Tests if `n` is a prime according to [Theorem 2.7] [LukPatWil1996]_.
    # We first factor `N` using trial division up to some limit `B`.
    # In fact, the number of primes used in the trial factoring is at
    # most ``FLINT_PSEUDOSQUARES_CUTOFF``.
    # Next we compute `N/B` and find the next pseudosquare `L_p` above
    # this value, using a static table as per
    # https://oeis.org/A002189/b002189.txt.
    # As noted in the text, if `p` is prime then Step 3 will pass. This
    # test rejects many composites, and so by this time we suspect
    # that `p` is prime. If `N` is `3` or `7` modulo `8`, we are done,
    # and `N` is prime.
    # We now run a probable prime test, for which no known
    # counterexamples are known, to reject any composites. We then
    # proceed to prove `N` prime by executing Step 4. In the case that
    # `N` is `1` modulo `8`, if Step 4 fails, we extend the number of primes
    # `p_i` at Step 3 and hope to find one which passes Step 4. We take
    # the test one past the largest `p` for which we have pseudosquares
    # `L_p` tabulated, as this already corresponds to the next `L_p` which
    # is bigger than `2^{64}` and hence larger than any prime we might be
    # testing.
    # As explained in the text, Condition 4 cannot fail if `N` is prime.
    # The possibility exists that the probable prime test declares a
    # composite prime. However in that case an error is printed, as
    # that would be of independent interest.

    int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, long num_pm1)
    # Applies the Pocklington primality test. The test computes a product
    # `F` of prime powers which divide `n - 1`.
    # The function then returns either `0` if `n` is definitely composite
    # or it returns `1` if all factors of `n` are `1 \pmod{F}`. Also in
    # that case, `R` is set to `(n - 1)/F`.
    # NB: a return value of `1` only proves `n` prime if `F \ge \sqrt{n}`.
    # The function does not compute which primes divide `n - 1`. Instead,
    # these must be supplied as an array ``pm1`` of length ``num_pm1``.
    # It does not matter how many prime factors are supplied, but the more
    # that are supplied, the larger F will be.
    # There is a balance between the amount of time spent looking for
    # factors of `n - 1` and the usefulness of the output (`F` may be as low
    # as `2` in some cases).
    # A reasonable heuristic seems to be to choose ``limit`` to be some
    # small multiple of `\log^3(n)/10` (e.g. `1, 2, 5` or `10`) depending
    # on how long one is prepared to wait, then to trial factor up to the
    # limit. (See ``_fmpz_nm1_trial_factors``.)
    # Requires `n` to be odd.

    void _fmpz_nm1_trial_factors(const fmpz_t n, mp_ptr pm1, long * num_pm1, unsigned long limit)
    # Trial factors `n - 1` up to the given limit (approximately) and stores
    # the factors in an array ``pm1`` whose length is written out to
    # ``num_pm1``.
    # One can use `\log(n) + 2` as a bound on the number of factors which might
    # be produced (and hence on the length of the array that needs to be
    # supplied).

    int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pp1, long num_pp1)
    # Applies the Morrison `p + 1` primality test. The test computes a
    # product `F` of primes which divide `n + 1`.
    # The function then returns either `0` if `n` is definitely composite
    # or it returns `1` if all factors of `n` are `\pm 1 \pmod{F}`. Also in
    # that case, `R` is set to `(n + 1)/F`.
    # NB: a return value of `1` only proves `n` prime if
    # `F > \sqrt{n} + 1`.
    # The function does not compute which primes divide `n + 1`. Instead,
    # these must be supplied as an array ``pp1`` of length ``num_pp1``.
    # It does not matter how many prime factors are supplied, but the more
    # that are supplied, the larger `F` will be.
    # There is a balance between the amount of time spent looking for
    # factors of `n + 1` and the usefulness of the output (`F` may be as low
    # as `2` in some cases).
    # A reasonable heuristic seems to be to choose ``limit`` to be some
    # small multiple of `\log^3(n)/10` (e.g. `1, 2, 5` or `10`) depending
    # on how long one is prepared to wait, then to trial factor up to the
    # limit. (See ``_fmpz_np1_trial_factors``.)
    # Requires `n` to be odd and non-square.

    void _fmpz_np1_trial_factors(const fmpz_t n, mp_ptr pp1, long * num_pp1, unsigned long limit)
    # Trial factors `n + 1` up to the given limit (approximately) and stores
    # the factors in an array ``pp1`` whose length is written out to
    # ``num_pp1``.
    # One can use `\log(n) + 2` as a bound on the number of factors which might
    # be produced (and hence on the length of the array that needs to be
    # supplied).

    int fmpz_is_prime(const fmpz_t n)
    # Attempts to prove `n` prime.  If `n` is proven prime, the function
    # returns `1`. If `n` is definitely composite, the function returns `0`.
    # This function calls :func:`n_is_prime` for `n` that fits in a single word.
    # For `n` larger than one word, it tests divisibility by a few small primes
    # and whether `n` is a perfect square to rule out trivial composites.
    # For `n` up to about 81 bits, it then uses a strong probable prime test
    # (Miller-Rabin test) with the first 13 primes as witnesses. This has
    # been shown to prove primality [SorWeb2016]_.
    # For larger `n`, it does a single base-2 strong probable prime test
    # to eliminate most composite numbers. If `n` passes, it does a
    # combination of Pocklington, Morrison and Brillhart, Lehmer, Selfridge
    # tests. If any of these tests fails to give a proof, it falls back to
    # performing an APRCL test.
    # The APRCL test could theoretically fail to prove that `n` is prime
    # or composite. In that case, the program aborts. This is not expected to
    # occur in practice.

    void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t m, const fmpz_t n)
    # Given `V_0 = 2`, `V_1 = A` compute `V_m, V_{m + 1} \pmod{n}` from the
    # recurrences `V_j = AV_{j - 1} - V_{j - 2} \pmod{n}`.
    # This is computed efficiently using `V_{2j} = V_j^2 - 2 \pmod{n}` and
    # `V_{2j + 1} = V_jV_{j + 1} - A \pmod{n}`.
    # No aliasing is permitted.

    void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t m, const fmpz_t n)
    # Given `V_0 = 2`, `V_1 = A` compute `V_m, V_{m + 1} \pmod{n}` from the
    # recurrences `V_j = AV_{j - 1} - BV_{j - 2} \pmod{n}`.
    # This is computed efficiently using double and add formulas.
    # No aliasing is permitted.

    void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t n)
    # Given `U_m, U_{m + 1} \pmod{n}` compute `U_{2m}, U_{2m + 1} \pmod{n}`.
    # Aliasing of `U_{2m}` and `U_m` and aliasing of `U_{2m + 1}` and `U_{m + 1}`
    # is permitted. No other aliasing is allowed.

    void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t Un, const fmpz_t Un1, const fmpz_t A, const fmpz_t B, const fmpz_t n)
    # Given `U_m, U_{m + 1} \pmod{n}` and `U_n, U_{n + 1} \pmod{n}` compute
    # `U_{m + n}, U_{m + n + 1} \pmod{n}`.
    # Aliasing of `U_{m + n}` with `U_m` or `U_n` and aliasing of `U_{m + n + 1}`
    # with `U_{m + 1}` or `U_{n + 1}` is permitted. No other aliasing is allowed.

    void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t k, const fmpz_t n)
    # Given `U_m, U_{m + 1} \pmod{n}` compute `U_{km}, U_{km + 1} \pmod{n}`.
    # Aliasing of `U_{km}` and `U_m` and aliasing of `U_{km + 1}` and `U_{m + 1}`
    # is permitted. No other aliasing is allowed.

    void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, const fmpz_t Vm, const fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, const fmpz_t n)
    # Given `V_m, V_{m + 1} \pmod{n}` compute `U_m, U_{m + 1} \pmod{n}`.
    # Aliasing of `V_m` and `U_m` and aliasing of `V_{m + 1}` and `U_{m + 1}`
    # is permitted. No other aliasing is allowed.

    int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, const fmpz_t r, const fmpz_t s)
    # If there exists a proper divisor of `n` which is `r \pmod{s}` for
    # `0 < r < s < n`, this function returns `1` and sets ``fac`` to such a
    # divisor. Otherwise the function returns `0` and the value of ``fac`` is
    # undefined.
    # We require `\gcd(r, s) = 1`.
    # This is efficient if `s^3 > n`.

    void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved)
    # Finds the next prime number larger than `n`.
    # If ``proved`` is nonzero, then the integer returned is
    # guaranteed to actually be prime. Otherwise if `n` fits in
    # ``FLINT_BITS - 3`` bits ``n_nextprime`` is called, and if not then
    # the GMP ``mpz_nextprime`` function is called. Up to and including
    # GMP 6.1.2 this used Miller-Rabin iterations, and thereafter uses
    # a BPSW test.

    void fmpz_primorial(fmpz_t res, unsigned long n)
    # Sets ``res`` to ``n`` primorial or `n \#`, the product of all prime
    # numbers less than or equal to `n`.

    void fmpz_factor_euler_phi(fmpz_t res, const fmpz_factor_t fac)
    void fmpz_euler_phi(fmpz_t res, const fmpz_t n)
    # Sets ``res`` to the Euler totient function `\phi(n)`, counting the
    # number of positive integers less than or equal to `n` that are coprime
    # to `n`. The factor version takes a precomputed
    # factorisation of `n`.

    int fmpz_factor_moebius_mu(const fmpz_factor_t fac)
    int fmpz_moebius_mu(const fmpz_t n)
    # Computes the Moebius function `\mu(n)`, which is defined as `\mu(n) = 0`
    # if `n` has a prime factor of multiplicity greater than `1`, `\mu(n) = -1`
    # if `n` has an odd number of distinct prime factors, and `\mu(n) = 1` if
    # `n` has an even number of distinct prime factors.  By convention,
    # `\mu(0) = 0`. The factor version takes a precomputed
    # factorisation of `n`.

    void fmpz_factor_divisor_sigma(fmpz_t res, unsigned long k, const fmpz_factor_t fac)
    void fmpz_divisor_sigma(fmpz_t res, unsigned long k, const fmpz_t n)
    # Sets ``res`` to `\sigma_k(n)`, the sum of `k`\th powers of all
    # divisors of `n`. The factor version takes a precomputed
    # factorisation of `n`.
