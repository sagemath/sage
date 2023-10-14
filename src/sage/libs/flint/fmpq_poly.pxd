# distutils: libraries = flint
# distutils: depends = flint/fmpq_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpq_poly_init(fmpq_poly_t poly)
    # Initialises the polynomial for use.  The length is set to zero.

    void fmpq_poly_init2(fmpq_poly_t poly, slong alloc)
    # Initialises the polynomial with space for at least ``alloc``
    # coefficients and sets the length to zero. The ``alloc`` coefficients
    # are all set to zero.

    void fmpq_poly_realloc(fmpq_poly_t poly, slong alloc)
    # Reallocates the given polynomial to have space for ``alloc``
    # coefficients. If ``alloc`` is zero then the polynomial is cleared
    # and then reinitialised.  If the current length is greater than
    # ``alloc`` then ``poly`` is first truncated to length
    # ``alloc``. Note that this might leave the rational polynomial in
    # non-canonical form.

    void fmpq_poly_fit_length(fmpq_poly_t poly, slong len)
    # If ``len`` is greater than the number of coefficients currently
    # allocated, then the polynomial is reallocated to have space for at
    # least ``len`` coefficients. No data is lost when calling this
    # function. The function efficiently deals with the case where
    # :func:`fit_length` is called many times in small increments by at
    # least doubling the number of allocated coefficients when ``len``
    # is larger than the number of coefficients currently allocated.

    void _fmpq_poly_set_length(fmpq_poly_t poly, slong len)
    # Sets the length of the numerator polynomial to ``len``, demoting
    # coefficients beyond the new length.  Note that this method does
    # not guarantee that the rational polynomial is in canonical form.

    void fmpq_poly_clear(fmpq_poly_t poly)
    # Clears the given polynomial, releasing any memory used. The polynomial
    # must be reinitialised in order to be used again.

    void _fmpq_poly_normalise(fmpq_poly_t poly)
    # Sets the length of ``poly`` so that the top coefficient is
    # non-zero. If all coefficients are zero, the length is set to zero.
    # Note that this function does not guarantee the coprimality of the
    # numerator polynomial and the integer denominator.

    void _fmpq_poly_canonicalise(fmpz * poly, fmpz_t den, slong len)
    # Puts ``(poly, den)`` of length ``len`` into canonical form.
    # It is assumed that the array ``poly`` contains a non-zero entry in
    # position ``len - 1`` whenever ``len > 0``.  Assumes that ``den``
    # is non-zero.

    void fmpq_poly_canonicalise(fmpq_poly_t poly)
    # Puts the polynomial ``poly`` into canonical form.  Firstly, the length
    # is set to the actual length of the numerator polynomial.  For non-zero
    # polynomials, it is then ensured that the numerator and denominator are
    # coprime and that the denominator is positive.  The canonical form of the
    # zero polynomial is a zero numerator polynomial and a one denominator.

    int _fmpq_poly_is_canonical(const fmpz * poly, const fmpz_t den, slong len)
    # Returns whether the polynomial is in canonical form.

    int fmpq_poly_is_canonical(const fmpq_poly_t poly)
    # Returns whether the polynomial is in canonical form.

    slong fmpq_poly_degree(const fmpq_poly_t poly)
    # Returns the degree of ``poly``, which is one less than its length, as
    # a ``slong``.

    slong fmpq_poly_length(const fmpq_poly_t poly)
    # Returns the length of ``poly``.

    fmpz * fmpq_poly_numref(fmpq_poly_t poly)
    # Returns a reference to the numerator polynomial as an array.
    # Note that, because of a delayed initialisation approach, this might
    # be ``NULL`` for zero polynomials.  This situation can be salvaged
    # by calling either :func:`fmpq_poly_fit_length` or
    # :func:`fmpq_poly_realloc`.
    # This function is implemented as a macro returning ``(poly)->coeffs``.

    fmpz_t fmpq_poly_denref(fmpq_poly_t poly)
    # Returns a reference to the denominator as a ``fmpz_t``.  The integer
    # is guaranteed to be properly initialised.
    # This function is implemented as a macro returning ``(poly)->den``.

    void fmpq_poly_get_numerator(fmpz_poly_t res, const fmpq_poly_t poly)
    # Sets ``res`` to the numerator of ``poly``, e.g. the primitive part
    # as an ``fmpz_poly_t`` if it is in canonical form.

    void fmpq_poly_get_denominator(fmpz_t den, const fmpq_poly_t poly)
    # Sets ``res`` to the denominator of ``poly``.

    void fmpq_poly_randtest(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # Sets `f` to a random polynomial with coefficients up to the given
    # length and where each coefficient has up to the given number of bits.
    # The coefficients are signed randomly.  One must call
    # :func:`flint_randinit` before calling this function.

    void fmpq_poly_randtest_unsigned(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # Sets `f` to a random polynomial with coefficients up to the given length
    # and where each coefficient has up to the given number of bits.  One must
    # call :func:`flint_randinit` before calling this function.

    void fmpq_poly_randtest_not_zero(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # As for :func:`fmpq_poly_randtest` except that ``len`` and ``bits``
    # may not be zero and the polynomial generated is guaranteed not to be the
    # zero polynomial.  One must call :func:`flint_randinit` before calling
    # this function.

    void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``poly1`` to equal ``poly2``.

    void fmpq_poly_set_si(fmpq_poly_t poly, slong x)
    # Sets ``poly`` to the integer `x`.

    void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x)
    # Sets ``poly`` to the integer `x`.

    void fmpq_poly_set_fmpz(fmpq_poly_t poly, const fmpz_t x)
    # Sets ``poly`` to the integer `x`.

    void fmpq_poly_set_fmpq(fmpq_poly_t poly, const fmpq_t x)
    # Sets ``poly`` to the rational `x`, which is assumed to be
    # given in lowest terms.

    void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, const fmpz_poly_t op)
    # Sets the rational polynomial ``rop`` to the same value
    # as the integer polynomial ``op``.

    void fmpq_poly_set_nmod_poly(fmpq_poly_t rop, const nmod_poly_t op)
    # Sets the coefficients of ``rop`` to the residues in ``op``,
    # normalised to the interval `-m/2 \le r < m/2` where `m` is the modulus.

    void fmpq_poly_get_nmod_poly(nmod_poly_t rop, const fmpq_poly_t op)
    # Sets the coefficients of ``rop`` to the coefficients in the denominator of ``op``,
    # reduced by the modulus of ``rop``. The result is multiplied by the inverse of the
    # denominator of ``op``. It is assumed that the reduction of the denominator of ``op``
    # is invertible.

    void fmpq_poly_get_nmod_poly_den(nmod_poly_t rop, const fmpq_poly_t op, int den)
    # Sets the coefficients of ``rop`` to the coefficients in the denominator
    # of ``op``, reduced by the modulus of ``rop``. If ``den == 1``, the result is
    # multiplied by the inverse of the denominator of ``op``. In this case it is
    # assumed that the reduction of the denominator of ``op`` is invertible.

    int _fmpq_poly_set_str(fmpz * poly, fmpz_t den, const char * str, slong len)
    # Sets ``(poly, den)`` to the polynomial specified by the
    # null-terminated string ``str`` of ``len`` coefficients. The input
    # format is a sequence of coefficients separated by one space.
    # The result is only guaranteed to be in lowest terms if all
    # coefficients in the input string are in lowest terms.
    # Returns `0` if no error occurred. Otherwise, returns -1
    # in which case the resulting value of ``(poly, den)`` is undefined.
    # If ``str`` is not null-terminated, calling this method might result
    # in a segmentation fault.

    int fmpq_poly_set_str(fmpq_poly_t poly, const char * str)
    # Sets ``poly`` to the polynomial specified by the null-terminated
    # string ``str``. The input format is the same as the output format
    # of ``fmpq_poly_get_str``: the length given as a decimal integer,
    # then two spaces, then the list of coefficients separated by one space.
    # The result is only guaranteed to be in canonical form if all
    # coefficients in the input string are in lowest terms.
    # Returns `0` if no error occurred.  Otherwise, returns -1 in which case
    # the resulting value of ``poly`` is set to zero. If ``str`` is not
    # null-terminated, calling this method might result in a segmentation fault.

    char * fmpq_poly_get_str(const fmpq_poly_t poly)
    # Returns the string representation of ``poly``.

    char * fmpq_poly_get_str_pretty(const fmpq_poly_t poly, const char * var)
    # Returns the pretty representation of ``poly``, using the
    # null-terminated string ``var`` not equal to ``"\0"`` as
    # the variable name.

    void fmpq_poly_zero(fmpq_poly_t poly)
    # Sets ``poly`` to zero.

    void fmpq_poly_one(fmpq_poly_t poly)
    # Sets ``poly`` to the constant polynomial `1`.

    void fmpq_poly_neg(fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``poly1`` to the additive inverse of ``poly2``.

    void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``poly1`` to the multiplicative inverse of ``poly2``
    # if possible.  Otherwise, if ``poly2`` is not a unit, leaves
    # ``poly1`` unmodified and calls :func:`abort`.

    void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2)
    # Efficiently swaps the polynomials ``poly1`` and ``poly2``.

    void fmpq_poly_truncate(fmpq_poly_t poly, slong n)
    # If the current length of ``poly`` is greater than `n`, it is
    # truncated to the given length.  Discarded coefficients are demoted,
    # but they are not necessarily set to zero.

    void fmpq_poly_set_trunc(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Sets ``res`` to a copy of ``poly``, truncated to length ``n``.

    void fmpq_poly_get_slice(fmpq_poly_t rop, const fmpq_poly_t op, slong i, slong j)
    # Returns the slice with coefficients from `x^i` (including) to
    # `x^j` (excluding).

    void fmpq_poly_reverse(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # This function considers the polynomial ``poly`` to be of length `n`,
    # notionally truncating and zero padding if required, and reverses
    # the result.  Since the function normalises its result ``res`` may be
    # of length less than `n`.

    void fmpq_poly_get_coeff_fmpz(fmpz_t x, const fmpq_poly_t poly, slong n)
    # Retrieves the `n`\th coefficient of the numerator of ``poly``.

    void fmpq_poly_get_coeff_fmpq(fmpq_t x, const fmpq_poly_t poly, slong n)
    # Retrieves the `n`\th coefficient of ``poly``, in lowest terms.

    void fmpq_poly_set_coeff_si(fmpq_poly_t poly, slong n, slong x)
    # Sets the `n`\th coefficient in ``poly`` to the integer `x`.

    void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, slong n, ulong x)
    # Sets the `n`\th coefficient in ``poly`` to the integer `x`.

    void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, slong n, const fmpz_t x)
    # Sets the `n`\th coefficient in ``poly`` to the integer `x`.

    void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, slong n, const fmpq_t x)
    # Sets the `n`\th coefficient in ``poly`` to the rational `x`.

    int fmpq_poly_equal(const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Returns `1` if ``poly1`` is equal to ``poly2``,
    # otherwise returns `0`.

    int _fmpq_poly_equal_trunc(const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # Returns `1` if ``poly1`` and ``poly2`` notionally truncated to length
    # `n` are equal, otherwise returns `0`.

    int fmpq_poly_equal_trunc(const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # Returns `1` if ``poly1`` and ``poly2`` notionally truncated to length
    # `n` are equal, otherwise returns `0`.

    int _fmpq_poly_cmp(const fmpz * lpoly, const fmpz_t lden, const fmpz * rpoly, const fmpz_t rden, slong len)
    # Compares two non-zero polynomials, assuming they have the same length
    # ``len > 0``.
    # The polynomials are expected to be provided in canonical form.

    int fmpq_poly_cmp(const fmpq_poly_t left, const fmpq_poly_t right)
    # Compares the two polynomials ``left`` and ``right``.
    # Compares the two polynomials ``left`` and ``right``, returning
    # `-1`, `0`, or `1` as ``left`` is less than, equal to, or greater
    # than ``right``.  The comparison is first done by the degree, and
    # then, in case of a tie, by the individual coefficients from highest
    # to lowest.

    int fmpq_poly_is_one(const fmpq_poly_t poly)
    # Returns `1` if ``poly`` is the constant polynomial `1`, otherwise
    # returns `0`.

    int fmpq_poly_is_zero(const fmpq_poly_t poly)
    # Returns `1` if ``poly`` is the zero polynomial, otherwise returns `0`.

    int fmpq_poly_is_gen(const fmpq_poly_t poly)
    # Returns `1` if ``poly`` is the degree `1` polynomial `x`, otherwise returns
    # `0`.

    void _fmpq_poly_add(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
    # Forms the sum ``(rpoly, rden)`` of ``(poly1, den1, len1)`` and
    # ``(poly2, den2, len2)``, placing the result into canonical form.
    # Assumes that ``rpoly`` is an array of length the maximum of
    # ``len1`` and ``len2``.  The input operands are assumed to
    # be in canonical form and are also allowed to be of length `0`.
    # ``(rpoly, rden)`` and ``(poly1, den1)`` may be aliased,
    # but ``(rpoly, rden)`` and ``(poly2, den2)`` may *not*
    # be aliased.

    void _fmpq_poly_add_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, int can)
    # As per ``_fmpq_poly_add`` except that one can specify whether to
    # canonicalise the output or not. This function is intended to be used with
    # weak canonicalisation to prevent explosion in memory usage. It exists for
    # performance reasons.

    void fmpq_poly_add(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``res`` to the sum of ``poly1`` and ``poly2``, using
    # Henrici's algorithm.

    void fmpq_poly_add_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can)
    # As per ``mpq_poly_add`` except that one can specify whether to
    # canonicalise the output or not. This function is intended to be used with
    # weak canonicalisation to prevent explosion in memory usage. It exists for
    # performance reasons.

    void _fmpq_poly_add_series(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # As per ``_fmpq_poly_add`` but the inputs are first notionally truncated
    # to length `n`. If `n` is less than ``len1`` or ``len2`` then the
    # output only needs space for `n` coefficients. We require `n \geq 0`.

    void _fmpq_poly_add_series_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n, int can)
    # As per ``_fmpq_poly_add_can`` but the inputs are first notionally
    # truncated to length `n`. If `n` is less than ``len1`` or ``len2``
    # then the output only needs space for `n` coefficients. We require
    # `n \geq 0`.

    void fmpq_poly_add_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # As per ``fmpq_poly_add`` but the inputs are first notionally
    # truncated to length `n`.

    void fmpq_poly_add_series_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can)
    # As per ``fmpq_poly_add_can`` but the inputs are first notionally
    # truncated to length `n`.

    void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
    # Forms the difference ``(rpoly, rden)`` of ``(poly1, den1, len1)``
    # and ``(poly2, den2, len2)``, placing the result into canonical form.
    # Assumes that ``rpoly`` is an array of length the maximum of
    # ``len1`` and ``len2``.  The input operands are assumed to be in
    # canonical form and are also allowed to be of length `0`.
    # ``(rpoly, rden)`` and ``(poly1, den1, len1)`` may be aliased,
    # but ``(rpoly, rden)`` and ``(poly2, den2, len2)`` may *not* be
    # aliased.

    void _fmpq_poly_sub_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, int can)
    # As per ``_fmpq_poly_sub`` except that one can specify whether to
    # canonicalise the output or not. This function is intended to be used with
    # weak canonicalisation to prevent explosion in memory usage. It exists for
    # performance reasons.

    void fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``res`` to the difference of ``poly1`` and ``poly2``,
    # using Henrici's algorithm.

    void fmpq_poly_sub_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can)
    # As per ``_fmpq_poly_sub`` except that one can specify whether to
    # canonicalise the output or not. This function is intended to be used with
    # weak canonicalisation to prevent explosion in memory usage. It exists for
    # performance reasons.

    void _fmpq_poly_sub_series(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # As per ``_fmpq_poly_sub`` but the inputs are first notionally truncated
    # to length `n`. If `n` is less than ``len1`` or ``len2`` then the
    # output only needs space for `n` coefficients. We require `n \geq 0`.

    void _fmpq_poly_sub_series_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n, int can)
    # As per ``_fmpq_poly_sub_can`` but the inputs are first notionally
    # truncated to length `n`. If `n` is less than ``len1`` or ``len2``
    # then the output only needs space for `n` coefficients. We require
    # `n \geq 0`.

    void fmpq_poly_sub_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # As per ``fmpq_poly_sub`` but the inputs are first notionally
    # truncated to length `n`.

    void fmpq_poly_sub_series_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can)
    # As per ``fmpq_poly_sub_can`` but the inputs are first notionally
    # truncated to length `n`.

    void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong c)
    # Sets ``(rpoly, rden, len)`` to the product of `c` of
    # ``(poly, den, len)``.
    # If the input is normalised, then so is the output, provided it is
    # non-zero.  If the input is in lowest terms, then so is the output.
    # However, even if neither of these conditions are met, the result
    # will be (mathematically) correct.
    # Supports exact aliasing between ``(rpoly, den)``
    # and ``(poly, den)``.

    void _fmpq_poly_scalar_mul_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, ulong c)
    # Sets ``(rpoly, rden, len)`` to the product of `c` of
    # ``(poly, den, len)``.
    # If the input is normalised, then so is the output, provided it is
    # non-zero.  If the input is in lowest terms, then so is the output.
    # However, even if neither of these conditions are met, the result
    # will be (mathematically) correct.
    # Supports exact aliasing between ``(rpoly, den)``
    # and ``(poly, den)``.

    void _fmpq_poly_scalar_mul_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t c)
    # Sets ``(rpoly, rden, len)`` to the product of `c` of
    # ``(poly, den, len)``.
    # If the input is normalised, then so is the output, provided it is
    # non-zero.  If the input is in lowest terms, then so is the output.
    # However, even if neither of these conditions are met, the result
    # will be (mathematically) correct.
    # Supports exact aliasing between ``(rpoly, den)``
    # and ``(poly, den)``.

    void _fmpq_poly_scalar_mul_fmpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s)
    # Sets ``(rpoly, rden)`` to the product of `r/s` and
    # ``(poly, den, len)``, in lowest terms.
    # Assumes that ``(poly, den, len)`` and `r/s` are provided in lowest
    # terms.  Assumes that ``rpoly`` is an array of length ``len``.
    # Supports aliasing of ``(rpoly, den)`` and ``(poly, den)``.
    # The ``fmpz_t``'s `r` and `s` may not be part of ``(rpoly, rden)``.

    void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
    # Sets ``rop`` to `c` times ``op``.

    void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
    # Sets ``rop`` to `c` times ``op``.

    void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
    # Sets ``rop`` to `c` times ``op``.  Assumes that the ``fmpz_t c``
    # is not part of ``rop``.

    void fmpq_poly_scalar_mul_mpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
    # Sets ``rop`` to `c` times ``op``.

    void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t c)
    # Sets ``(rpoly, rden, len)`` to ``(poly, den, len)`` divided by `c`,
    # in lowest terms.
    # Assumes that ``len`` is positive.  Assumes that `c` is non-zero.
    # Supports aliasing between ``(rpoly, rden)`` and ``(poly, den)``.
    # Assumes that `c` is not part of ``(rpoly, rden)``.

    void _fmpq_poly_scalar_div_si(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong c)
    # Sets ``(rpoly, rden, len)`` to ``(poly, den, len)`` divided by `c`,
    # in lowest terms.
    # Assumes that ``len`` is positive.  Assumes that `c` is non-zero.
    # Supports aliasing between ``(rpoly, rden)`` and ``(poly, den)``.

    void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, ulong c)
    # Sets ``(rpoly, rden, len)`` to ``(poly, den, len)`` divided by `c`,
    # in lowest terms.
    # Assumes that ``len`` is positive.  Assumes that `c` is non-zero.
    # Supports aliasing between ``(rpoly, rden)`` and ``(poly, den)``.

    void _fmpq_poly_scalar_div_fmpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s)
    # Sets ``(rpoly, rden, len)`` to ``(poly, den, len)`` divided by `r/s`,
    # in lowest terms.
    # Assumes that ``len`` is positive.  Assumes that `r/s` is non-zero and
    # in lowest terms.  Supports aliasing between ``(rpoly, rden)`` and
    # ``(poly, den)``. The ``fmpz_t``'s `r` and `s` may not be part of
    # ``(rpoly, poly)``.

    void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
    void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
    void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
    void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
    # Sets ``rop`` to ``op`` divided by the scalar ``c``.

    void _fmpq_poly_mul(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
    # Sets ``(rpoly, rden, len1 + len2 - 1)`` to the product of
    # ``(poly1, den1, len1)`` and ``(poly2, den2, len2)``. If the
    # input is provided in canonical form, then so is the output.
    # Assumes ``len1 >= len2 > 0``.  Allows zero-padding in the input.
    # Does not allow aliasing between the inputs and outputs.

    void fmpq_poly_mul(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``res`` to the product of ``poly1`` and ``poly2``.

    void _fmpq_poly_mullow(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # Sets ``(rpoly, rden, n)`` to the low `n` coefficients of
    # ``(poly1, den1)`` and ``(poly2, den2)``.  The output is
    # not guaranteed to be in canonical form.
    # Assumes ``len1 >= len2 > 0`` and ``0 < n <= len1 + len2 - 1``.
    # Allows for zero-padding in the inputs.  Does not allow aliasing between
    # the inputs and outputs.

    void fmpq_poly_mullow(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # Sets ``res`` to the product of ``poly1`` and ``poly2``,
    # truncated to length `n`.

    void fmpq_poly_addmul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)
    # Adds the product of ``op1`` and ``op2`` to ``rop``.

    void fmpq_poly_submul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)
    # Subtracts the product of ``op1`` and ``op2`` from ``rop``.

    void _fmpq_poly_pow(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, ulong e)
    # Sets ``(rpoly, rden)`` to ``(poly, den)^e``, assuming
    # ``e, len > 0``.  Assumes that ``rpoly`` is an array of
    # length at least ``e * (len - 1) + 1``.  Supports aliasing
    # of ``(rpoly, den)`` and ``(poly, den)``.

    void fmpq_poly_pow(fmpq_poly_t res, const fmpq_poly_t poly, ulong e)
    # Sets ``res`` to ``poly^e``, where the only special case `0^0` is
    # defined as `1`.

    void _fmpq_poly_pow_trunc(fmpz * res, fmpz_t rden, const fmpz * f, const fmpz_t fden, slong flen, ulong exp, slong len)
    # Sets ``(rpoly, rden, len)`` to ``(poly, den)^e`` truncated to length ``len``,
    # where ``len`` is at most ``e * (flen - 1) + 1``.

    void fmpq_poly_pow_trunc(fmpq_poly_t res, const fmpq_poly_t poly, ulong e, slong n)
    # Sets ``res`` to ``poly^e`` truncated to length ``n``.

    void fmpq_poly_shift_left(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Set ``res`` to ``poly`` shifted left by `n` coefficients. Zero
    # coefficients are inserted.

    void fmpq_poly_shift_right(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Set ``res`` to ``poly`` shifted right by `n` coefficients.
    # If `n` is equal to or greater than the current length of ``poly``,
    # ``res`` is set to the zero polynomial.

    void _fmpq_poly_divrem(fmpz * Q, fmpz_t q, fmpz * R, fmpz_t r, const fmpz * A, const fmpz_t a, slong lenA, const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
    # Finds the quotient ``(Q, q)`` and remainder ``(R, r)`` of the
    # Euclidean division of ``(A, a)`` by ``(B, b)``.
    # Assumes that ``lenA >= lenB > 0``.  Assumes that `R` has space for
    # ``lenA`` coefficients, although only the bottom ``lenB - 1`` will
    # carry meaningful data on exit.  Supports no aliasing between the two
    # outputs, or between the inputs and the outputs.
    # An optional precomputed inverse of the leading coefficient of `B` from
    # ``fmpz_preinvn_init`` can be supplied. Otherwise ``inv`` should be
    # ``NULL``.
    # Note: ``fmpz.h`` has to be included before ``fmpq_poly.h`` in order for the
    # latter to declare this function.

    void fmpq_poly_divrem(fmpq_poly_t Q, fmpq_poly_t R, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Finds the quotient `Q` and remainder `R` of the Euclidean division of
    # ``poly1`` by ``poly2``.

    void _fmpq_poly_div(fmpz * Q, fmpz_t q, const fmpz * A, const fmpz_t a, slong lenA, const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
    # Finds the quotient ``(Q, q)`` of the Euclidean division
    # of ``(A, a)`` by ``(B, b)``.
    # Assumes that ``lenA >= lenB > 0``.  Supports no aliasing
    # between the inputs and the outputs.
    # An optional precomputed inverse of the leading coefficient of `B` from
    # ``fmpz_preinvn_init`` can be supplied. Otherwise ``inv`` should be
    # ``NULL``.
    # Note: ``fmpz.h`` has to be included before ``fmpq_poly.h`` in order for the
    # latter to declare this function.

    void fmpq_poly_div(fmpq_poly_t Q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Finds the quotient `Q` and remainder `R` of the Euclidean division
    # of ``poly1`` by ``poly2``.

    void _fmpq_poly_rem(fmpz * R, fmpz_t r, const fmpz * A, const fmpz_t a, slong lenA, const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
    # Finds the remainder ``(R, r)`` of the Euclidean division
    # of ``(A, a)`` by ``(B, b)``.
    # Assumes that ``lenA >= lenB > 0``.  Supports no aliasing between
    # the inputs and the outputs.
    # An optional precomputed inverse of the leading coefficient of `B` from
    # ``fmpz_preinvn_init`` can be supplied. Otherwise ``inv`` should be
    # ``NULL``.
    # Note: ``fmpz.h`` has to be included before ``fmpq_poly.h`` in order for the
    # latter to declare this function.

    void fmpq_poly_rem(fmpq_poly_t R, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Finds the remainder `R` of the Euclidean division
    # of ``poly1`` by ``poly2``.

    fmpq_poly_struct * _fmpq_poly_powers_precompute(const fmpz * B, const fmpz_t denB, slong len)
    # Computes ``2*len - 1`` powers of `x` modulo the polynomial `B` of
    # the given length. This is used as a kind of precomputed inverse in
    # the remainder routine below.

    void fmpq_poly_powers_precompute(fmpq_poly_powers_precomp_t pinv, fmpq_poly_t poly)
    # Computes ``2*len - 1`` powers of `x` modulo the polynomial `B` of the given length. This is used as a kind of precomputed inverse in the remainder routine below.

    void _fmpq_poly_powers_clear(fmpq_poly_struct * powers, slong len)
    # Clean up resources used by precomputed powers which have been computed
    # by ``_fmpq_poly_powers_precompute``.

    void fmpq_poly_powers_clear(fmpq_poly_powers_precomp_t pinv)
    # Clean up resources used by precomputed powers which have been computed
    # by ``fmpq_poly_powers_precompute``.

    void _fmpq_poly_rem_powers_precomp(fmpz * A, fmpz_t denA, slong m, const fmpz * B, const fmpz_t denB, slong n, fmpq_poly_struct * const powers)
    # Set `A` to the remainder of `A` divide `B` given precomputed powers mod `B`
    # provided by ``_fmpq_poly_powers_precompute``. No aliasing is allowed.
    # This function is only faster if `m \leq 2\cdot n - 1`.
    # The output of this function is *not* canonicalised.

    void fmpq_poly_rem_powers_precomp(fmpq_poly_t R, const fmpq_poly_t A, const fmpq_poly_t B, const fmpq_poly_powers_precomp_t B_inv)
    # Set `R` to the remainder of `A` divide `B` given precomputed powers mod `B`
    # provided by ``fmpq_poly_powers_precompute``.
    # This function is only faster if ``A->length <= 2*B->length - 1``.
    # The output of this function is *not* canonicalised.

    int _fmpq_poly_divides(fmpz * qpoly, fmpz_t qden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
    # Return `1` if ``(poly2, den2, len2)`` divides ``(poly1, den1, len1)`` and
    # set ``(qpoly, qden, len1 - len2 + 1)`` to the quotient. Otherwise return
    # `0`. Requires that ``qpoly`` has space for ``len1 - len2 + 1``
    # coefficients and that ``len1 >= len2 > 0``.

    int fmpq_poly_divides(fmpq_poly_t q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Return `1` if ``poly2`` divides ``poly1`` and set ``q`` to the quotient.
    # Otherwise return `0`.

    slong fmpq_poly_remove(fmpq_poly_t q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``q`` to the quotient of ``poly1`` by the highest power of ``poly2``
    # which divides it, and returns the power. The divisor ``poly2`` must not be
    # constant or an exception is raised.

    void _fmpq_poly_inv_series_newton(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong n)
    # Computes the first `n` terms of the inverse power series of
    # ``(poly, den, len)`` using Newton iteration.
    # The result is produced in canonical form.
    # Assumes that `n \geq 1` and that ``poly`` has non-zero constant term.
    # Does not support aliasing.

    void fmpq_poly_inv_series_newton(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Computes the first `n` terms of the inverse power series
    # of ``poly`` using Newton iteration, assuming that ``poly``
    # has non-zero constant term and `n \geq 1`.

    void _fmpq_poly_inv_series(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong den_len, slong n)
    # Computes the first `n` terms of the inverse power series of
    # ``(poly, den, len)``.
    # The result is produced in canonical form.
    # Assumes that `n \geq 1` and that ``poly`` has non-zero constant term.
    # Does not support aliasing.

    void fmpq_poly_inv_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Computes the first `n` terms of the inverse power series of ``poly``,
    # assuming that ``poly`` has non-zero constant term and `n \geq 1`.

    void _fmpq_poly_div_series(fmpz * Q, fmpz_t denQ, const fmpz * A, const fmpz_t denA, slong lenA, const fmpz * B, const fmpz_t denB, slong lenB, slong n)
    # Divides ``(A, denA, lenA)`` by ``(B, denB, lenB)`` as power series
    # over `\mathbb{Q}`, assuming `B` has non-zero constant term and that
    # all lengths are positive.
    # Aliasing is not supported.
    # This function ensures that the numerator and denominator
    # are coprime on exit.

    void fmpq_poly_div_series(fmpq_poly_t Q, const fmpq_poly_t A, const fmpq_poly_t B, slong n)
    # Performs power series division in `\mathbb{Q}[[x]] / (x^n)`.  The function
    # considers the polynomials `A` and `B` as power series of length `n`
    # starting with the constant terms.  The function assumes that `B` has
    # non-zero constant term and `n \geq 1`.

    void _fmpq_poly_gcd(fmpz *G, fmpz_t denG, const fmpz *A, slong lenA, const fmpz *B, slong lenB)
    # Computes the monic greatest common divisor `G` of `A` and `B`.
    # Assumes that `G` has space for `\operatorname{len}(B)` coefficients,
    # where `\operatorname{len}(A) \geq \operatorname{len}(B) > 0`.
    # Aliasing between the output and input arguments is not supported.
    # Does not support zero-padding.

    void fmpq_poly_gcd(fmpq_poly_t G, const fmpq_poly_t A, const fmpq_poly_t B)
    # Computes the monic greatest common divisor `G` of `A` and `B`.
    # In the special case when `A = B = 0`, sets `G = 0`.

    void _fmpq_poly_xgcd(fmpz *G, fmpz_t denG, fmpz *S, fmpz_t denS, fmpz *T, fmpz_t denT, const fmpz *A, const fmpz_t denA, slong lenA, const fmpz *B, const fmpz_t denB, slong lenB)
    # Computes polynomials `G`, `S`, and `T` such that
    # `G = \gcd(A, B) = S A + T B`, where `G` is the monic
    # greatest common divisor of `A` and `B`.
    # Assumes that `G`, `S`, and `T` have space for `\operatorname{len}(B)`,
    # `\operatorname{len}(B)`, and `\operatorname{len}(A)` coefficients, respectively,
    # where it is also assumed that `\operatorname{len}(A) \geq \operatorname{len}(B) > 0`.
    # Does not support zero padding of the input arguments.

    void fmpq_poly_xgcd(fmpq_poly_t G, fmpq_poly_t S, fmpq_poly_t T, const fmpq_poly_t A, const fmpq_poly_t B)
    # Computes polynomials `G`, `S`, and `T` such that
    # `G = \gcd(A, B) = S A + T B`, where `G` is the monic
    # greatest common divisor of `A` and `B`.
    # Corner cases are handled as follows.  If `A = B = 0`, returns
    # `G = S = T = 0`.  If `A \neq 0`, `B = 0`, returns the suitable
    # scalar multiple of `G = A`, `S = 1`, and `T = 0`.  The case
    # when `A = 0`, `B \neq 0` is handled similarly.

    void _fmpq_poly_lcm(fmpz *L, fmpz_t denL, const fmpz *A, slong lenA, const fmpz *B, slong lenB)
    # Computes the monic least common multiple `L` of `A` and `B`.
    # Assumes that `L` has space for `\operatorname{len}(A) + \operatorname{len}(B) - 1` coefficients,
    # where `\operatorname{len}(A) \geq \operatorname{len}(B) > 0`.
    # Aliasing between the output and input arguments is not supported.
    # Does not support zero-padding.

    void fmpq_poly_lcm(fmpq_poly_t L, const fmpq_poly_t A, const fmpq_poly_t B)
    # Computes the monic least common multiple `L` of `A` and `B`.
    # In the special case when `A = B = 0`, sets `L = 0`.

    void _fmpq_poly_resultant(fmpz_t rnum, fmpz_t rden, const fmpz *poly1, const fmpz_t den1, slong len1, const fmpz *poly2, const fmpz_t den2, slong len2)
    # Sets ``(rnum, rden)`` to the resultant of the two input
    # polynomials.
    # Assumes that ``len1 >= len2 > 0``.  Does not support zero-padding
    # of the input polynomials.  Does not support aliasing of the input and
    # output arguments.

    void fmpq_poly_resultant(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g)
    # Returns the resultant of `f` and `g`.
    # Enumerating the roots of `f` and `g` over `\bar{\mathbf{Q}}` as
    # `r_1, \dotsc, r_m` and `s_1, \dotsc, s_n`, respectively, and
    # letting `x` and `y` denote the leading coefficients, the resultant
    # is defined as
    # .. math ::
    # x^{\deg(f)} y^{\deg(g)} \prod_{1 \leq i, j \leq n} (r_i - s_j).
    # We handle special cases as follows:  if one of the polynomials is zero,
    # the resultant is zero.  Note that otherwise if one of the polynomials is
    # constant, the last term in the above expression is the empty product.

    void fmpq_poly_resultant_div(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g, const fmpz_t div, slong nbits)
    # Returns the resultant of `f` and `g` divided by ``div`` under the
    # assumption that the result has at most ``nbits`` bits. The result must
    # be an integer.

    void _fmpq_poly_derivative(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)
    # Sets ``(rpoly, rden, len - 1)`` to the derivative of
    # ``(poly, den, len)``.  Does nothing if ``len <= 1``.
    # Supports aliasing between the two polynomials.

    void fmpq_poly_derivative(fmpq_poly_t res, const fmpq_poly_t poly)
    # Sets ``res`` to the derivative of ``poly``.

    void _fmpq_poly_nth_derivative(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, ulong n, slong len)
    # Sets ``(rpoly, rden, len - n)`` to the nth derivative of
    # ``(poly, den, len)``.  Does nothing if ``len <= n``.
    # Supports aliasing between the two polynomials.

    void fmpq_poly_nth_derivative(fmpq_poly_t res, const fmpq_poly_t poly, ulong n)
    # Sets ``res`` to the nth derivative of ``poly``.

    void _fmpq_poly_integral(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)
    # Sets ``(rpoly, rden, len)`` to the integral of
    # ``(poly, den, len - 1)``.  Assumes ``len >= 0``.
    # Supports aliasing between the two polynomials.
    # The output will be in canonical form if the input is
    # in canonical form.

    void fmpq_poly_integral(fmpq_poly_t res, const fmpq_poly_t poly)
    # Sets ``res`` to the integral of ``poly``. The constant
    # term is set to zero. In particular, the integral of the zero
    # polynomial is the zero polynomial.

    void _fmpq_poly_sqrt_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # square root of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 1.
    # Does not support aliasing between the input and output polynomials.

    void fmpq_poly_sqrt_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the square root of ``f``
    # to order ``n > 1``. Requires ``f`` to have constant term 1.

    void _fmpq_poly_invsqrt_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the inverse
    # square root of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 1.
    # Does not support aliasing between the input and output polynomials.

    void fmpq_poly_invsqrt_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the inverse square root of
    # ``f`` to order ``n > 0``. Requires ``f`` to have constant term 1.

    void _fmpq_poly_power_sums(fmpz * res, fmpz_t rden, const fmpz * poly, slong len, slong n)
    # Compute the (truncated) power sums series of the polynomial
    # ``(poly,len)`` up to length `n` using Newton identities.

    void fmpq_poly_power_sums(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Compute the (truncated) power sum series of the monic polynomial
    # ``poly`` up to length `n` using Newton identities. That is the power
    # series whose coefficient of degree `i` is the sum of the `i`-th power of
    # all (complex) roots of the polynomial ``poly``.

    void _fmpq_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, const fmpz_t den, slong len)
    # Compute an integer polynomial given by its power sums series ``(poly,den,len)``.

    void fmpq_poly_power_sums_to_fmpz_poly(fmpz_poly_t res, const fmpq_poly_t Q)
    # Compute the integer polynomial with content one and positive leading
    # coefficient given by its power sums series ``Q``.

    void fmpq_poly_power_sums_to_poly(fmpq_poly_t res, const fmpq_poly_t Q)
    # Compute the monic polynomial from its power sums series ``Q``.

    void _fmpq_poly_log_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # logarithm of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 1.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_log_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the logarithm of ``f``
    # to order ``n > 0``. Requires ``f`` to have constant term 1.

    void _fmpq_poly_exp_series(fmpz * g, fmpz_t gden, const fmpz * h, const fmpz_t hden, slong hlen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # exponential function of ``(h, hden, hlen)``.  Assumes
    # ``n > 0, hlen > 0`` and
    # that ``(h, hden, hlen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t h, slong n)
    # Sets ``res`` to the series expansion of the exponential function
    # of ``h`` to order ``n > 0``. Requires ``f`` to have
    # constant term 0.

    void _fmpq_poly_exp_expinv_series(fmpz * res1, fmpz_t res1den, fmpz * res2, fmpz_t res2den, const fmpz * h, const fmpz_t hden, slong hlen, slong n)
    # The same as ``fmpq_poly_exp_series``, but simultaneously computes
    # the exponential (in ``res1``, ``res1den``) and its multiplicative inverse
    # (in ``res2``, ``res2den``).
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_exp_expinv_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t h, slong n)
    # The same as ``fmpq_poly_exp_series``, but simultaneously computes
    # the exponential (in ``res1``) and its multiplicative inverse
    # (in ``res2``).

    void _fmpq_poly_atan_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # inverse tangent of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_atan_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the inverse tangent of ``f``
    # to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_atanh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the inverse
    # hyperbolic tangent of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_atanh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the inverse hyperbolic
    # tangent of ``f`` to order ``n > 0``. Requires ``f`` to have
    # constant term 0.

    void _fmpq_poly_asin_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # inverse sine of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_asin_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the inverse sine of ``f``
    # to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_asinh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the inverse
    # hyperbolic sine of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_asinh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the inverse hyperbolic
    # sine of ``f`` to order ``n > 0``. Requires ``f`` to have
    # constant term 0.

    void _fmpq_poly_tan_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # tangent function of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Does not support aliasing between the input and output polynomials.

    void fmpq_poly_tan_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the tangent function
    # of ``f`` to order ``n > 0``. Requires ``f`` to have
    # constant term 0.

    void _fmpq_poly_sin_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # sine of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_sin_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the sine of ``f``
    # to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_cos_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # cosine of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the cosine of ``f``
    # to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_sin_cos_series(fmpz * s, fmpz_t sden, fmpz * c, fmpz_t cden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(s, sden, n)`` to the series expansion of the
    # sine of ``(f, fden, flen)``, and ``(c, cden, n)`` to the series
    # expansion of the cosine.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_sin_cos_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t f, slong n)
    # Sets ``res1`` to the series expansion of the sine of ``f``
    # to order ``n > 0``, and ``res2`` to the series expansion
    # of the cosine. Requires ``f`` to have constant term 0.

    void _fmpq_poly_sinh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # hyperbolic sine of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Does not support aliasing between the input and output polynomials.

    void fmpq_poly_sinh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the hyperbolic sine of ``f``
    # to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_cosh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the hyperbolic
    # cosine of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Does not support aliasing between the input and output polynomials.

    void fmpq_poly_cosh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the hyperbolic cosine of
    # ``f`` to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_sinh_cosh_series(fmpz * s, fmpz_t sden, fmpz * c, fmpz_t cden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(s, sden, n)`` to the series expansion of the hyperbolic
    # sine of ``(f, fden, flen)``, and ``(c, cden, n)`` to the series
    # expansion of the hyperbolic cosine.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Supports aliasing between the input and output polynomials.

    void fmpq_poly_sinh_cosh_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t f, slong n)
    # Sets ``res1`` to the series expansion of the hyperbolic sine of ``f``
    # to order ``n > 0``, and ``res2`` to the series expansion
    # of the hyperbolic cosine. Requires ``f`` to have constant term 0.

    void _fmpq_poly_tanh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)
    # Sets ``(g, gden, n)`` to the series expansion of the
    # hyperbolic tangent of ``(f, fden, flen)``.  Assumes ``n > 0`` and
    # that ``(f, fden, flen)`` has constant term 0.
    # Does not support aliasing between the input and output polynomials.

    void fmpq_poly_tanh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    # Sets ``res`` to the series expansion of the hyperbolic tangent of
    # ``f`` to order ``n > 0``. Requires ``f`` to have constant term 0.

    void _fmpq_poly_legendre_p(fmpz * coeffs, fmpz_t den, ulong n)
    # Sets ``coeffs`` to the coefficient array of the Legendre polynomial
    # `P_n(x)`, defined by `(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)`,
    # for `n\ge0`. Sets ``den`` to the overall denominator.
    # The coefficients are calculated using a hypergeometric recurrence.
    # The length of the array will be ``n+1``.
    # To improve performance, the common denominator is computed in one step
    # and the coefficients are evaluated using integer arithmetic. The
    # denominator is given by
    # `\gcd(n!,2^n) = 2^{\lfloor n/2 \rfloor + \lfloor n/4 \rfloor + \ldots}.`
    # See ``fmpz_poly`` for the shifted Legendre polynomials.

    void fmpq_poly_legendre_p(fmpq_poly_t poly, ulong n)
    # Sets ``poly`` to the Legendre polynomial `P_n(x)`, defined
    # by `(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)`, for `n\ge0`.
    # The coefficients are calculated using a hypergeometric recurrence.
    # To improve performance, the common denominator is computed in one step
    # and the coefficients are evaluated using integer arithmetic. The
    # denominator is given by
    # `\gcd(n!,2^n) = 2^{\lfloor n/2 \rfloor + \lfloor n/4 \rfloor + \ldots}.`
    # See ``fmpz_poly`` for the shifted Legendre polynomials.

    void _fmpq_poly_laguerre_l(fmpz * coeffs, fmpz_t den, ulong n)
    # Sets ``coeffs`` to the coefficient array of the Laguerre polynomial
    # `L_n(x)`, defined by `(n+1) L_{n+1}(x) = (2n+1-x) L_n(x) - n L_{n-1}(x)`,
    # for `n\ge0`. Sets ``den`` to the overall denominator.
    # The coefficients are calculated using a hypergeometric recurrence.
    # The length of the array will be ``n+1``.

    void fmpq_poly_laguerre_l(fmpq_poly_t poly, ulong n)
    # Sets ``poly`` to the Laguerre polynomial `L_n(x)`, defined by
    # `(n+1) L_{n+1}(x) = (2n+1-x) L_n(x) - n L_{n-1}(x)`, for `n\ge0`.
    # The coefficients are calculated using a hypergeometric recurrence.

    void _fmpq_poly_gegenbauer_c(fmpz * coeffs, fmpz_t den, ulong n, const fmpq_t a)
    # Sets ``coeffs`` to the coefficient array of the Gegenbauer
    # (ultraspherical) polynomial
    # `C^{(\alpha)}_n(x) = \frac{(2\alpha)_n}{n!}{}_2F_1\left(-n,2\alpha+n;
    # \alpha+\frac12;\frac{1-x}{2}\right)`, for integer `n\ge0` and rational
    # `\alpha>0`. Sets ``den`` to the overall denominator.
    # The coefficients are calculated using a hypergeometric recurrence.

    void fmpq_poly_gegenbauer_c(fmpq_poly_t poly, ulong n, const fmpq_t a)
    # Sets ``poly`` to the Gegenbauer (ultraspherical) polynomial
    # `C^{(\alpha)}_n(x) = \frac{(2\alpha)_n}{n!}{}_2F_1\left(-n,2\alpha+n;
    # \alpha+\frac12;\frac{1-x}{2}\right)`, for integer `n\ge0` and rational
    # `\alpha>0`.
    # The coefficients are calculated using a hypergeometric recurrence.

    void _fmpq_poly_evaluate_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t a)
    # Evaluates the polynomial ``(poly, den, len)`` at the integer `a` and
    # sets ``(rnum, rden)`` to the result in lowest terms.

    void fmpq_poly_evaluate_fmpz(fmpq_t res, const fmpq_poly_t poly, const fmpz_t a)
    # Evaluates the polynomial ``poly`` at the integer `a` and sets
    # ``res`` to the result.

    void _fmpq_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t anum, const fmpz_t aden)
    # Evaluates the polynomial ``(poly, den, len)`` at the rational
    # ``(anum, aden)`` and sets ``(rnum, rden)`` to the result in
    # lowest terms.  Aliasing between ``(rnum, rden)`` and
    # ``(anum, aden)`` is not supported.

    void fmpq_poly_evaluate_fmpq(fmpq_t res, const fmpq_poly_t poly, const fmpq_t a)
    # Evaluates the polynomial ``poly`` at the rational `a` and
    # sets ``res`` to the result.

    void _fmpq_poly_interpolate_fmpz_vec(fmpz * poly, fmpz_t den, const fmpz * xs, const fmpz * ys, slong n)
    # Sets ``poly`` / ``den`` to the unique interpolating polynomial of
    # degree at most `n - 1` satisfying `f(x_i) = y_i` for every pair `x_i, y_i`
    # in ``xs`` and ``ys``.
    # The vector ``poly`` must have room for ``n+1`` coefficients,
    # even if the interpolating polynomial is shorter.
    # Aliasing of ``poly`` or ``den`` with any other argument is not
    # allowed.
    # It is assumed that the `x` values are distinct.
    # This function uses a simple `O(n^2)` implementation of Lagrange
    # interpolation, clearing denominators to avoid working with fractions.
    # It is currently not designed to be efficient for large `n`.

    void fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)
    # Sets ``poly`` to the unique interpolating polynomial of degree
    # at most `n - 1` satisfying `f(x_i) = y_i` for every pair `x_i, y_i`
    # in ``xs`` and ``ys``. It is assumed that the `x` values are distinct.

    void _fmpq_poly_compose(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)
    # Sets ``(res, den)`` to the composition of ``(poly1, den1, len1)``
    # and ``(poly2, den2, len2)``, assuming ``len1, len2 > 0``.
    # Assumes that ``res`` has space for ``(len1 - 1) * (len2 - 1) + 1``
    # coefficients.  Does not support aliasing.

    void fmpq_poly_compose(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    # Sets ``res`` to the composition of ``poly1`` and ``poly2``.

    void _fmpq_poly_rescale(fmpz * res, fmpz_t denr, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t anum, const fmpz_t aden)
    # Sets ``(res, denr, len)`` to ``(poly, den, len)`` with the
    # indeterminate rescaled by ``(anum, aden)``.
    # Assumes that ``len > 0`` and that ``(anum, aden)`` is non-zero and
    # in lowest terms.  Supports aliasing between ``(res, denr, len)`` and
    # ``(poly, den, len)``.

    void fmpq_poly_rescale(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t a)
    # Sets ``res`` to ``poly`` with the indeterminate rescaled by `a`.

    void _fmpq_poly_compose_series_horner(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # Sets ``(res, den, n)`` to the composition of
    # ``(poly1, den1, len1)`` and ``(poly2, den2, len2)`` modulo `x^n`,
    # where the constant term of ``poly2`` is required to be zero.
    # Assumes that ``len1, len2, n > 0``, that ``len1, len2 <= n``,
    # that ``(len1-1) * (len2-1) + 1 <= n``, and that ``res`` has
    # space for ``n`` coefficients. Does not support aliasing between any
    # of the inputs and the output.
    # This implementation uses the Horner scheme.
    # The default ``fmpz_poly`` composition algorithm is automatically
    # used when the composition can be performed over the integers.

    void fmpq_poly_compose_series_horner(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # Sets ``res`` to the composition of ``poly1`` and ``poly2``
    # modulo `x^n`, where the constant term of ``poly2`` is required
    # to be zero.
    # This implementation uses the Horner scheme.
    # The default ``fmpz_poly`` composition algorithm is automatically
    # used when the composition can be performed over the integers.

    void _fmpq_poly_compose_series_brent_kung(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # Sets ``(res, den, n)`` to the composition of
    # ``(poly1, den1, len1)`` and ``(poly2, den2, len2)`` modulo `x^n`,
    # where the constant term of ``poly2`` is required to be zero.
    # Assumes that ``len1, len2, n > 0``, that ``len1, len2 <= n``,
    # that ``(len1-1) * (len2-1) + 1 <= n``, and that ``res`` has
    # space for ``n`` coefficients. Does not support aliasing between any
    # of the inputs and the output.
    # This implementation uses Brent-Kung algorithm 2.1 [BrentKung1978]_.
    # The default ``fmpz_poly`` composition algorithm is automatically
    # used when the composition can be performed over the integers.

    void fmpq_poly_compose_series_brent_kung(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # Sets ``res`` to the composition of ``poly1`` and ``poly2``
    # modulo `x^n`, where the constant term of ``poly2`` is required
    # to be zero.
    # This implementation uses Brent-Kung algorithm 2.1 [BrentKung1978]_.
    # The default ``fmpz_poly`` composition algorithm is automatically
    # used when the composition can be performed over the integers.

    void _fmpq_poly_compose_series(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
    # Sets ``(res, den, n)`` to the composition of
    # ``(poly1, den1, len1)`` and ``(poly2, den2, len2)`` modulo `x^n`,
    # where the constant term of ``poly2`` is required to be zero.
    # Assumes that ``len1, len2, n > 0``, that ``len1, len2 <= n``,
    # that ``(len1-1) * (len2-1) + 1 <= n``, and that ``res`` has
    # space for ``n`` coefficients. Does not support aliasing between any
    # of the inputs and the output.
    # This implementation automatically switches between the Horner scheme
    # and Brent-Kung algorithm 2.1 depending on the size of the inputs.
    # The default ``fmpz_poly`` composition algorithm is automatically
    # used when the composition can be performed over the integers.

    void fmpq_poly_compose_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    # Sets ``res`` to the composition of ``poly1`` and ``poly2``
    # modulo `x^n`, where the constant term of ``poly2`` is required
    # to be zero.
    # This implementation automatically switches between the Horner scheme
    # and Brent-Kung algorithm 2.1 depending on the size of the inputs.
    # The default ``fmpz_poly`` composition algorithm is automatically
    # used when the composition can be performed over the integers.

    void _fmpq_poly_revert_series_lagrange(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
    # Sets ``(res, den)`` to the power series reversion of
    # ``(poly1, den1, len1)`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero. Assumes that `n > 0`.
    # Does not support aliasing between any of the inputs and the output.
    # This implementation uses the Lagrange inversion formula.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void fmpq_poly_revert_series_lagrange(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Sets ``res`` to the power series reversion of ``poly1`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero.
    # This implementation uses the Lagrange inversion formula.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void _fmpq_poly_revert_series_lagrange_fast(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
    # Sets ``(res, den)`` to the power series reversion of
    # ``(poly1, den1, len1)`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero. Assumes that `n > 0`.
    # Does not support aliasing between any of the inputs and the output.
    # This implementation uses a reduced-complexity implementation
    # of the Lagrange inversion formula.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void fmpq_poly_revert_series_lagrange_fast(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Sets ``res`` to the power series reversion of ``poly1`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero.
    # This implementation uses a reduced-complexity implementation
    # of the Lagrange inversion formula.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void _fmpq_poly_revert_series_newton(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
    # Sets ``(res, den)`` to the power series reversion of
    # ``(poly1, den1, len1)`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero. Assumes that `n > 0`.
    # Does not support aliasing between any of the inputs and the output.
    # This implementation uses Newton iteration.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void fmpq_poly_revert_series_newton(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Sets ``res`` to the power series reversion of ``poly1`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero.
    # This implementation uses Newton iteration.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void _fmpq_poly_revert_series(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
    # Sets ``(res, den)`` to the power series reversion of
    # ``(poly1, den1, len1)`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero. Assumes that `n > 0`.
    # Does not support aliasing between any of the inputs and the output.
    # This implementation defaults to using Newton iteration.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void fmpq_poly_revert_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    # Sets ``res`` to the power series reversion of ``poly1`` modulo `x^n`.
    # The constant term of ``poly2`` is required to be zero and
    # the linear term is required to be nonzero.
    # This implementation defaults to using Newton iteration.
    # The default ``fmpz_poly`` reversion algorithm is automatically
    # used when the reversion can be performed over the integers.

    void _fmpq_poly_content(fmpq_t res, const fmpz * poly, const fmpz_t den, slong len)
    # Sets ``res`` to the content of ``(poly, den, len)``.
    # If ``len == 0``, sets ``res`` to zero.

    void fmpq_poly_content(fmpq_t res, const fmpq_poly_t poly)
    # Sets ``res`` to the content of ``poly``.  The content of the zero
    # polynomial is defined to be zero.

    void _fmpq_poly_primitive_part(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)
    # Sets ``(rpoly, rden, len)`` to the primitive part, with non-negative
    # leading coefficient, of ``(poly, den, len)``.  Assumes that
    # ``len > 0``.  Supports aliasing between the two polynomials.

    void fmpq_poly_primitive_part(fmpq_poly_t res, const fmpq_poly_t poly)
    # Sets ``res`` to the primitive part, with non-negative leading
    # coefficient, of ``poly``.

    int _fmpq_poly_is_monic(const fmpz * poly, const fmpz_t den, slong len)
    # Returns whether the polynomial ``(poly, den, len)`` is monic.
    # The zero polynomial is not monic by definition.

    int fmpq_poly_is_monic(const fmpq_poly_t poly)
    # Returns whether the polynomial ``poly`` is monic. The zero
    # polynomial is not monic by definition.

    void _fmpq_poly_make_monic(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)
    # Sets ``(rpoly, rden, len)`` to the monic scalar multiple of
    # ``(poly, den, len)``.  Assumes that ``len > 0``.  Supports
    # aliasing between the two polynomials.

    void fmpq_poly_make_monic(fmpq_poly_t res, const fmpq_poly_t poly)
    # Sets ``res`` to the monic scalar multiple of ``poly`` whenever
    # ``poly`` is non-zero.  If ``poly`` is the zero polynomial, sets
    # ``res`` to zero.

    int fmpq_poly_is_squarefree(const fmpq_poly_t poly)
    # Returns whether the polynomial ``poly`` is square-free.  A non-zero
    # polynomial is defined to be square-free if it has no non-unit square
    # factors.  We also define the zero polynomial to be square-free.

    int _fmpq_poly_print(const fmpz * poly, const fmpz_t den, slong len)
    # Prints the polynomial ``(poly, den, len)`` to ``stdout``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fmpq_poly_print(const fmpq_poly_t poly)
    # Prints the polynomial to ``stdout``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fmpq_poly_print_pretty(const fmpz *poly, const fmpz_t den, slong len, const char * x)

    int fmpq_poly_print_pretty(const fmpq_poly_t poly, const char * var)
    # Prints the pretty representation of ``poly`` to ``stdout``, using
    # the null-terminated string ``var`` not equal to ``"\0"`` as the
    # variable name.
    # In the current implementation always returns `1`.

    int _fmpq_poly_fprint(FILE * file, const fmpz * poly, const fmpz_t den, slong len)
    # Prints the polynomial ``(poly, den, len)`` to the stream ``file``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fmpq_poly_fprint(FILE * file, const fmpq_poly_t poly)
    # Prints the polynomial to the stream ``file``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fmpq_poly_fprint_pretty(FILE * file, const fmpz *poly, const fmpz_t den, slong len, const char * x)

    int fmpq_poly_fprint_pretty(FILE * file, const fmpq_poly_t poly, const char * var)
    # Prints the pretty representation of ``poly`` to ``stdout``, using
    # the null-terminated string ``var`` not equal to ``"\0"`` as the
    # variable name.
    # In the current implementation, always returns `1`.

    int fmpq_poly_read(fmpq_poly_t poly)
    # Reads a polynomial from ``stdin``, storing the result
    # in ``poly``.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive value.

    int fmpq_poly_fread(FILE * file, fmpq_poly_t poly)
    # Reads a polynomial from the stream ``file``, storing the result
    # in ``poly``.
    # In case of success, returns a positive number.  In case of failure,
    # returns a non-positive value.
