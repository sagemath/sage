# distutils: libraries = flint
# distutils: depends = flint/fmpz_vec.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    fmpz * _fmpz_vec_init(slong len)
    # Returns an initialised vector of ``fmpz``'s of given length.

    void _fmpz_vec_clear(fmpz * vec, slong len)
    # Clears the entries of ``(vec, len)`` and frees the space allocated
    # for ``vec``.

    void _fmpz_vec_randtest(fmpz * f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # Sets the entries of a vector of the given length to random integers with
    # up to the given number of bits per entry.

    void _fmpz_vec_randtest_unsigned(fmpz * f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    # Sets the entries of a vector of the given length to random unsigned
    # integers with up to the given number of bits per entry.

    slong _fmpz_vec_max_bits(const fmpz * vec, slong len)
    # If `b` is the maximum number of bits of the absolute value of any
    # coefficient of ``vec``, then if any coefficient of ``vec`` is
    # negative, `-b` is returned, else `b` is returned.

    slong _fmpz_vec_max_bits_ref(const fmpz * vec, slong len)
    # If `b` is the maximum number of bits of the absolute value of any
    # coefficient of ``vec``, then if any coefficient of ``vec`` is
    # negative, `-b` is returned, else `b` is returned.
    # This is a slower reference implementation of ``_fmpz_vec_max_bits``.

    void _fmpz_vec_sum_max_bits(slong * sumabs, slong * maxabs, const fmpz * vec, slong len)
    # Sets ``sumabs`` to the bit count of the sum of the absolute values of
    # the elements of ``vec``. Sets ``maxabs`` to the bit count of the
    # maximum of the absolute values of the elements of ``vec``.

    mp_size_t _fmpz_vec_max_limbs(const fmpz * vec, slong len)
    # Returns the maximum number of limbs needed to store the absolute value
    # of any entry in ``(vec, len)``.  If all entries are zero, returns
    # zero.

    void _fmpz_vec_height(fmpz_t height, const fmpz * vec, slong len)
    # Computes the height of ``(vec, len)``, defined as the largest of the
    # absolute values the coefficients. Equivalently, this gives the infinity
    # norm of the vector. If ``len`` is zero, the height is `0`.

    slong _fmpz_vec_height_index(const fmpz * vec, slong len)
    # Returns the index of an entry of maximum absolute value in the vector.
    # The length must be at least 1.

    int _fmpz_vec_fread(FILE * file, fmpz ** vec, slong * len)
    # Reads a vector from the stream ``file`` and stores it at
    # ``*vec``.  The format is the same as the output format of
    # ``_fmpz_vec_fprint()``, followed by either any character
    # or the end of the file.
    # The interpretation of the various input arguments depends on whether
    # or not ``*vec`` is ``NULL``:
    # If ``*vec == NULL``, the value of ``*len`` on input is ignored.
    # Once the length has been read from ``file``, ``*len`` is set
    # to that value and a vector of this length is allocated at ``*vec``.
    # Finally, ``*len`` coefficients are read from the input stream.  In
    # case of a file or parsing error, clears the vector and sets ``*vec``
    # and ``*len`` to ``NULL`` and ``0``, respectively.
    # Otherwise, if ``*vec != NULL``, it is assumed that ``(*vec, *len)``
    # is a properly initialised vector.  If the length on the input stream
    # does not match ``*len``, a parsing error is raised.  Attempts to read
    # the right number of coefficients from the input stream.  In case of a
    # file or parsing error, leaves the vector ``(*vec, *len)`` in its
    # current state.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fmpz_vec_read(fmpz ** vec, slong * len)
    # Reads a vector from ``stdin`` and stores it at ``*vec``.
    # For further details, see ``_fmpz_vec_fread()``.

    int _fmpz_vec_fprint(FILE * file, const fmpz * vec, slong len)
    # Prints the vector of given length to the stream ``file``. The
    # format is the length followed by two spaces, then a space separated
    # list of coefficients. If the length is zero, only `0` is printed.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fmpz_vec_print(const fmpz * vec, slong len)
    # Prints the vector of given length to ``stdout``.
    # For further details, see ``_fmpz_vec_fprint()``.

    void _fmpz_vec_get_nmod_vec(mp_ptr res, const fmpz * poly, slong len, nmod_t mod)
    # Reduce the coefficients of ``(poly, len)`` modulo the given
    # modulus and set ``(res, len)`` to the result.

    void _fmpz_vec_set_nmod_vec(fmpz * res, mp_srcptr poly, slong len, nmod_t mod)
    # Set the coefficients of ``(res, len)`` to the symmetric modulus
    # of the coefficients of ``(poly, len)``, i.e. convert the given
    # coefficients modulo the given modulus `n` to their signed integer
    # representatives in the range `[-n/2, n/2)`.

    void _fmpz_vec_get_fft(mp_limb_t ** coeffs_f, const fmpz * coeffs_m, slong l, slong length)
    # Convert the vector of coeffs ``coeffs_m`` to an fft vector
    # ``coeffs_f`` of the given ``length`` with ``l`` limbs per
    # coefficient with an additional limb for overflow.

    void _fmpz_vec_set_fft(fmpz * coeffs_m, slong length, const mp_ptr * coeffs_f, slong limbs, slong sign)
    # Convert an fft vector ``coeffs_f`` of fully reduced Fermat numbers of the
    # given ``length`` to a vector of ``fmpz``'s. Each is assumed to be the given
    # number of limbs in length with an additional limb for overflow. If the
    # output coefficients are to be signed then set ``sign``, otherwise clear it.
    # The resulting ``fmpz``s will be in the range `[-n,n]` in the signed case
    # and in the range `[0,2n]` in the unsigned case where
    # ``n = 2^(FLINT_BITS*limbs - 1)``.

    slong _fmpz_vec_get_d_vec_2exp(double * appv, const fmpz * vec, slong len)
    # Export the array of ``len`` entries starting at the pointer ``vec``
    # to an array of doubles ``appv``, each entry of which is notionally
    # multiplied by a single returned exponent to give the original entry. The
    # returned exponent is set to be the maximum exponent of all the original
    # entries so that all the doubles in ``appv`` have a maximum absolute
    # value of 1.0.

    void _fmpz_vec_set(fmpz * vec1, const fmpz * vec2, slong len2)
    # Makes a copy of ``(vec2, len2)`` into ``vec1``.

    void _fmpz_vec_swap(fmpz * vec1, fmpz * vec2, slong len2)
    # Swaps the integers in ``(vec1, len2)`` and ``(vec2, len2)``.

    void _fmpz_vec_zero(fmpz * vec, slong len)
    # Zeros the entries of ``(vec, len)``.

    void _fmpz_vec_neg(fmpz * vec1, const fmpz * vec2, slong len2)
    # Negates ``(vec2, len2)`` and places it into ``vec1``.

    void _fmpz_vec_scalar_abs(fmpz * vec1, const fmpz * vec2, slong len2)
    # Takes the absolute value of entries in ``(vec2, len2)`` and places the
    # result into ``vec1``.

    int _fmpz_vec_equal(const fmpz * vec1, const fmpz * vec2, slong len)
    # Compares two vectors of the given length and returns `1` if they are
    # equal, otherwise returns `0`.

    int _fmpz_vec_is_zero(const fmpz * vec, slong len)
    # Returns `1` if ``(vec, len)`` is zero, and `0` otherwise.

    void _fmpz_vec_max(fmpz * vec1, const fmpz * vec2, const fmpz * vec3, slong len)
    # Sets ``vec1`` to the pointwise maximum of ``vec2`` and ``vec3``.

    void _fmpz_vec_max_inplace(fmpz * vec1, const fmpz * vec2, slong len)
    # Sets ``vec1`` to the pointwise maximum of ``vec1`` and ``vec2``.

    void _fmpz_vec_sort(fmpz * vec, slong len)
    # Sorts the coefficients of ``vec`` in ascending order.

    void _fmpz_vec_add(fmpz * res, const fmpz * vec1, const fmpz * vec2, slong len2)
    # Sets ``(res, len2)`` to the sum of ``(vec1, len2)``
    # and ``(vec2, len2)``.

    void _fmpz_vec_sub(fmpz * res, const fmpz * vec1, const fmpz * vec2, slong len2)
    # Sets ``(res, len2)`` to ``(vec1, len2)`` minus ``(vec2, len2)``.

    void _fmpz_vec_scalar_mul_fmpz(fmpz * vec1, const fmpz * vec2, slong len2, const fmpz_t x)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` multiplied by `c`,
    # where `c` is an ``fmpz_t``.

    void _fmpz_vec_scalar_mul_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` multiplied by `c`,
    # where `c` is a ``slong``.

    void _fmpz_vec_scalar_mul_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` multiplied by `c`,
    # where `c` is an ``ulong``.

    void _fmpz_vec_scalar_mul_2exp(fmpz * vec1, const fmpz * vec2, slong len2, ulong exp)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` multiplied by ``2^exp``.

    void _fmpz_vec_scalar_divexact_fmpz(fmpz * vec1, const fmpz * vec2, slong len2, const fmpz_t x)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `x`, where the
    # division is assumed to be exact for every entry in ``vec2``.

    void _fmpz_vec_scalar_divexact_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `x`, where the
    # division is assumed to be exact for every entry in ``vec2``.

    void _fmpz_vec_scalar_divexact_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `x`, where the
    # division is assumed to be exact for every entry in ``vec2``.

    void _fmpz_vec_scalar_fdiv_q_fmpz(fmpz * vec1, const fmpz * vec2, slong len2, const fmpz_t c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `c`, rounding
    # down towards minus infinity whenever the division is not exact.

    void _fmpz_vec_scalar_fdiv_q_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `c`, rounding
    # down towards minus infinity whenever the division is not exact.

    void _fmpz_vec_scalar_fdiv_q_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `c`, rounding
    # down towards minus infinity whenever the division is not exact.

    void _fmpz_vec_scalar_fdiv_q_2exp(fmpz * vec1, const fmpz * vec2, slong len2, ulong exp)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by ``2^exp``,
    # rounding down towards minus infinity whenever the division is not exact.

    void _fmpz_vec_scalar_fdiv_r_2exp(fmpz * vec1, const fmpz * vec2, slong len2, ulong exp)
    # Sets ``(vec1, len2)`` to the remainder of ``(vec2, len2)``
    # divided by ``2^exp``, rounding down the quotient towards minus
    # infinity whenever the division is not exact.

    void _fmpz_vec_scalar_tdiv_q_fmpz(fmpz * vec1, const fmpz * vec2, slong len2, const fmpz_t c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `c`, rounding
    # towards zero whenever the division is not exact.

    void _fmpz_vec_scalar_tdiv_q_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `c`, rounding
    # towards zero whenever the division is not exact.

    void _fmpz_vec_scalar_tdiv_q_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by `c`, rounding
    # towards zero whenever the division is not exact.

    void _fmpz_vec_scalar_tdiv_q_2exp(fmpz * vec1, const fmpz * vec2, slong len2, ulong exp)
    # Sets ``(vec1, len2)`` to ``(vec2, len2)`` divided by ``2^exp``,
    # rounding down towards zero whenever the division is not exact.

    void _fmpz_vec_scalar_addmul_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)

    void _fmpz_vec_scalar_addmul_ui(fmpz * vec1, const fmpz * vec2, slong len2, ulong c)

    void _fmpz_vec_scalar_addmul_fmpz(fmpz * vec1, const fmpz * vec2, slong len2, const fmpz_t c)
    # Adds ``(vec2, len2)`` times `c` to ``(vec1, len2)``.

    void _fmpz_vec_scalar_addmul_si_2exp(fmpz * vec1, const fmpz * vec2, slong len2, slong c, ulong exp)
    # Adds ``(vec2, len2)`` times ``c * 2^exp`` to ``(vec1, len2)``,
    # where `c` is a ``slong``.

    void _fmpz_vec_scalar_submul_fmpz(fmpz * vec1, const fmpz * vec2, slong len2, const fmpz_t x)
    # Subtracts ``(vec2, len2)`` times `c` from ``(vec1, len2)``,
    # where `c` is a ``fmpz_t``.

    void _fmpz_vec_scalar_submul_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
    # Subtracts ``(vec2, len2)`` times `c` from ``(vec1, len2)``,
    # where `c` is a ``slong``.

    void _fmpz_vec_scalar_submul_si_2exp(fmpz * vec1, const fmpz * vec2, slong len2, slong c, ulong e)
    # Subtracts ``(vec2, len2)`` times `c \times 2^e`
    # from ``(vec1, len2)``, where `c` is a ``slong``.

    void _fmpz_vec_sum(fmpz_t res, const fmpz * vec, slong len)
    # Sets ``res`` to the sum of the entries in ``(vec, len)``.
    # Aliasing of ``res`` with the entries in ``vec`` is not permitted.

    void _fmpz_vec_prod(fmpz_t res, const fmpz * vec, slong len)
    # Sets ``res`` to the product of the entries in ``(vec, len)``.
    # Aliasing of ``res`` with the entries in ``vec`` is not permitted.
    # Uses binary splitting.

    void _fmpz_vec_scalar_mod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
    # Reduces all entries in ``(vec, len)`` modulo `p > 0`.

    void _fmpz_vec_scalar_smod_fmpz(fmpz *res, const fmpz *vec, slong len, const fmpz_t p)
    # Reduces all entries in ``(vec, len)`` modulo `p > 0`, choosing
    # the unique representative in `(-p/2, p/2]`.

    void _fmpz_vec_content(fmpz_t res, const fmpz * vec, slong len)
    # Sets ``res`` to the non-negative content of the entries in ``vec``.
    # The content of a zero vector, including the case when the length is zero,
    # is defined to be zero.

    void _fmpz_vec_content_chained(fmpz_t res, const fmpz * vec, slong len, const fmpz_t input)
    # Sets ``res`` to the non-negative content of ``input`` and the entries in ``vec``.
    # This is useful for calculating the common content of several vectors.

    void _fmpz_vec_lcm(fmpz_t res, const fmpz * vec, slong len)
    # Sets ``res`` to the nonnegative least common multiple of the entries
    # in ``vec``. The least common multiple is zero if any entry in
    # the vector is zero. The least common multiple of a length zero vector is
    # defined to be one.

    void _fmpz_vec_dot(fmpz_t res, const fmpz * vec1, const fmpz * vec2, slong len2)
    # Sets ``res`` to the dot product of ``(vec1, len2)`` and
    # ``(vec2, len2)``.

    void _fmpz_vec_dot_ptr(fmpz_t res, const fmpz * vec1, fmpz ** const vec2, slong offset, slong len)
    # Sets ``res`` to the dot product of ``len`` values at ``vec1`` and the
    # ``len`` values ``vec2[i] + offset`` for `0 \leq i < len`.
