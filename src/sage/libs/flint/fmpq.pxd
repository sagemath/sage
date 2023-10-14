# distutils: libraries = flint
# distutils: depends = flint/fmpq.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    fmpz * fmpq_numref(const fmpq_t x)
    fmpz * fmpq_denref(const fmpq_t x)
    # Returns respectively a pointer to the numerator and denominator of x.

    void fmpq_init(fmpq_t x)
    # Initialises the ``fmpq_t`` variable ``x`` for use. Its value
    # is set to 0.

    void fmpq_clear(fmpq_t x)
    # Clears the ``fmpq_t`` variable ``x``. To use the variable again,
    # it must be re-initialised with ``fmpq_init``.

    void fmpq_canonicalise(fmpq_t res)
    # Puts ``res`` in canonical form: the numerator and denominator are
    # reduced to lowest terms, and the denominator is made positive.
    # If the numerator is zero, the denominator is set to one.
    # If the denominator is zero, the outcome of calling this function is
    # undefined, regardless of the value of the numerator.

    void _fmpq_canonicalise(fmpz_t num, fmpz_t den)
    # Does the same thing as ``fmpq_canonicalise``, but for numerator
    # and denominator given explicitly as ``fmpz_t`` variables. Aliasing
    # of ``num`` and ``den`` is not allowed.

    bint fmpq_is_canonical(const fmpq_t x)
    # Returns nonzero if ``fmpq_t`` x is in canonical form
    # (as produced by ``fmpq_canonicalise``), and zero otherwise.

    bint _fmpq_is_canonical(const fmpz_t num, const fmpz_t den)
    # Does the same thing as ``fmpq_is_canonical``, but for numerator
    # and denominator given explicitly as ``fmpz_t`` variables.

    void fmpq_set(fmpq_t dest, const fmpq_t src)
    # Sets ``dest`` to a copy of ``src``. No canonicalisation
    # is performed.

    void fmpq_swap(fmpq_t op1, fmpq_t op2)
    # Swaps the two rational numbers ``op1`` and ``op2``.

    void fmpq_neg(fmpq_t dest, const fmpq_t src)
    # Sets ``dest`` to the additive inverse of ``src``.

    void fmpq_abs(fmpq_t dest, const fmpq_t src)
    # Sets ``dest`` to the absolute value of ``src``.

    void fmpq_zero(fmpq_t res)
    # Sets the value of ``res`` to 0.

    void fmpq_one(fmpq_t res)
    # Sets the value of ``res`` to `1`.

    bint fmpq_is_zero(const fmpq_t res)
    # Returns nonzero if ``res`` has value 0, and returns zero otherwise.

    bint fmpq_is_one(const fmpq_t res)
    # Returns nonzero if ``res`` has value `1`, and returns zero otherwise.

    bint fmpq_is_pm1(const fmpq_t res)
    # Returns nonzero if ``res`` has value `\pm{1}` and zero otherwise.

    bint fmpq_equal(const fmpq_t x, const fmpq_t y)
    # Returns nonzero if ``x`` and ``y`` are equal, and zero otherwise.
    # Assumes that ``x`` and ``y`` are both in canonical form.

    int fmpq_sgn(const fmpq_t x)
    # Returns the sign of the rational number `x`.

    int fmpq_cmp(const fmpq_t x, const fmpq_t y)
    int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
    int fmpq_cmp_ui(const fmpq_t x, ulong y)
    # Returns negative if `x < y`, zero if `x = y`, and positive if `x > y`.

    int fmpq_cmp_si(const fmpq_t x, slong y)
    # Returns negative if `x < y`, zero if `x = y`, and positive if `x > y`.

    bint fmpq_equal_ui(fmpq_t x, ulong y)
    # Returns `1` if `x = y`, otherwise returns `0`.

    bint fmpq_equal_si(fmpq_t x, slong y)
    # Returns `1` if `x = y`, otherwise returns `0`.

    void fmpq_height(fmpz_t height, const fmpq_t x)
    # Sets ``height`` to the height of `x`, defined as the larger of
    # the absolute values of the numerator and denominator of `x`.

    flint_bitcnt_t fmpq_height_bits(const fmpq_t x)
    # Returns the number of bits in the height of `x`.

    void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)
    # Sets ``res`` to the canonical form of the fraction ``p / q``.
    # This is equivalent to assigning the numerator and denominator
    # separately and calling ``fmpq_canonicalise``.

    void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c)
    # Sets ``a``, ``b`` to the numerator and denominator of ``c``
    # respectively.

    void fmpq_set_si(fmpq_t res, slong p, ulong q)
    # Sets ``res`` to the canonical form of the fraction ``p / q``.

    void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q)
    # Sets ``(rnum, rden)`` to the canonical form of the fraction
    # ``p / q``. ``rnum`` and ``rden`` may not be aliased.

    void fmpq_set_ui(fmpq_t res, ulong p, ulong q)
    # Sets ``res`` to the canonical form of the fraction ``p / q``.

    void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q)
    # Sets ``(rnum, rden)`` to the canonical form of the fraction
    # ``p / q``. ``rnum`` and ``rden`` may not be aliased.

    void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
    # Sets the value of ``dest`` to that of the ``mpq_t`` variable
    # ``src``.

    int fmpq_set_str(fmpq_t dest, const char * s, int base)
    # Sets the value of ``dest`` to the value represented in the string
    # ``s`` in base ``base``.
    # Returns 0 if no error occurs. Otherwise returns -1 and ``dest`` is
    # set to zero.

    void fmpq_init_set_mpz_frac_readonly(fmpq_t z, const mpz_t p, const mpz_t q)
    # Assuming ``z`` is an ``fmpz_t`` which will not be cleaned up,
    # this temporarily copies ``p`` and ``q`` into the numerator and
    # denominator of ``z`` for read only operations only. The user must not
    # run ``fmpq_clear`` on ``z``.

    double fmpq_get_d(const fmpq_t f)
    # Returns `f` as a ``double``, rounding towards zero if ``f`` cannot be represented exactly. The return is system dependent if ``f`` is too large or too small to fit in a ``double``.

    void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
    # Sets the value of ``dest``

    int fmpq_get_mpfr(mpfr_t dest, const fmpq_t src, mpfr_rnd_t rnd)
    # Sets the MPFR variable ``dest`` to the value of ``src``,
    # rounded to the nearest representable binary floating-point value
    # in direction ``rnd``. Returns the sign of the rounding,
    # according to MPFR conventions.
    # **Note:** Requires that ``mpfr.h`` has been included before any FLINT
    # header is included.

    char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den)
    char * fmpq_get_str(char * str, int b, const fmpq_t x)
    # Prints the string representation of `x` in base `b \in [2, 36]`
    # to a suitable buffer.
    # If ``str`` is not ``NULL``, this is used as the buffer and
    # also the return value.  If ``str`` is ``NULL``, allocates
    # sufficient space and returns a pointer to the string.

    void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f)
    # Sets the uninitialised ``mpq_t`` `z` to the value of the
    # readonly ``fmpq_t`` `f`.
    # Note that it is assumed that `f` does not change during
    # the lifetime of `z`.
    # The rational `z` has to be cleared by a call to
    # :func:`flint_mpq_clear_readonly`.
    # The suggested use of the two functions is as follows::
    # fmpq_t f;
    # ...
    # {
    # mpq_t z;
    # flint_mpq_init_set_readonly(z, f);
    # foo(..., z);
    # flint_mpq_clear_readonly(z);
    # }
    # This provides a convenient function for user code, only
    # requiring to work with the types ``fmpq_t`` and ``mpq_t``.

    void flint_mpq_clear_readonly(mpq_t z)
    # Clears the readonly ``mpq_t`` `z`.

    void fmpq_init_set_readonly(fmpq_t f, const mpq_t z)
    # Sets the uninitialised ``fmpq_t`` `f` to a readonly
    # version of the rational `z`.
    # Note that the value of `z` is assumed to remain constant
    # throughout the lifetime of `f`.
    # The ``fmpq_t`` `f` has to be cleared by calling the
    # function :func:`fmpq_clear_readonly`.
    # The suggested use of the two functions is as follows::
    # mpq_t z;
    # ...
    # {
    # fmpq_t f;
    # fmpq_init_set_readonly(f, z);
    # foo(..., f);
    # fmpq_clear_readonly(f);
    # }

    void fmpq_clear_readonly(fmpq_t f)
    # Clears the readonly ``fmpq_t`` `f`.

    int fmpq_fprint(FILE * file, const fmpq_t x)
    # Prints ``x`` as a fraction to the stream ``file``.
    # The numerator and denominator are printed verbatim as integers,
    # with a forward slash (/) printed in between.
    # In case of success, returns a positive number. In case of failure,
    # returns a non-positive number.

    int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den)
    # Does the same thing as ``fmpq_fprint``, but for numerator
    # and denominator given explicitly as ``fmpz_t`` variables.
    # In case of success, returns a positive number. In case of failure,
    # returns a non-positive number.

    int fmpq_print(const fmpq_t x)
    # Prints ``x`` as a fraction. The numerator and denominator are
    # printed verbatim as integers, with a forward slash (/) printed in
    # between.
    # In case of success, returns a positive number. In case of failure,
    # returns a non-positive number.

    int _fmpq_print(const fmpz_t num, const fmpz_t den)
    # Does the same thing as ``fmpq_print``, but for numerator
    # and denominator given explicitly as ``fmpz_t`` variables.
    # In case of success, returns a positive number. In case of failure,
    # returns a non-positive number.

    void fmpq_randtest(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
    # Sets ``res`` to a random value, with numerator and denominator
    # having up to ``bits`` bits. The fraction will be in canonical
    # form. This function has an increased probability of generating
    # special values which are likely to trigger corner cases.

    void _fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)
    # Does the same thing as ``fmpq_randtest``, but for numerator
    # and denominator given explicitly as ``fmpz_t`` variables. Aliasing
    # of ``num`` and ``den`` is not allowed.

    void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
    # As per ``fmpq_randtest``, but the result will not be `0`.
    # If ``bits`` is set to `0`, an exception will result.

    void fmpq_randbits(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
    # Sets ``res`` to a random value, with numerator and denominator
    # both having exactly ``bits`` bits before canonicalisation,
    # and then puts ``res`` in canonical form. Note that as a result
    # of the canonicalisation, the resulting numerator and denominator can
    # be slightly smaller than ``bits`` bits.

    void _fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)
    # Does the same thing as ``fmpq_randbits``, but for numerator
    # and denominator given explicitly as ``fmpz_t`` variables. Aliasing
    # of ``num`` and ``den`` is not allowed.

    void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    # Sets ``res`` respectively to ``op1 + op2``, ``op1 - op2``,
    # ``op1 * op2``, or ``op1 / op2``. Assumes that the inputs
    # are in canonical form, and produces output in canonical form.
    # Division by zero results in an error.
    # Aliasing between any combination of the variables is allowed.

    void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    # Sets ``(rnum, rden)`` to the canonical form of the sum,
    # difference, product or quotient respectively of the fractions
    # represented by ``(op1num, op1den)`` and ``(op2num, op2den)``.
    # Aliasing between any combination of the variables is allowed,
    # whilst no numerator is aliased with a denominator.

    void _fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
    void _fmpq_sub_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
    void _fmpq_add_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
    void _fmpq_sub_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
    void _fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)
    void _fmpq_sub_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)
    # Sets ``(rnum, rden)`` to the canonical form of the sum or difference
    # respectively of the fractions represented by ``(p, q)`` and
    # ``(r, 1)``. Numerators may not be aliased with denominators.

    void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c)
    void fmpq_sub_si(fmpq_t res, const fmpq_t op1, slong c)
    void fmpq_add_ui(fmpq_t res, const fmpq_t op1, ulong c)
    void fmpq_sub_ui(fmpq_t res, const fmpq_t op1, ulong c)
    void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
    void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)

    void _fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)

    void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c)

    void _fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)

    void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c)

    void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    # Sets ``res`` to ``res + op1 * op2`` or ``res - op1 * op2``
    # respectively, placing the result in canonical form. Aliasing
    # between any combination of the variables is allowed.

    void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    # Sets ``(rnum, rden)`` to the canonical form of the fraction
    # ``(rnum, rden)`` + ``(op1num, op1den)`` * ``(op2num, op2den)`` or
    # ``(rnum, rden)`` - ``(op1num, op1den)`` * ``(op2num, op2den)``
    # respectively. Aliasing between any combination of the variables is allowed,
    # whilst no numerator is aliased with a denominator.

    void fmpq_inv(fmpq_t dest, const fmpq_t src)
    # Sets ``dest`` to ``1 / src``. The result is placed in canonical
    # form, assuming that ``src`` is already in canonical form.

    void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, const fmpz_t opnum, const fmpz_t opden, slong e)
    void fmpq_pow_si(fmpq_t res, const fmpq_t op, slong e)
    # Sets ``res`` to ``op`` raised to the power `e`, where `e`
    # is a ``slong``.  If `e` is `0` and ``op`` is `0`, then
    # ``res`` will be set to `1`.

    int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)
    # Set ``res`` to ``op`` raised to the power `e`.
    # Return `1` for success and `0` for failure.

    void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
    # Sets ``res`` to the product of the rational number ``op``
    # and the integer ``x``.

    void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
    # Sets ``res`` to the quotient of the rational number ``op``
    # and the integer ``x``.

    void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp)
    # Sets ``res`` to ``x`` multiplied by ``2^exp``.

    void fmpq_div_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp)
    # Sets ``res`` to ``x`` divided by ``2^exp``.

    void _fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)
    # Set ``(rnum, rden)`` to the gcd of ``(p, q)`` and ``(r, s)``
    # which we define to be the canonicalisation of `\operatorname{gcd}(ps, qr)/(qs)`.
    # (This is apparently Euclid's original definition and is stable under scaling of
    # numerator and denominator. It also agrees with the gcd on the integers.
    # Note that it does not agree with gcd as defined in ``fmpq_poly``.)
    # This definition agrees with the result as output by Sage and Pari/GP.

    void fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    # Set ``res`` to the gcd of ``op1`` and ``op2``. See the low
    # level function ``_fmpq_gcd`` for our definition of gcd.

    void _fmpq_gcd_cofactors(fmpz_t gnum, fmpz_t gden, fmpz_t abar, fmpz_t bbar, const fmpz_t anum, const fmpz_t aden, const fmpz_t bnum, const fmpz_t bden)
    void fmpq_gcd_cofactors(fmpq_t g, fmpz_t abar, fmpz_t bbar, const fmpq_t a, const fmpq_t b)
    # Set `g` to `\operatorname{gcd}(a,b)` as per :func:`fmpq_gcd` and also compute `\overline{a} = a/g` and `\overline{b} = b/g`.
    # Unlike :func:`fmpq_gcd`, this function requires canonical inputs.

    void _fmpq_add_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)
    # Sets ``(rnum, rden)`` to the sum of ``(p1, q1)`` and ``(p2, q2)``.
    # Assumes that ``(p1, q1)`` and ``(p2, q2)`` are in canonical form
    # and that all inputs are between ``COEFF_MIN`` and ``COEFF_MAX``.

    void _fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)
    # Sets ``(rnum, rden)`` to the product of ``(p1, q1)`` and ``(p2, q2)``.
    # Assumes that ``(p1, q1)`` and ``(p2, q2)`` are in canonical form
    # and that all inputs are between ``COEFF_MIN`` and ``COEFF_MAX``.

    int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod)
    int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)
    # Sets the integer ``res`` to the residue `a` of
    # `x = n/d` = ``(num, den)`` modulo the positive integer `m` = ``mod``,
    # defined as the `0 \le a < m` satisfying `n \equiv a d \pmod m`.
    # If such an `a` exists, 1 will be returned, otherwise 0 will
    # be returned.

    int _fmpq_reconstruct_fmpz_2_naive(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    # Reconstructs a rational number from its residue `a` modulo `m`.
    # Given a modulus `m > 2`, a residue `0 \le a < m`, and positive `N, D`
    # satisfying `2ND < m`, this function attempts to find a fraction `n/d` with
    # `0 \le |n| \le N` and `0 < d \le D` such that `\gcd(n,d) = 1` and
    # `n \equiv ad \pmod m`. If a solution exists, then it is also unique.
    # The function returns 1 if successful, and 0 to indicate that no solution
    # exists.

    int _fmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m)
    int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
    # Reconstructs a rational number from its residue `a` modulo `m`,
    # returning 1 if successful and 0 if no solution exists.
    # Uses the balanced bounds `N = D = \lfloor\sqrt{\frac{m-1}{2}}\rfloor`.

    void _fmpq_next_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_minimal(fmpq_t res, const fmpq_t x)
    # Given `x` which is assumed to be nonnegative and in canonical form, sets
    # ``res`` to the next rational number in the sequence obtained by
    # enumerating all positive denominators `q`, for each `q` enumerating
    # the numerators `1 \le p < q` in order and generating both `p/q` and `q/p`,
    # but skipping all `\gcd(p,q) \ne 1`. Starting with zero, this generates
    # every nonnegative rational number once and only once, with the first
    # few entries being:
    # `0, 1, 1/2, 2, 1/3, 3, 2/3, 3/2, 1/4, 4, 3/4, 4/3, 1/5, 5, 2/5, \ldots.`
    # This enumeration produces the rational numbers in order of
    # minimal height. It has the disadvantage of being somewhat slower to
    # compute than the Calkin-Wilf enumeration.

    void _fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x)
    # Given a signed rational number `x` assumed to be in canonical form, sets
    # ``res`` to the next element in the minimal-height sequence
    # generated by ``fmpq_next_minimal`` but with negative numbers
    # interleaved:
    # `0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots.`
    # Starting with zero, this generates every rational number once
    # and only once, in order of minimal height.

    void _fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)
    # Given `x` which is assumed to be nonnegative and in canonical form, sets
    # ``res`` to the next number in the breadth-first traversal of the
    # Calkin-Wilf tree. Starting with zero, this generates every nonnegative
    # rational number once and only once, with the first few entries being:
    # `0, 1, 1/2, 2, 1/3, 3/2, 2/3, 3, 1/4, 4/3, 3/5, 5/2, 2/5, \ldots.`
    # Despite the appearance of the initial entries, the Calkin-Wilf
    # enumeration does not produce the rational numbers in order of height:
    # some small fractions will appear late in the sequence. This order
    # has the advantage of being faster to produce than the minimal-height
    # order.

    void _fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x)
    # Given a signed rational number `x` assumed to be in canonical form, sets
    # ``res`` to the next element in the Calkin-Wilf sequence with
    # negative numbers interleaved:
    # `0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots.`
    # Starting with zero, this generates every rational number once
    # and only once, but not in order of minimal height.

    void fmpq_farey_neighbors(fmpq_t l, fmpq_t r, const fmpq_t x, const fmpz_t Q)
    # Set `l` and `r` to the fractions directly below and above `x` in the Farey sequence of order `Q`.
    # This function will throw if `x` is not canonical or `Q` is less than the denominator of `x`.

    void fmpq_simplest_between(fmpq_t x, const fmpq_t l, const fmpq_t r)
    void _fmpq_simplest_between(fmpz_t x_num, fmpz_t x_den, const fmpz_t l_num, const fmpz_t l_den, const fmpz_t r_num, const fmpz_t r_den)
    # Set `x` to the simplest fraction in the closed interval `[l, r]`. The underscore version makes the additional assumption that `l \le r`.
    # The endpoints `l` and `r` do not need to be reduced, but their denominators do need to be positive.
    # `x` will always be returned in canonical form. A canonical fraction `a_1/b_1` is defined to be simpler than `a_2/b_2` iff `b_1<b_2` or `b_1=b_2` and `a_1<a_2`.

    slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)
    slong fmpq_get_cfrac_naive(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)
    # Generates up to `n` terms of the (simple) continued fraction expansion
    # of `x`, writing the coefficients to the vector `c` and the remainder `r`
    # to the ``rem`` variable. The return value is the number `k` of
    # generated terms. The output satisfies
    # .. math ::
    # x = c_0 + \cfrac{1}{c_1 + \cfrac{1}{c_2 +
    # \cfrac{1}{ \ddots + \cfrac{1}{c_{k-1} + r }}}}
    # If `r` is zero, the continued fraction expansion is complete.
    # If `r` is nonzero, `1/r` can be passed back as input to generate
    # `c_k, c_{k+1}, \ldots`. Calls to ``fmpq_get_cfrac`` can therefore
    # be chained to generate the continued fraction incrementally,
    # extracting any desired number of coefficients at a time.
    # In general, a rational number has exactly two continued fraction
    # expansions. By convention, we generate the shorter one. The longer
    # expansion can be obtained by replacing the last coefficient
    # `a_{k-1}` by the pair of coefficients `a_{k-1} - 1, 1`.
    # The behaviour of this function in corner cases is as follows:
    # - if `x` is infinite (anything over 0), ``rem`` will be zero and the return is `k=0` regardless of `n`.
    # - else (if `x` is finite),
    # - if `n \le 0`, ``rem`` will be `1/x` (allowing for infinite in the case `x=0`) and the return is `k=0`
    # - else (if `n > 0`), ``rem`` will finite and the return is `0 < k \le n`.
    # Essentially, if this function is called with canonical `x` and `n > 0`, then ``rem`` will be canonical.
    # Therefore, applications relying on canonical ``fmpq_t``'s should not call this function with `n \le 0`.

    void fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n)
    # Sets `x` to the value of the continued fraction
    # .. math ::
    # x = c_0 + \cfrac{1}{c_1 + \cfrac{1}{c_2 +
    # \cfrac{1}{ \ddots + \cfrac{1}{c_{n-1}}}}}
    # where all `c_i` except `c_0` should be nonnegative.
    # It is assumed that `n > 0`.
    # For large `n`, this function implements a subquadratic algorithm.
    # The convergents are given by a chain product of 2 by 2 matrices.
    # This product is split in half recursively to balance the size
    # of the coefficients.

    slong fmpq_cfrac_bound(const fmpq_t x)
    # Returns an upper bound for the number of terms in the continued
    # fraction expansion of `x`. The computed bound is not necessarily sharp.
    # We use the fact that the smallest denominator
    # that can give a continued fraction of length `n` is the Fibonacci
    # number `F_{n+1}`.

    void _fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n)
    void fmpq_harmonic_ui(fmpq_t x, ulong n)
    # Computes the harmonic number `H_n = 1 + 1/2 + 1/3 + \dotsb + 1/n`.
    # Table lookup is used for `H_n` whose numerator and denominator
    # fit in single limb. For larger `n`, a divide and conquer strategy is used.

    void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
    void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)
    # Computes `s(h,k)` for arbitrary `h` and `k`. The naive version uses a straightforward
    # implementation of the defining sum using ``fmpz`` arithmetic and is slow for large `k`.
