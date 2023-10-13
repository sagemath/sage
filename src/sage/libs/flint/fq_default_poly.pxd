# distutils: libraries = flint
# distutils: depends = flint/fq_default_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fq_default_poly_init(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Initialises ``poly`` for use, with context ctx, and setting its
    # length to zero. A corresponding call to :func:`fq_default_poly_clear`
    # must be made after finishing with the ``fq_default_poly_t`` to free the
    # memory used by the polynomial.

    void fq_default_poly_init2(fq_default_poly_t poly, long alloc, const fq_default_ctx_t ctx)
    # Initialises ``poly`` with space for at least ``alloc``
    # coefficients and sets the length to zero.  The allocated
    # coefficients are all set to zero.  A corresponding call to
    # :func:`fq_default_poly_clear` must be made after finishing with the
    # ``fq_default_poly_t`` to free the memory used by the polynomial.

    void fq_default_poly_realloc(fq_default_poly_t poly, long alloc, const fq_default_ctx_t ctx)
    # Reallocates the given polynomial to have space for ``alloc``
    # coefficients.  If ``alloc`` is zero the polynomial is cleared
    # and then reinitialised.  If the current length is greater than
    # ``alloc`` the polynomial is first truncated to length
    # ``alloc``.

    void fq_default_poly_fit_length(fq_default_poly_t poly, long len, const fq_default_ctx_t ctx)
    # If ``len`` is greater than the number of coefficients currently
    # allocated, then the polynomial is reallocated to have space for at
    # least ``len`` coefficients.  No data is lost when calling this
    # function.
    # The function efficiently deals with the case where
    # ``fit_length`` is called many times in small increments by at
    # least doubling the number of allocated coefficients when length is
    # larger than the number of coefficients currently allocated.

    void fq_default_poly_clear(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Clears the given polynomial, releasing any memory used.  It must
    # be reinitialised in order to be used again.

    void _fq_default_poly_set_length(fq_default_poly_t poly, long len, const fq_default_ctx_t ctx)
    # Set the length of ``poly`` to ``len``.

    void fq_default_poly_truncate(fq_default_poly_t poly, long newlen, const fq_default_ctx_t ctx)
    # Truncates the polynomial to length at most `n`.

    void fq_default_poly_set_trunc(fq_default_poly_t poly1, fq_default_poly_t poly2, long newlen, const fq_default_ctx_t ctx)
    # Sets ``poly1`` to ``poly2`` truncated to length `n`.

    void fq_default_poly_reverse(fq_default_poly_t output, const fq_default_poly_t input, long m, const fq_default_ctx_t ctx)
    # Sets ``output`` to the reverse of ``input``, thinking of it
    # as a polynomial of length ``m``, notionally zero-padded if
    # necessary).  The length ``m`` must be non-negative, but there
    # are no other restrictions. The output polynomial will be set to
    # length ``m`` and then normalised.

    long fq_default_poly_degree(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Returns the degree of the polynomial ``poly``.

    long fq_default_poly_length(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Returns the length of the polynomial ``poly``.

    void fq_default_poly_randtest(fq_default_poly_t f, flint_rand_t state, long len, const fq_default_ctx_t ctx)
    # Sets `f` to a random polynomial of length at most ``len``
    # with entries in the field described by ``ctx``.

    void fq_default_poly_randtest_not_zero(fq_default_poly_t f, flint_rand_t state, long len, const fq_default_ctx_t ctx)
    # Same as ``fq_default_poly_randtest`` but guarantees that the polynomial
    # is not zero.

    void fq_default_poly_randtest_monic(fq_default_poly_t f, flint_rand_t state, long len, const fq_default_ctx_t ctx)
    # Sets `f` to a random monic polynomial of length ``len`` with
    # entries in the field described by ``ctx``.

    void fq_default_poly_randtest_irreducible(fq_default_poly_t f, flint_rand_t state, long len, const fq_default_ctx_t ctx)
    # Sets `f` to a random monic, irreducible polynomial of length
    # ``len`` with entries in the field described by ``ctx``.

    void fq_default_poly_set(fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    # Sets the polynomial ``poly1`` to the polynomial ``poly2``.

    void fq_default_poly_set_fq_default(fq_default_poly_t poly, const fq_default_t c, const fq_default_ctx_t ctx)
    # Sets the polynomial ``poly`` to ``c``.

    void fq_default_poly_swap(fq_default_poly_t op1, fq_default_poly_t op2, const fq_default_ctx_t ctx)
    # Swaps the two polynomials ``op1`` and ``op2``.

    void fq_default_poly_zero(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Sets ``poly`` to the zero polynomial.

    void fq_default_poly_one(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Sets ``poly`` to the constant polynomial `1`.

    void fq_default_poly_gen(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Sets ``poly`` to the polynomial `x`.

    void fq_default_poly_make_monic(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Sets ``rop`` to ``op``, normed to have leading coefficient 1.

    void fq_default_poly_set_nmod_poly(fq_default_poly_t rop, const nmod_poly_t op, const fq_default_ctx_t ctx)
    # Sets the polynomial ``rop`` to the polynomial ``op``.

    void fq_default_poly_set_fmpz_mod_poly(fq_default_poly_t rop, const fmpz_mod_poly_t op, const fq_default_ctx_t ctx)
    # Sets the polynomial ``rop`` to the polynomial ``op``.

    void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop, const fmpz_poly_t op, const fq_default_ctx_t ctx)
    # Sets the polynomial ``rop`` to the polynomial ``op``.

    void fq_default_poly_get_coeff(fq_default_t x, const fq_default_poly_t poly, long n, const fq_default_ctx_t ctx)
    # Sets `x` to the coefficient of `X^n` in ``poly``.

    void fq_default_poly_set_coeff(fq_default_poly_t poly, long n, const fq_default_t x, const fq_default_ctx_t ctx)
    # Sets the coefficient of `X^n` in ``poly`` to `x`.

    void fq_default_poly_set_coeff_fmpz(fq_default_poly_t poly, long n, const fmpz_t x, const fq_default_ctx_t ctx)
    # Sets the coefficient of `X^n` in the polynomial to `x`,
    # assuming `n \geq 0`.

    int fq_default_poly_equal(const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    # Returns nonzero if the two polynomials ``poly1`` and ``poly2``
    # are equal, otherwise returns zero.

    int fq_default_poly_equal_trunc(const fq_default_poly_t poly1, const fq_default_poly_t poly2, long n, const fq_default_ctx_t ctx)
    # Notionally truncate ``poly1`` and ``poly2`` to length `n` and
    # return nonzero if they are equal, otherwise return zero.

    int fq_default_poly_is_zero(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is the zero polynomial.

    int fq_default_poly_is_one(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is equal
    # to the constant polynomial `1`.

    int fq_default_poly_is_gen(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is equal
    # to the polynomial `x`.

    int fq_default_poly_is_unit(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is a unit in the polynomial
    # ring `\mathbf{F}_q[X]`, i.e. if it has degree `0` and is non-zero.

    int fq_default_poly_equal_fq_default(const fq_default_poly_t poly, const fq_default_t c, const fq_default_ctx_t ctx)
    # Returns whether the polynomial ``poly`` is equal the (constant)
    # `\mathbf{F}_q` element ``c``

    void fq_default_poly_add(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    # Sets ``res`` to the sum of ``poly1`` and ``poly2``.

    void fq_default_poly_add_si(fq_default_poly_t res, const fq_default_poly_t poly1, long c, const fq_default_ctx_t ctx)
    # Sets ``res`` to the sum of ``poly1`` and ``c``.

    void fq_default_poly_add_series(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, long n, const fq_default_ctx_t ctx)
    # Notionally truncate ``poly1`` and ``poly2`` to length ``n`` and set
    # ``res`` to the sum.

    void fq_default_poly_sub(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    # Sets ``res`` to the difference of ``poly1`` and ``poly2``.

    void fq_default_poly_sub_series(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, long n, const fq_default_ctx_t ctx)
    # Notionally truncate ``poly1`` and ``poly2`` to length ``n`` and set
    # ``res`` to the difference.

    void fq_default_poly_neg(fq_default_poly_t res, const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Sets ``res`` to the additive inverse of ``poly``.

    void fq_default_poly_scalar_mul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the product of ``op`` by the scalar ``x``, in the context
    # defined by ``ctx``.

    void fq_default_poly_scalar_addmul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    # Adds to ``rop`` the product of ``op`` by the
    # scalar ``x``, in the context defined by ``ctx``.

    void fq_default_poly_scalar_submul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    # Subtracts from ``rop`` the product of ``op`` by the
    # scalar ``x``, in the context defined by ``ctx``.

    void fq_default_poly_scalar_div_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the quotient of ``op`` by the scalar ``x``, in the context
    # defined by ``ctx``. An exception is raised if ``x`` is zero.

    void fq_default_poly_mul(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the product of ``op1`` and ``op2``,
    # choosing an appropriate algorithm.

    void fq_default_poly_mullow(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, long n, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the lowest `n` coefficients of the product of
    # ``op1`` and ``op2``.

    void fq_default_poly_mulhigh(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, long start, const fq_default_ctx_t ctx)
    # Computes the product of ``poly1`` and ``poly2`` and writes the
    # coefficients from ``start`` onwards into the high coefficients of
    # ``res``, the remaining coefficients being arbitrary but reduced.

    void fq_default_poly_mulmod(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    # Sets ``res`` to the remainder of the product of ``poly1``
    # and ``poly2`` upon polynomial division by ``f``.

    void fq_default_poly_sqr(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the square of ``op``,
    # choosing an appropriate algorithm.

    void fq_default_poly_pow(fq_default_poly_t rop, const fq_default_poly_t op, unsigned long e, const fq_default_ctx_t ctx)
    # Computes ``rop = op^e``.  If `e` is zero, returns one,
    # so that in particular ``0^0 = 1``.

    void fq_default_poly_powmod_ui_binexp(fq_default_poly_t res, const fq_default_poly_t poly, unsigned long e, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e >= 0``.

    void fq_default_poly_powmod_fmpz_binexp(fq_default_poly_t res, const fq_default_poly_t poly, const fmpz_t e, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    # Sets ``res`` to ``poly`` raised to the power ``e`` modulo
    # ``f``, using binary exponentiation. We require ``e >= 0``.

    void fq_default_poly_pow_trunc(fq_default_poly_t res, const fq_default_poly_t poly, unsigned long e, long trunc, const fq_default_ctx_t ctx)
    # Sets ``res`` to the low ``trunc`` coefficients of ``poly``
    # to the power ``e``. This is equivalent to doing a powering
    # followed by a truncation.

    void fq_default_poly_shift_left(fq_default_poly_t rop, const fq_default_poly_t op, long n, const fq_default_ctx_t ctx)
    # Sets ``rop`` to ``op`` shifted left by `n` coeffs.  Zero
    # coefficients are inserted.

    void fq_default_poly_shift_right(fq_default_poly_t rop, const fq_default_poly_t op, long n, const fq_default_ctx_t ctx)
    # Sets ``rop`` to ``op`` shifted right by `n` coefficients.
    # If `n` is equal to or greater than the current length of
    # ``op``, ``rop`` is set to the zero polynomial.

    long fq_default_poly_hamming_weight(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Returns the number of non-zero entries in the polynomial ``op``.

    void fq_default_poly_divrem(fq_default_poly_t Q, fq_default_poly_t R, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    # Computes `Q`, `R` such that `A = B Q + R` with
    # `0 \leq \operatorname{len}(R) < \operatorname{len}(B)`.
    # Assumes that the leading coefficient of `B` is invertible.  This can
    # be taken for granted the context is for a finite field, that is, when
    # `p` is prime and `f(X)` is irreducible.

    void fq_default_poly_rem(fq_default_poly_t R, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    # Sets ``R`` to the remainder of the division of ``A`` by
    # ``B`` in the context described by ``ctx``.

    void fq_default_poly_inv_series(fq_default_poly_t Qinv, const fq_default_poly_t Q, long n, const fq_default_ctx_t ctx)
    # Given ``Q`` find ``Qinv`` such that ``Q * Qinv`` is
    # ``1`` modulo `x^n`. The constant coefficient of ``Q`` must
    # be invertible modulo the modulus of ``Q``. An exception is
    # raised if this is not the case or if ``n = 0``.

    void fq_default_poly_div_series(fq_default_poly_t Q, const fq_default_poly_t A, const fq_default_poly_t B, long n, const fq_default_ctx_t ctx)
    # Set `Q` to the quotient of the series `A` by `B`, thinking of the series as
    # though they were of length `n`. We assume that the bottom coefficient of
    # `B` is invertible.

    void fq_default_poly_gcd(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the greatest common divisor of ``op1`` and
    # ``op2``, using the either the Euclidean or HGCD algorithm. The
    # GCD of zero polynomials is defined to be zero, whereas the GCD of
    # the zero polynomial and some other polynomial `P` is defined to be
    # `P`. Except in the case where the GCD is zero, the GCD `G` is made
    # monic.

    void fq_default_poly_xgcd(fq_default_poly_t G, fq_default_poly_t S, fq_default_poly_t T, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    # Computes the GCD of `A` and `B`. The GCD of zero polynomials is
    # defined to be zero, whereas the GCD of the zero polynomial and some other
    # polynomial `P` is defined to be `P`. Except in the case where
    # the GCD is zero, the GCD `G` is made monic.
    # Polynomials ``S`` and ``T`` are computed such that
    # ``S*A + T*B = G``. The length of ``S`` will be at most
    # ``lenB`` and the length of ``T`` will be at most ``lenA``.

    int fq_default_poly_divides(fq_default_poly_t Q, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    # Returns `1` if `B` divides `A` exactly and sets `Q` to the quotient,
    # otherwise returns `0`.
    # This function is currently unoptimised and provided for convenience
    # only.

    void fq_default_poly_derivative(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the derivative of ``op``.

    void fq_default_poly_invsqrt_series(fq_default_poly_t g, const fq_default_poly_t h, long n, fq_default_ctx_t ctx)
    # Set `g` to the series expansion of `1/\sqrt{h}` to order `O(x^n)`.
    # It is assumed that `h` has constant term 1.

    void fq_default_poly_sqrt_series(fq_default_poly_t g, const fq_default_poly_t h, long n, fq_default_ctx_t ctx)
    # Set `g` to the series expansion of `\sqrt{h}` to order `O(x^n)`.
    # It is assumed that `h` has constant term 1.

    int fq_default_poly_sqrt(fq_default_poly_t s, const fq_default_poly_t p, fq_default_ctx_t mod)
    # If `p` is a perfect square, sets `s` to a square root of `p`
    # and returns 1. Otherwise returns 0.

    void fq_default_poly_evaluate_fq_default(fq_default_t rop, const fq_default_poly_t f, const fq_default_t a, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the value of `f(a)`.
    # As the coefficient ring `\mathbf{F}_q` is finite, Horner's method
    # is sufficient.

    void fq_default_poly_compose(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
    # Sets ``rop`` to the composition of ``op1`` and ``op2``.
    # To be precise about the order of composition, denoting ``rop``,
    # ``op1``, and ``op2`` by `f`, `g`, and `h`, respectively,
    # sets `f(t) = g(h(t))`.

    void fq_default_poly_compose_mod(fq_default_poly_t res, const fq_default_poly_t f, const fq_default_poly_t g, const fq_default_poly_t h, const fq_default_ctx_t ctx)
    # Sets ``res`` to the composition `f(g)` modulo `h`. We require
    # that `h` is nonzero.

    int fq_default_poly_fprint_pretty(FILE * file, const fq_default_poly_t poly, const char *x, const fq_default_ctx_t ctx)
    # Prints the pretty representation of ``poly`` to the stream
    # ``file``, using the string ``x`` to represent the indeterminate.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_default_poly_print_pretty(const fq_default_poly_t poly, const char *x, const fq_default_ctx_t ctx)
    # Prints the pretty representation of ``poly`` to ``stdout``,
    # using the string ``x`` to represent the indeterminate.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_default_poly_fprint(FILE * file, const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Prints the pretty representation of ``poly`` to the stream
    # ``file``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int fq_default_poly_print(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Prints the representation of ``poly`` to ``stdout``.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    char * fq_default_poly_get_str(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    # Returns the plain FLINT string representation of the polynomial
    # ``poly``.

    char * fq_default_poly_get_str_pretty(const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)
    # Returns a pretty representation of the polynomial ``poly`` using the
    # null-terminated string ``x`` as the variable name

    void fq_default_poly_inflate(fq_default_poly_t result, const fq_default_poly_t input, unsigned long inflation, const fq_default_ctx_t ctx)
    # Sets ``result`` to the inflated polynomial `p(x^n)` where
    # `p` is given by ``input`` and `n` is given by ``inflation``.

    void fq_default_poly_deflate(fq_default_poly_t result, const fq_default_poly_t input, unsigned long deflation, const fq_default_ctx_t ctx)
    # Sets ``result`` to the deflated polynomial `p(x^{1/n})` where
    # `p` is given by ``input`` and `n` is given by ``deflation``.
    # Requires `n > 0`.

    unsigned long fq_default_poly_deflation(const fq_default_poly_t input, const fq_default_ctx_t ctx)
    # Returns the largest integer by which ``input`` can be deflated.
    # As special cases, returns 0 if ``input`` is the zero polynomial
    # and 1 of ``input`` is a constant polynomial.
