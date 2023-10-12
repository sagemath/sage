# distutils: libraries = flint
# distutils: depends = flint/ca_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void ca_poly_init(ca_poly_t poly, ca_ctx_t ctx)
    # Initializes the polynomial for use, setting it to the zero polynomial.

    void ca_poly_clear(ca_poly_t poly, ca_ctx_t ctx)
    # Clears the polynomial, deallocating all coefficients and the
    # coefficient array.

    void ca_poly_fit_length(ca_poly_t poly, long len, ca_ctx_t ctx)
    # Makes sure that the coefficient array of the polynomial contains at
    # least *len* initialized coefficients.

    void _ca_poly_set_length(ca_poly_t poly, long len, ca_ctx_t ctx)
    # Directly changes the length of the polynomial, without allocating or
    # deallocating coefficients. The value should not exceed the allocation length.

    void _ca_poly_normalise(ca_poly_t poly, ca_ctx_t ctx)
    # Strips any top coefficients which can be proved identical to zero.

    void ca_poly_zero(ca_poly_t poly, ca_ctx_t ctx)
    # Sets *poly* to the zero polynomial.

    void ca_poly_one(ca_poly_t poly, ca_ctx_t ctx)
    # Sets *poly* to the constant polynomial 1.

    void ca_poly_x(ca_poly_t poly, ca_ctx_t ctx)
    # Sets *poly* to the monomial *x*.

    void ca_poly_set_ca(ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
    void ca_poly_set_si(ca_poly_t poly, long c, ca_ctx_t ctx)
    # Sets *poly* to the constant polynomial *c*.

    void ca_poly_set(ca_poly_t res, const ca_poly_t src, ca_ctx_t ctx)
    void ca_poly_set_fmpz_poly(ca_poly_t res, const fmpz_poly_t src, ca_ctx_t ctx)
    void ca_poly_set_fmpq_poly(ca_poly_t res, const fmpq_poly_t src, ca_ctx_t ctx)
    # Sets *poly* the polynomial *src*.

    void ca_poly_set_coeff_ca(ca_poly_t poly, long n, const ca_t x, ca_ctx_t ctx)
    # Sets the coefficient at position *n* in *poly* to *x*.

    void ca_poly_transfer(ca_poly_t res, ca_ctx_t res_ctx, const ca_poly_t src, ca_ctx_t src_ctx)
    # Sets *res* to *src* where the corresponding context objects *res_ctx* and
    # *src_ctx* may be different.
    # This operation preserves the mathematical value represented by *src*,
    # but may result in a different internal representation depending on the
    # settings of the context objects.

    void ca_poly_randtest(ca_poly_t poly, flint_rand_t state, long len, long depth, long bits, ca_ctx_t ctx)
    # Sets *poly* to a random polynomial of length up to *len* and with entries having complexity up to
    # *depth* and *bits* (see :func:`ca_randtest`).

    void ca_poly_randtest_rational(ca_poly_t poly, flint_rand_t state, long len, long bits, ca_ctx_t ctx)
    # Sets *poly* to a random rational polynomial of length up to *len* and with entries up to *bits* bits in size.

    void ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx)
    # Prints *poly* to standard output. The coefficients are printed on separate lines.

    void ca_poly_printn(const ca_poly_t poly, long digits, ca_ctx_t ctx)
    # Prints a decimal representation of *poly* with precision specified by *digits*.
    # The coefficients are comma-separated and the whole list is enclosed in square brackets.

    int ca_poly_is_proper(const ca_poly_t poly, ca_ctx_t ctx)
    # Checks that *poly* represents an element of `\mathbb{C}[X]` with
    # well-defined degree. This returns 1 if the leading coefficient
    # of *poly* is nonzero and all coefficients of *poly* are
    # numbers (not special values). It returns 0 otherwise.
    # It returns 1 when *poly* is precisely the zero polynomial (which
    # does not have a leading coefficient).

    int ca_poly_make_monic(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
    # Makes *poly* monic by dividing by the leading coefficient if possible
    # and returns 1. Returns 0 if the leading coefficient cannot be
    # certified to be nonzero, or if *poly* is the zero polynomial.

    void _ca_poly_reverse(ca_ptr res, ca_srcptr poly, long len, long n, ca_ctx_t ctx)

    void ca_poly_reverse(ca_poly_t res, const ca_poly_t poly, long n, ca_ctx_t ctx)
    # Sets *res* to the reversal of *poly* considered as a polynomial
    # of length *n*, zero-padding if needed. The underscore method
    # assumes that *len* is positive and less than or equal to *n*.

    truth_t _ca_poly_check_equal(ca_srcptr poly1, long len1, ca_srcptr poly2, long len2, ca_ctx_t ctx)
    truth_t ca_poly_check_equal(const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
    # Checks if *poly1* and *poly2* represent the same polynomial.
    # The underscore method assumes that *len1* is at least as
    # large as *len2*.

    truth_t ca_poly_check_is_zero(const ca_poly_t poly, ca_ctx_t ctx)
    # Checks if *poly* is the zero polynomial.

    truth_t ca_poly_check_is_one(const ca_poly_t poly, ca_ctx_t ctx)
    # Checks if *poly* is the constant polynomial 1.

    void _ca_poly_shift_left(ca_ptr res, ca_srcptr poly, long len, long n, ca_ctx_t ctx)
    void ca_poly_shift_left(ca_poly_t res, const ca_poly_t poly, long n, ca_ctx_t ctx)
    # Sets *res* to *poly* shifted *n* coefficients to the left; that is,
    # multiplied by `x^n`.

    void _ca_poly_shift_right(ca_ptr res, ca_srcptr poly, long len, long n, ca_ctx_t ctx)
    void ca_poly_shift_right(ca_poly_t res, const ca_poly_t poly, long n, ca_ctx_t ctx)
    # Sets *res* to *poly* shifted *n* coefficients to the right; that is,
    # divided by `x^n`.

    void ca_poly_neg(ca_poly_t res, const ca_poly_t src, ca_ctx_t ctx)
    # Sets *res* to the negation of *src*.

    void _ca_poly_add(ca_ptr res, ca_srcptr poly1, long len1, ca_srcptr poly2, long len2, ca_ctx_t ctx)
    void ca_poly_add(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
    # Sets *res* to the sum of *poly1* and *poly2*.

    void _ca_poly_sub(ca_ptr res, ca_srcptr poly1, long len1, ca_srcptr poly2, long len2, ca_ctx_t ctx)
    void ca_poly_sub(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
    # Sets *res* to the difference of *poly1* and *poly2*.

    void _ca_poly_mul(ca_ptr res, ca_srcptr poly1, long len1, ca_srcptr poly2, long len2, ca_ctx_t ctx)
    void ca_poly_mul(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
    # Sets *res* to the product of *poly1* and *poly2*.

    void _ca_poly_mullow(ca_ptr C, ca_srcptr poly1, long len1, ca_srcptr poly2, long len2, long n, ca_ctx_t ctx)
    void ca_poly_mullow(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, long n, ca_ctx_t ctx)
    # Sets *res* to the product of *poly1* and *poly2* truncated to length *n*.

    void ca_poly_mul_ca(ca_poly_t res, const ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
    # Sets *res* to *poly* multiplied by the scalar *c*.

    void ca_poly_div_ca(ca_poly_t res, const ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
    # Sets *res* to *poly* divided by the scalar *c*.

    void _ca_poly_divrem_basecase(ca_ptr Q, ca_ptr R, ca_srcptr A, long lenA, ca_srcptr B, long lenB, const ca_t invB, ca_ctx_t ctx)
    int ca_poly_divrem_basecase(ca_poly_t Q, ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
    void _ca_poly_divrem(ca_ptr Q, ca_ptr R, ca_srcptr A, long lenA, ca_srcptr B, long lenB, const ca_t invB, ca_ctx_t ctx)
    int ca_poly_divrem(ca_poly_t Q, ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
    int ca_poly_div(ca_poly_t Q, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
    int ca_poly_rem(ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
    # If the leading coefficient of *B* can be proved invertible, sets *Q* and *R*
    # to the quotient and remainder of polynomial division of *A* by *B*
    # and returns 1. If the leading coefficient cannot be proved
    # invertible, returns 0.
    # The underscore method takes a precomputed inverse of the leading coefficient of *B*.

    void _ca_poly_pow_ui_trunc(ca_ptr res, ca_srcptr f, long flen, unsigned long exp, long len, ca_ctx_t ctx)
    void ca_poly_pow_ui_trunc(ca_poly_t res, const ca_poly_t poly, unsigned long exp, long len, ca_ctx_t ctx)
    # Sets *res* to *poly* raised to the power *exp*, truncated to length *len*.

    void _ca_poly_pow_ui(ca_ptr res, ca_srcptr f, long flen, unsigned long exp, ca_ctx_t ctx)
    void ca_poly_pow_ui(ca_poly_t res, const ca_poly_t poly, unsigned long exp, ca_ctx_t ctx)
    # Sets *res* to *poly* raised to the power *exp*.

    void _ca_poly_evaluate_horner(ca_t res, ca_srcptr f, long len, const ca_t x, ca_ctx_t ctx)
    void ca_poly_evaluate_horner(ca_t res, const ca_poly_t f, const ca_t a, ca_ctx_t ctx)
    void _ca_poly_evaluate(ca_t res, ca_srcptr f, long len, const ca_t x, ca_ctx_t ctx)
    void ca_poly_evaluate(ca_t res, const ca_poly_t f, const ca_t a, ca_ctx_t ctx)
    # Sets *res* to *f* evaluated at the point *a*.

    void _ca_poly_compose(ca_ptr res, ca_srcptr poly1, long len1, ca_srcptr poly2, long len2, ca_ctx_t ctx)
    void ca_poly_compose(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
    # Sets *res* to the composition of *poly1* with *poly2*.

    void _ca_poly_derivative(ca_ptr res, ca_srcptr poly, long len, ca_ctx_t ctx)
    void ca_poly_derivative(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
    # Sets *res* to the derivative of *poly*. The underscore method needs one less
    # coefficient than *len* for the output array.

    void _ca_poly_integral(ca_ptr res, ca_srcptr poly, long len, ca_ctx_t ctx)
    void ca_poly_integral(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
    # Sets *res* to the integral of *poly*. The underscore method needs one more
    # coefficient than *len* for the output array.

    void _ca_poly_inv_series(ca_ptr res, ca_srcptr f, long flen, long len, ca_ctx_t ctx)
    void ca_poly_inv_series(ca_poly_t res, const ca_poly_t f, long len, ca_ctx_t ctx)
    # Sets *res* to the power series inverse of *f* truncated
    # to length *len*.

    void _ca_poly_div_series(ca_ptr res, ca_srcptr f, long flen, ca_srcptr g, long glen, long len, ca_ctx_t ctx)
    void ca_poly_div_series(ca_poly_t res, const ca_poly_t f, const ca_poly_t g, long len, ca_ctx_t ctx)
    # Sets *res* to the power series quotient of *f* and *g* truncated
    # to length *len*.
    # This function divides by zero if *g* has constant term zero;
    # the user should manually remove initial zeros when an
    # exact cancellation is required.

    void _ca_poly_exp_series(ca_ptr res, ca_srcptr f, long flen, long len, ca_ctx_t ctx)
    void ca_poly_exp_series(ca_poly_t res, const ca_poly_t f, long len, ca_ctx_t ctx)
    # Sets *res* to the power series exponential of *f* truncated
    # to length *len*.

    void _ca_poly_log_series(ca_ptr res, ca_srcptr f, long flen, long len, ca_ctx_t ctx)
    void ca_poly_log_series(ca_poly_t res, const ca_poly_t f, long len, ca_ctx_t ctx)
    # Sets *res* to the power series logarithm of *f* truncated
    # to length *len*.

    long _ca_poly_gcd_euclidean(ca_ptr res, ca_srcptr A, long lenA, ca_srcptr B, long lenB, ca_ctx_t ctx)
    int ca_poly_gcd_euclidean(ca_poly_t res, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
    long _ca_poly_gcd(ca_ptr res, ca_srcptr A, long lenA, ca_srcptr B, long lenB, ca_ctx_t ctx)
    int ca_poly_gcd(ca_poly_t res, const ca_poly_t A, const ca_poly_t g, ca_ctx_t ctx)
    # Sets *res* to the GCD of *A* and *B* and returns 1 on success.
    # On failure, returns 0 leaving the value of *res* arbitrary.
    # The computation can fail if testing a leading coefficient
    # for zero fails in the execution of the GCD algorithm.
    # The output is normalized to be monic if it is not the zero polynomial.
    # The underscore methods assume `\text{lenA} \ge \text{lenB} \ge 1`, and that
    # both *A* and *B* have nonzero leading coefficient.
    # They return the length of the GCD, or 0 if the computation fails.
    # The *euclidean* version implements the standard Euclidean algorithm.
    # The default version first checks for rational polynomials or
    # attempts to certify numerically that the polynomials are coprime
    # and otherwise falls back to an automatic choice
    # of algorithm (currently only the Euclidean algorithm).

    int ca_poly_factor_squarefree(ca_t c, ca_poly_vec_t fac, unsigned long * exp, const ca_poly_t F, ca_ctx_t ctx)
    # Computes the squarefree factorization of *F*, giving a product
    # `F = c f_1 f_2^2 \ldots f_n^n` where all `f_i` with `f_i \ne 1`
    # are squarefree and pairwise coprime. The nontrivial factors
    # `f_i` are written to *fac* and the corresponding exponents
    # are written to *exp*. This algorithm can fail if GCD computation
    # fails internally. Returns 1 on success and 0 on failure.

    int ca_poly_squarefree_part(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
    # Sets *res* to the squarefree part of *poly*, normalized to be monic.
    # This algorithm can fail if GCD computation fails internally.
    # Returns 1 on success and 0 on failure.

    void _ca_poly_set_roots(ca_ptr poly, ca_srcptr roots, const unsigned long * exp, long n, ca_ctx_t ctx)
    void ca_poly_set_roots(ca_poly_t poly, ca_vec_t roots, const unsigned long * exp, ca_ctx_t ctx)
    # Sets *poly* to the monic polynomial with the *n* roots
    # given in the vector *roots*, with multiplicities given
    # in the vector *exp*. In other words, this constructs
    # the polynomial
    # `(x-r_0)^{e_0} (x-r_1)^{e_1} \cdots (x-r_{n-1})^{e_{n-1}}`.
    # Uses binary splitting.

    int _ca_poly_roots(ca_ptr roots, ca_srcptr poly, long len, ca_ctx_t ctx)
    int ca_poly_roots(ca_vec_t roots, unsigned long * exp, const ca_poly_t poly, ca_ctx_t ctx)
    # Attempts to compute all complex roots of the given polynomial *poly*.
    # On success, returns 1 and sets *roots* to a vector containing all
    # the distinct roots with corresponding multiplicities in *exp*.
    # On failure, returns 0 and leaves the values in *roots* arbitrary.
    # The roots are returned in arbitrary order.
    # Failure will occur if the leading coefficient of *poly* cannot
    # be proved to be nonzero, if determining the correct multiplicities
    # fails, or if the builtin algorithms do not have a means to
    # represent the roots symbolically.
    # The underscore method assumes that the polynomial is squarefree.
    # The non-underscore method performs a squarefree factorization.

    ca_poly_struct * _ca_poly_vec_init(long len, ca_ctx_t ctx)
    void ca_poly_vec_init(ca_poly_vec_t res, long len, ca_ctx_t ctx)
    # Initializes a vector with *len* polynomials.

    void _ca_poly_vec_fit_length(ca_poly_vec_t vec, long len, ca_ctx_t ctx)
    # Allocates space for *len* polynomials in *vec*.

    void ca_poly_vec_set_length(ca_poly_vec_t vec, long len, ca_ctx_t ctx)
    # Resizes *vec* to length *len*, zero-extending if needed.

    void _ca_poly_vec_clear(ca_poly_struct * vec, long len, ca_ctx_t ctx)
    void ca_poly_vec_clear(ca_poly_vec_t vec, ca_ctx_t ctx)
    # Clears the vector *vec*.

    void ca_poly_vec_append(ca_poly_vec_t vec, const ca_poly_t poly, ca_ctx_t ctx)
    # Appends *poly* to the end of the vector *vec*.
