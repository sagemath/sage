# distutils: libraries = flint
# distutils: depends = flint/fmpz_mpoly_q.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_mpoly_q_init(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
    # Initializes *res* for use, and sets its value to zero.

    void fmpz_mpoly_q_clear(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
    # Clears *res*, freeing or recycling its allocated memory.

    void fmpz_mpoly_q_swap(fmpz_mpoly_q_t x, fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
    # Swaps the values of *x* and *y* efficiently.

    void fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_set_fmpq(fmpz_mpoly_q_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_set_fmpz(fmpz_mpoly_q_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_set_si(fmpz_mpoly_q_t res, slong x, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the value *x*.

    void fmpz_mpoly_q_canonicalise(fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Puts the numerator and denominator of *x* in canonical form by removing
    # common content and making the leading term of the denominator positive.

    int fmpz_mpoly_q_is_canonical(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Returns whether *x* is in canonical form.
    # In addition to verifying that the numerator and denominator
    # have no common content and that the leading term of the denominator
    # is positive, this function checks that the denominator is nonzero and that
    # the numerator and denominator have correctly sorted terms
    # (these properties should normally hold; verifying them
    # provides an extra consistency check for test code).

    int fmpz_mpoly_q_is_zero(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Returns whether *x* is the constant 0.

    int fmpz_mpoly_q_is_one(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Returns whether *x* is the constant 1.

    void fmpz_mpoly_q_used_vars(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_used_vars_num(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_used_vars_den(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
    # For each variable, sets the corresponding entry in *used* to the
    # boolean flag indicating whether that variable appears in the
    # rational function (respectively its numerator or denominator).

    void fmpz_mpoly_q_zero(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the constant 0.

    void fmpz_mpoly_q_one(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the constant 1.

    void fmpz_mpoly_q_gen(fmpz_mpoly_q_t res, slong i, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the generator `x_{i+1}`.
    # Requires `0 \le i < n` where *n* is the number of variables of *ctx*.

    void fmpz_mpoly_q_print_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)
    # Prints *res* to standard output. If *x* is not *NULL*, the strings in
    # *x* are used as the symbols for the variables.

    void fmpz_mpoly_q_randtest(fmpz_mpoly_q_t res, flint_rand_t state, slong length, mp_limb_t coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to a random rational function where both numerator and denominator
    # have up to *length* terms, coefficients up to size *coeff_bits*, and
    # exponents strictly smaller than *exp_bound*.

    int fmpz_mpoly_q_equal(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
    # Returns whether *x* and *y* are equal.

    void fmpz_mpoly_q_neg(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the negation of *x*.

    void fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_add_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_add_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the sum of *x* and *y*.

    void fmpz_mpoly_q_sub(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_sub_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_sub_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the difference of *x* and *y*.

    void fmpz_mpoly_q_mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_mul_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_mul_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the product of *x* and *y*.

    void fmpz_mpoly_q_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_div_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_div_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the quotient of *x* and *y*.
    # Division by zero calls *flint_abort*.

    void fmpz_mpoly_q_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the inverse of *x*. Division by zero
    # calls *flint_abort*.

    void _fmpz_mpoly_q_content(fmpz_t num, fmpz_t den, const fmpz_mpoly_t xnum, const fmpz_mpoly_t xden, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_q_content(fmpq_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the content of the coefficients of *x*.
