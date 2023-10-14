# distutils: libraries = flint
# distutils: depends = flint/qqbar.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void qqbar_init(qqbar_t res)
    # Initializes the variable *res* for use, and sets its value to zero.

    void qqbar_clear(qqbar_t res)
    # Clears the variable *res*, freeing or recycling its allocated memory.

    qqbar_ptr _qqbar_vec_init(slong len)
    # Returns a pointer to an array of *len* initialized *qqbar_struct*:s.

    void _qqbar_vec_clear(qqbar_ptr vec, slong len)
    # Clears all *len* entries in the vector *vec* and frees the
    # vector itself.

    void qqbar_swap(qqbar_t x, qqbar_t y)
    # Swaps the values of *x* and *y* efficiently.

    void qqbar_set(qqbar_t res, const qqbar_t x)
    void qqbar_set_si(qqbar_t res, slong x)
    void qqbar_set_ui(qqbar_t res, ulong x)
    void qqbar_set_fmpz(qqbar_t res, const fmpz_t x)
    void qqbar_set_fmpq(qqbar_t res, const fmpq_t x)
    # Sets *res* to the value *x*.

    void qqbar_set_re_im(qqbar_t res, const qqbar_t x, const qqbar_t y)
    # Sets *res* to the value `x + yi`.

    int qqbar_set_d(qqbar_t res, double x)
    int qqbar_set_re_im_d(qqbar_t res, double x, double y)
    # Sets *res* to the value *x* or `x + yi` respectively. These functions
    # performs error handling: if *x* and *y* are finite, the conversion succeeds
    # and the return flag is 1. If *x* or *y* is non-finite (infinity or NaN),
    # the conversion fails and the return flag is 0.

    slong qqbar_degree(const qqbar_t x)
    # Returns the degree of *x*, i.e. the degree of the minimal polynomial.

    int qqbar_is_rational(const qqbar_t x)
    # Returns whether *x* is a rational number.

    int qqbar_is_integer(const qqbar_t x)
    # Returns whether *x* is an integer (an element of `\mathbb{Z}`).

    int qqbar_is_algebraic_integer(const qqbar_t x)
    # Returns whether *x* is an algebraic integer, i.e. whether its minimal
    # polynomial has leading coefficient 1.

    int qqbar_is_zero(const qqbar_t x)
    int qqbar_is_one(const qqbar_t x)
    int qqbar_is_neg_one(const qqbar_t x)
    # Returns whether *x* is the number `0`, `1`, `-1`.

    int qqbar_is_i(const qqbar_t x)
    int qqbar_is_neg_i(const qqbar_t x)
    # Returns whether *x* is the imaginary unit `i` (respectively `-i`).

    int qqbar_is_real(const qqbar_t x)
    # Returns whether *x* is a real number.

    void qqbar_height(fmpz_t res, const qqbar_t x)
    # Sets *res* to the height of *x* (the largest absolute value of the
    # coefficients of the minimal polynomial of *x*).

    slong qqbar_height_bits(const qqbar_t x)
    # Returns the height of *x* (the largest absolute value of the
    # coefficients of the minimal polynomial of *x*) measured in bits.

    int qqbar_within_limits(const qqbar_t x, slong deg_limit, slong bits_limit)
    # Checks if *x* has degree bounded by *deg_limit* and height
    # bounded by *bits_limit* bits, returning 0 (false) or 1 (true).
    # If *deg_limit* is set to 0, the degree check is skipped,
    # and similarly for *bits_limit*.

    int qqbar_binop_within_limits(const qqbar_t x, const qqbar_t y, slong deg_limit, slong bits_limit)
    # Checks if `x + y`, `x - y`, `x \cdot y` and `x / y` certainly have
    # degree bounded by *deg_limit* (by multiplying the degrees for *x* and *y*
    # to obtain a trivial bound). For *bits_limits*, the sum of the bit heights
    # of *x* and *y* is checked against the bound (this is only a heuristic).
    # If *deg_limit* is set to 0, the degree check is skipped,
    # and similarly for *bits_limit*.

    void _qqbar_get_fmpq(fmpz_t num, fmpz_t den, const qqbar_t x)
    # Sets *num* and *den* to the numerator and denominator of *x*.
    # Aborts if *x* is not a rational number.

    void qqbar_get_fmpq(fmpq_t res, const qqbar_t x)
    # Sets *res* to *x*. Aborts if *x* is not a rational number.

    void qqbar_get_fmpz(fmpz_t res, const qqbar_t x)
    # Sets *res* to *x*. Aborts if *x* is not an integer.

    void qqbar_zero(qqbar_t res)
    # Sets *res* to the number 0.

    void qqbar_one(qqbar_t res)
    # Sets *res* to the number 1.

    void qqbar_i(qqbar_t res)
    # Sets *res* to the imaginary unit `i`.

    void qqbar_phi(qqbar_t res)
    # Sets *res* to the golden ratio `\varphi = \tfrac{1}{2}(\sqrt{5} + 1)`.

    void qqbar_print(const qqbar_t x)
    # Prints *res* to standard output. The output shows the degree
    # and the list of coefficients
    # of the minimal polynomial followed by a decimal representation of
    # the enclosing interval. This function is mainly intended for debugging.

    void qqbar_printn(const qqbar_t x, slong n)
    # Prints *res* to standard output. The output shows a decimal
    # approximation to *n* digits.

    void qqbar_printnd(const qqbar_t x, slong n)
    # Prints *res* to standard output. The output shows a decimal
    # approximation to *n* digits, followed by the degree of the number.

    void qqbar_randtest(qqbar_t res, flint_rand_t state, slong deg, slong bits)
    # Sets *res* to a random algebraic number with degree up to *deg* and
    # with height (measured in bits) up to *bits*.

    void qqbar_randtest_real(qqbar_t res, flint_rand_t state, slong deg, slong bits)
    # Sets *res* to a random real algebraic number with degree up to *deg* and
    # with height (measured in bits) up to *bits*.

    void qqbar_randtest_nonreal(qqbar_t res, flint_rand_t state, slong deg, slong bits)
    # Sets *res* to a random nonreal algebraic number with degree up to *deg* and
    # with height (measured in bits) up to *bits*. Since all algebraic numbers
    # of degree 1 are real, *deg* must be at least 2.

    int qqbar_equal(const qqbar_t x, const qqbar_t y)
    # Returns whether *x* and *y* are equal.

    int qqbar_equal_fmpq_poly_val(const qqbar_t x, const fmpq_poly_t f, const qqbar_t y)
    # Returns whether *x* is equal to `f(y)`. This function is more efficient
    # than evaluating `f(y)` and comparing the results.

    int qqbar_cmp_re(const qqbar_t x, const qqbar_t y)
    # Compares the real parts of *x* and *y*, returning -1, 0 or +1.

    int qqbar_cmp_im(const qqbar_t x, const qqbar_t y)
    # Compares the imaginary parts of *x* and *y*, returning -1, 0 or +1.

    int qqbar_cmpabs_re(const qqbar_t x, const qqbar_t y)
    # Compares the absolute values of the real parts of *x* and *y*, returning -1, 0 or +1.

    int qqbar_cmpabs_im(const qqbar_t x, const qqbar_t y)
    # Compares the absolute values of the imaginary parts of *x* and *y*, returning -1, 0 or +1.

    int qqbar_cmpabs(const qqbar_t x, const qqbar_t y)
    # Compares the absolute values of *x* and *y*, returning -1, 0 or +1.

    int qqbar_cmp_root_order(const qqbar_t x, const qqbar_t y)
    # Compares *x* and *y* using an arbitrary but convenient ordering
    # defined on the complex numbers. This is useful for sorting the
    # roots of a polynomial in a canonical order.
    # We define the root order as follows: real roots come first, in
    # descending order. Nonreal roots are subsequently ordered first by
    # real part in descending order, then in ascending order by the
    # absolute value of the imaginary part, and then in descending
    # order of the sign. This implies that complex conjugate roots
    # are adjacent, with the root in the upper half plane first.

    ulong qqbar_hash(const qqbar_t x)
    # Returns a hash of *x*. As currently implemented, this function
    # only hashes the minimal polynomial of *x*. The user should
    # mix in some bits based on the numerical value if it is critical
    # to distinguish between conjugates of the same minimal polynomial.
    # This function is also likely to produce serial runs of values for
    # lexicographically close minimal polynomials. This is not
    # necessarily a problem for use in hash tables, but if it is
    # important that all bits in the output are random,
    # the user should apply an integer hash function to the output.

    void qqbar_conj(qqbar_t res, const qqbar_t x)
    # Sets *res* to the complex conjugate of *x*.

    void qqbar_re(qqbar_t res, const qqbar_t x)
    # Sets *res* to the real part of *x*.

    void qqbar_im(qqbar_t res, const qqbar_t x)
    # Sets *res* to the imaginary part of *x*.

    void qqbar_re_im(qqbar_t res1, qqbar_t res2, const qqbar_t x)
    # Sets *res1* to the real part of *x* and *res2* to the imaginary part of *x*.

    void qqbar_abs(qqbar_t res, const qqbar_t x)
    # Sets *res* to the absolute value of *x*:

    void qqbar_abs2(qqbar_t res, const qqbar_t x)
    # Sets *res* to the square of the absolute value of *x*.

    void qqbar_sgn(qqbar_t res, const qqbar_t x)
    # Sets *res* to the complex sign of *x*, defined as 0 if *x* is zero
    # and as `x / |x|` otherwise.

    int qqbar_sgn_re(const qqbar_t x)
    # Returns the sign of the real part of *x* (-1, 0 or +1).

    int qqbar_sgn_im(const qqbar_t x)
    # Returns the sign of the imaginary part of *x* (-1, 0 or +1).

    int qqbar_csgn(const qqbar_t x)
    # Returns the extension of the real sign function taking the
    # value 1 for *x* strictly in the right half plane, -1 for *x* strictly
    # in the left half plane, and the sign of the imaginary part when *x* is on
    # the imaginary axis. Equivalently, `\operatorname{csgn}(x) = x / \sqrt{x^2}`
    # except that the value is 0 when *x* is zero.

    void qqbar_floor(fmpz_t res, const qqbar_t x)
    # Sets *res* to the floor function of *x*. If *x* is not real, the
    # value is defined as the floor function of the real part of *x*.

    void qqbar_ceil(fmpz_t res, const qqbar_t x)
    # Sets *res* to the ceiling function of *x*. If *x* is not real, the
    # value is defined as the ceiling function of the real part of *x*.

    void qqbar_neg(qqbar_t res, const qqbar_t x)
    # Sets *res* to the negation of *x*.

    void qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_add_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
    void qqbar_add_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
    void qqbar_add_ui(qqbar_t res, const qqbar_t x, ulong y)
    void qqbar_add_si(qqbar_t res, const qqbar_t x, slong y)
    # Sets *res* to the sum of *x* and *y*.

    void qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_sub_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
    void qqbar_sub_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
    void qqbar_sub_ui(qqbar_t res, const qqbar_t x, ulong y)
    void qqbar_sub_si(qqbar_t res, const qqbar_t x, slong y)
    void qqbar_fmpq_sub(qqbar_t res, const fmpq_t x, const qqbar_t y)
    void qqbar_fmpz_sub(qqbar_t res, const fmpz_t x, const qqbar_t y)
    void qqbar_ui_sub(qqbar_t res, ulong x, const qqbar_t y)
    void qqbar_si_sub(qqbar_t res, slong x, const qqbar_t y)
    # Sets *res* to the difference of *x* and *y*.

    void qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_mul_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
    void qqbar_mul_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
    void qqbar_mul_ui(qqbar_t res, const qqbar_t x, ulong y)
    void qqbar_mul_si(qqbar_t res, const qqbar_t x, slong y)
    # Sets *res* to the product of *x* and *y*.

    void qqbar_mul_2exp_si(qqbar_t res, const qqbar_t x, slong e)
    # Sets *res* to *x* multiplied by `2^e`.

    void qqbar_sqr(qqbar_t res, const qqbar_t x)
    # Sets *res* to the square of *x*.

    void qqbar_inv(qqbar_t res, const qqbar_t x)
    # Sets *res* to the multiplicative inverse of *y*.
    # Division by zero calls *flint_abort*.

    void qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_div_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
    void qqbar_div_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
    void qqbar_div_ui(qqbar_t res, const qqbar_t x, ulong y)
    void qqbar_div_si(qqbar_t res, const qqbar_t x, slong y)
    void qqbar_fmpq_div(qqbar_t res, const fmpq_t x, const qqbar_t y)
    void qqbar_fmpz_div(qqbar_t res, const fmpz_t x, const qqbar_t y)
    void qqbar_ui_div(qqbar_t res, ulong x, const qqbar_t y)
    void qqbar_si_div(qqbar_t res, slong x, const qqbar_t y)
    # Sets *res* to the quotient of *x* and *y*.
    # Division by zero calls *flint_abort*.

    void qqbar_scalar_op(qqbar_t res, const qqbar_t x, const fmpz_t a, const fmpz_t b, const fmpz_t c)
    # Sets *res* to the rational affine transformation `(ax+b)/c`, performed as
    # a single operation. There are no restrictions on *a*, *b* and *c*
    # except that *c* must be nonzero. Division by zero calls *flint_abort*.

    void qqbar_sqrt(qqbar_t res, const qqbar_t x)
    void qqbar_sqrt_ui(qqbar_t res, ulong x)
    # Sets *res* to the principal square root of *x*.

    void qqbar_rsqrt(qqbar_t res, const qqbar_t x)
    # Sets *res* to the reciprocal of the principal square root of *x*.
    # Division by zero calls *flint_abort*.

    void qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong n)
    void qqbar_pow_si(qqbar_t res, const qqbar_t x, slong n)
    void qqbar_pow_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t n)
    void qqbar_pow_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t n)
    # Sets *res* to *x* raised to the *n*-th power.
    # Raising zero to a negative power aborts.

    void qqbar_root_ui(qqbar_t res, const qqbar_t x, ulong n)
    void qqbar_fmpq_root_ui(qqbar_t res, const fmpq_t x, ulong n)
    # Sets *res* to the principal *n*-th root of *x*. The order *n*
    # must be positive.

    void qqbar_fmpq_pow_si_ui(qqbar_t res, const fmpq_t x, slong m, ulong n)
    # Sets *res* to the principal branch of `x^{m/n}`. The order *n*
    # must be positive. Division by zero calls *flint_abort*.

    int qqbar_pow(qqbar_t res, const qqbar_t x, const qqbar_t y)
    # General exponentiation: if `x^y` is an algebraic number, sets *res*
    # to this value and returns 1. If `x^y` is transcendental or
    # undefined, returns 0. Note that this function returns 0 instead of
    # aborting on division zero.

    void qqbar_get_acb(acb_t res, const qqbar_t x, slong prec)
    # Sets *res* to an enclosure of *x* rounded to *prec* bits.

    void qqbar_get_arb(arb_t res, const qqbar_t x, slong prec)
    # Sets *res* to an enclosure of *x* rounded to *prec* bits, assuming that
    # *x* is a real number. If *x* is not real, *res* is set to
    # `[\operatorname{NaN} \pm \infty]`.

    void qqbar_get_arb_re(arb_t res, const qqbar_t x, slong prec)
    # Sets *res* to an enclosure of the real part of *x* rounded to *prec* bits.

    void qqbar_get_arb_im(arb_t res, const qqbar_t x, slong prec)
    # Sets *res* to an enclosure of the imaginary part of *x* rounded to *prec* bits.

    void qqbar_cache_enclosure(qqbar_t res, slong prec)
    # Polishes the internal enclosure of *res* to at least *prec* bits
    # of precision in-place. Normally, *qqbar* operations that need
    # high-precision enclosures compute them on the fly without caching the results;
    # if *res* will be used as an invariant operand for many operations,
    # calling this function as a precomputation step can improve performance.

    void qqbar_denominator(fmpz_t res, const qqbar_t y)
    # Sets *res* to the denominator of *y*, i.e. the leading coefficient
    # of the minimal polynomial of *y*.

    void qqbar_numerator(qqbar_t res, const qqbar_t y)
    # Sets *res* to the numerator of *y*, i.e. *y* multiplied by
    # its denominator.

    void qqbar_conjugates(qqbar_ptr res, const qqbar_t x)
    # Sets the entries of the vector *res* to the *d* algebraic conjugates of
    # *x*, including *x* itself, where *d* is the degree of *x*. The output
    # is sorted in a canonical order (as defined by :func:`qqbar_cmp_root_order`).

    void _qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpz * poly, const fmpz_t den, slong len, const qqbar_t x)
    void qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpq_poly_t poly, const qqbar_t x)
    void _qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz * poly, slong len, const qqbar_t x)
    void qqbar_evaluate_fmpz_poly(qqbar_t res, const fmpz_poly_t poly, const qqbar_t x)
    # Sets *res* to the value of the given polynomial *poly* evaluated at
    # the algebraic number *x*. These methods detect simple special cases and
    # automatically reduce *poly* if its degree is greater or equal
    # to that of the minimal polynomial of *x*. In the generic case, evaluation
    # is done by computing minimal polynomials of representation matrices.

    int qqbar_evaluate_fmpz_mpoly_iter(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
    int qqbar_evaluate_fmpz_mpoly_horner(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
    int qqbar_evaluate_fmpz_mpoly(qqbar_t res, const fmpz_mpoly_t poly, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
    # Sets *res* to the value of *poly* evaluated at the algebraic numbers
    # given in the vector *x*. The number of variables is defined by
    # the context object *ctx*.
    # The parameters *deg_limit* and *bits_limit*
    # define evaluation limits: if any temporary result exceeds these limits
    # (not necessarily the final value, in case of cancellation), the
    # evaluation is aborted and 0 (failure) is returned. If evaluation
    # succeeds, 1 is returned.
    # The *iter* version iterates over all terms in succession and computes
    # the powers that appear. The *horner* version uses a multivariate
    # implementation of the Horner scheme. The default algorithm currently
    # uses the Horner scheme.

    void qqbar_roots_fmpz_poly(qqbar_ptr res, const fmpz_poly_t poly, int flags)
    void qqbar_roots_fmpq_poly(qqbar_ptr res, const fmpq_poly_t poly, int flags)
    # Sets the entries of the vector *res* to the *d* roots of the polynomial
    # *poly*. Roots with multiplicity appear with repetition in the
    # output array. By default, the roots will be sorted in a
    # convenient canonical order (as defined by :func:`qqbar_cmp_root_order`).
    # Instances of a repeated root always appear consecutively.
    # The following *flags* are supported:
    # - QQBAR_ROOTS_IRREDUCIBLE - if set, *poly* is assumed to be
    # irreducible (it may still have constant content), and no polynomial
    # factorization is performed internally.
    # - QQBAR_ROOTS_UNSORTED - if set, the roots will not be guaranteed
    # to be sorted (except for repeated roots being listed
    # consecutively).

    void qqbar_eigenvalues_fmpz_mat(qqbar_ptr res, const fmpz_mat_t mat, int flags)
    void qqbar_eigenvalues_fmpq_mat(qqbar_ptr res, const fmpq_mat_t mat, int flags)
    # Sets the entries of the vector *res* to the eigenvalues of the
    # square matrix *mat*. These functions compute the characteristic polynomial
    # of *mat* and then call :func:`qqbar_roots_fmpz_poly` with the same
    # flags.

    void qqbar_root_of_unity(qqbar_t res, slong p, ulong q)
    # Sets *res* to the root of unity `e^{2 \pi i p / q}`.

    int qqbar_is_root_of_unity(slong * p, ulong * q, const qqbar_t x)
    # If *x* is not a root of unity, returns 0.
    # If *x* is a root of unity, returns 1.
    # If *p* and *q* are not *NULL* and *x* is a root of unity,
    # this also sets *p* and *q* to the minimal integers with `0 \le p < q`
    # such that `x = e^{2 \pi i p / q}`.

    void qqbar_exp_pi_i(qqbar_t res, slong p, ulong q)
    # Sets *res* to the root of unity `e^{\pi i p / q}`.

    void qqbar_cos_pi(qqbar_t res, slong p, ulong q)
    void qqbar_sin_pi(qqbar_t res, slong p, ulong q)
    int qqbar_tan_pi(qqbar_t res, slong p, ulong q)
    int qqbar_cot_pi(qqbar_t res, slong p, ulong q)
    int qqbar_sec_pi(qqbar_t res, slong p, ulong q)
    int qqbar_csc_pi(qqbar_t res, slong p, ulong q)
    # Sets *res* to the trigonometric function `\cos(\pi x)`,
    # `\sin(\pi x)`, etc., with `x = \tfrac{p}{q}`.
    # The functions tan, cot, sec and csc return the flag 1 if the value exists,
    # and return 0 if the evaluation point is a pole of the function.

    int qqbar_log_pi_i(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{log}(x) / (\pi i)` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `-1 < y \le 1` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_atan_pi(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{atan}(x) / \pi` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `|y| < \tfrac{1}{2}` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_asin_pi(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{asin}(x) / \pi` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `|y| \le \tfrac{1}{2}` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_acos_pi(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{acos}(x) / \pi` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `0 \le y \le 1` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_acot_pi(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{acot}(x) / \pi` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `-\tfrac{1}{2} < y \le \tfrac{1}{2}` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_asec_pi(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{asec}(x) / \pi` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `0 \le y \le 1` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_acsc_pi(slong * p, ulong * q, const qqbar_t x)
    # If `y = \operatorname{acsc}(x) / \pi` is algebraic, and hence
    # necessarily rational, sets `y = p / q` to the reduced such
    # fraction with `-\tfrac{1}{2} \le y \le \tfrac{1}{2}` and returns 1.
    # If *y* is not algebraic, returns 0.

    int qqbar_guess(qqbar_t res, const acb_t z, slong max_deg, slong max_bits, int flags, slong prec)
    # Attempts to find an algebraic number *res* of degree at most *max_deg* and
    # height at most *max_bits* bits matching the numerical enclosure *z*.
    # The return flag indicates success.
    # This is only a heuristic method, and the return flag neither implies a
    # rigorous proof that *res* is the correct result, nor a rigorous proof
    # that no suitable algebraic number with the given *max_deg* and *max_bits*
    # exists. (Proof of nonexistence could in principle be computed,
    # but this is not yet implemented.)
    # The working precision *prec* should normally be the same as the precision
    # used to compute *z*. It does not make much sense to run this algorithm
    # with precision smaller than O(*max_deg* · *max_bits*).
    # This function does a single iteration at the target *max_deg*, *max_bits*,
    # and *prec*. For best performance, one should invoke this function
    # repeatedly with successively larger parameters when the size of the
    # intended solution is unknown or may be much smaller than a worst-case bound.

    int qqbar_express_in_field(fmpq_poly_t res, const qqbar_t alpha, const qqbar_t x, slong max_bits, int flags, slong prec)
    # Attempts to express *x* in the number field generated by *alpha*, returning
    # success (0 or 1). On success, *res* is set to a polynomial *f* of degree
    # less than the degree of *alpha* and with height (counting both the numerator
    # and the denominator, when the coefficients of *g* are
    # put on a common denominator) bounded by *max_bits* bits, such that
    # `f(\alpha) = x`.
    # (Exception: the *max_bits* parameter is currently ignored if *x* is
    # rational, in which case *res* is just set to the value of *x*.)
    # This function looks for a linear relation heuristically using a working
    # precision of *prec* bits. If *x* is expressible in terms of *alpha*,
    # then this function is guaranteed to succeed when *prec* is taken
    # large enough. The identity `f(\alpha) = x` is checked
    # rigorously, i.e. a return value of 1 implies a proof of correctness.
    # In principle, choosing a sufficiently large *prec* can be used to
    # prove that *x* does not lie in the field generated by *alpha*,
    # but the present implementation does not support doing so automatically.
    # This function does a single iteration at the target *max_bits* and
    # and *prec*. For best performance, one should invoke this function
    # repeatedly with successively larger parameters when the size of the
    # intended solution is unknown or may be much smaller than a worst-case bound.

    void qqbar_get_quadratic(fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t q, const qqbar_t x, int factoring)
    # Assuming that *x* has degree 1 or 2, computes integers *a*, *b*, *c*
    # and *q* such that
    # .. math ::
    # x = \frac{a + b \sqrt{c}}{q}
    # and such that *c* is not a perfect square, *q* is positive, and
    # *q* has no content in common with both *a* and *b*. In other words,
    # this determines a quadratic field `\mathbb{Q}(\sqrt{c})` containing
    # *x*, and then finds the canonical reduced coefficients *a*, *b* and
    # *q* expressing *x* in this field.
    # For convenience, this function supports rational *x*,
    # for which *b* and *c* will both be set to zero.
    # The following remarks apply to irrationals.
    # The radicand *c* will not be a perfect square, but will not
    # automatically be squarefree since this would require factoring the
    # discriminant. As a special case, *c* will be set to `-1` if *x*
    # is a Gaussian rational number. Otherwise, behavior is controlled
    # by the *factoring* parameter.
    # * If *factoring* is 0, no factorization is performed apart from
    # removing powers of two.
    # * If *factoring* is 1, a complete factorization is performed (*c*
    # will be minimal). This can be very expensive if the discriminant
    # is large.
    # * If *factoring* is 2, a smooth factorization is performed to remove
    # small factors from *c*. This is a tradeoff that provides pretty
    # output in most cases while avoiding extreme worst-case slowdown.
    # The smooth factorization guarantees finding all small factors
    # (up to some trial division limit determined internally by Flint),
    # but large factors are only found heuristically.

    int qqbar_set_fexpr(qqbar_t res, const fexpr_t expr)
    # Sets *res* to the algebraic number represented by the symbolic
    # expression *expr*, returning 1 on success and 0 on failure.
    # This function performs a "static" evaluation using *qqbar* arithmetic,
    # supporting only closed-form expressions with explicitly algebraic
    # subexpressions. It can be used to recover values generated by
    # :func:`qqbar_get_expr_formula` and variants.
    # For evaluating more complex expressions involving other
    # types of values or requiring symbolic simplifications,
    # the user should preprocess *expr* so that it is in a form
    # which can be parsed by :func:`qqbar_set_fexpr`.
    # The following expressions are supported:
    # * Integer constants
    # * Arithmetic operations with algebraic operands
    # * Square roots of algebraic numbers
    # * Powers with algebraic base and exponent an explicit rational number
    # * NumberI, GoldenRatio, RootOfUnity
    # * Floor, Ceil, Abs, Sign, Csgn, Conjugate, Re, Im, Max, Min
    # * Trigonometric functions with argument an explicit rational number times Pi
    # * Exponentials with argument an explicit rational number times Pi * NumberI
    # * The Decimal() constructor
    # * AlgebraicNumberSerialized() (assuming valid data, which is not checked)
    # * PolynomialRootIndexed()
    # * PolynomialRootNearest()
    # Examples of formulas that are not supported, despite the value being
    # an algebraic number:
    # * ``Pi - Pi``                 (general transcendental simplifications are not performed)
    # * ``1 / Infinity``            (only numbers are handled)
    # * ``Sum(n, For(n, 1, 10))``   (only static evaluation is performed)

    void qqbar_get_fexpr_repr(fexpr_t res, const qqbar_t x)
    # Sets *res* to a symbolic expression reflecting the exact internal
    # representation of *x*. The output will have the form
    # ``AlgebraicNumberSerialized(List(coeffs), enclosure)``.
    # The output can be converted back to a ``qqbar_t`` value using
    # :func:`qqbar_set_fexpr`. This is the recommended format for
    # serializing algebraic numbers as it requires minimal computation,
    # but it has the disadvantage of not being human-readable.

    void qqbar_get_fexpr_root_nearest(fexpr_t res, const qqbar_t x)
    # Sets *res* to a symbolic expression unambiguously describing *x*
    # in the form ``PolynomialRootNearest(List(coeffs), point)``
    # where *point* is an approximation of *x* guaranteed to be closer
    # to *x* than any conjugate root.
    # The output can be converted back to a ``qqbar_t`` value using
    # :func:`qqbar_set_fexpr`. This is a useful format for human-readable
    # presentation, but serialization and deserialization can be expensive.

    void qqbar_get_fexpr_root_indexed(fexpr_t res, const qqbar_t x)
    # Sets *res* to a symbolic expression unambiguously describing *x*
    # in the form ``PolynomialRootIndexed(List(coeffs), index)``
    # where *index* is the index of *x* among its conjugate roots
    # in the builtin root sort order.
    # The output can be converted back to a ``qqbar_t`` value using
    # :func:`qqbar_set_fexpr`. This is a useful format for human-readable
    # presentation when the numerical value is important, but serialization
    # and deserialization can be expensive.

    int qqbar_get_fexpr_formula(fexpr_t res, const qqbar_t x, ulong flags)
    # Attempts to express the algebraic number *x* as a closed-form expression
    # using arithmetic operations, radicals, and possibly exponentials
    # or trigonometric functions, but without using ``PolynomialRootNearest``
    # or ``PolynomialRootIndexed``.
    # Returns 0 on failure and 1 on success.
    # The *flags* parameter toggles different methods for generating formulas.
    # It can be set to any combination of the following. If *flags* is 0,
    # only rational numbers will be handled.
    # .. macro:: QQBAR_FORMULA_ALL
    # Toggles all methods (potentially expensive).
    # .. macro:: QQBAR_FORMULA_GAUSSIANS
    # Detect Gaussian rational numbers `a + bi`.
    # .. macro:: QQBAR_FORMULA_QUADRATICS
    # Solve quadratics in the form `a + b \sqrt{d}`.
    # .. macro:: QQBAR_FORMULA_CYCLOTOMICS
    # Detect elements of cyclotomic fields. This works by trying plausible
    # cyclotomic fields (based on the degree of the input), using LLL
    # to find candidate number field elements, and certifying candidates
    # through an exact computation. Detection is heuristic and
    # is not guaranteed to find all cyclotomic numbers.
    # .. macro:: QQBAR_FORMULA_CUBICS
    # QQBAR_FORMULA_QUARTICS
    # QQBAR_FORMULA_QUINTICS
    # Solve polynomials of degree 3, 4 and (where applicable) 5 using
    # cubic, quartic and quintic formulas (not yet implemented).
    # .. macro:: QQBAR_FORMULA_DEPRESSION
    # Use depression to try to generate simpler numbers.
    # .. macro:: QQBAR_FORMULA_DEFLATION
    # Use deflation to try to generate simpler numbers.
    # This allows handling number of the form `a^{1/n}` where *a* can
    # be represented in closed form.
    # .. macro:: QQBAR_FORMULA_SEPARATION
    # Try separating real and imaginary parts or sign and magnitude of
    # complex numbers. This allows handling numbers of the form `a + bi`
    # or `m \cdot s` (with `m > 0, |s| = 1`) where *a* and *b* or *m* and *s* can be
    # represented in closed form. This is only attempted as a fallback after
    # other methods fail: if an explicit Cartesian or magnitude-sign
    # represented is desired, the user should manually separate the number
    # into complex parts before calling :func:`qqbar_get_fexpr_formula`.
    # .. macro:: QQBAR_FORMULA_EXP_FORM
    # QQBAR_FORMULA_TRIG_FORM
    # QQBAR_FORMULA_RADICAL_FORM
    # QQBAR_FORMULA_AUTO_FORM
    # Select output form for cyclotomic numbers. The *auto* form (equivalent
    # to no flags being set) results in radicals for numbers of low degree,
    # trigonometric functions for real numbers, and complex exponentials
    # for nonreal numbers. The other flags (not fully implemented) can be
    # used to force exponential form, trigonometric form, or radical form.

    void qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op)
    # Given nonconstant polynomials *A* and *B*, sets *res* to a polynomial
    # whose roots are `a+b`, `a-b`, `ab` or `a/b` for all roots *a* of *A*
    # and all roots *b* of *B*. The parameter *op* selects the arithmetic
    # operation: 0 for addition, 1 for subtraction, 2 for multiplication
    # and 3 for division. If *op* is 3, *B* must not have zero as a root.

    void qqbar_binary_op(qqbar_t res, const qqbar_t x, const qqbar_t y, int op)
    # Performs a binary operation using a generic algorithm. This does not
    # check for special cases.

    int _qqbar_validate_uniqueness(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)
    # Given *z* known to be an enclosure of at least one root of *poly*,
    # certifies that the enclosure contains a unique root, and in that
    # case sets *res* to a new (possibly improved) enclosure for the same
    # root, returning 1. Returns 0 if uniqueness cannot be certified.
    # The enclosure is validated by performing a single step with the
    # interval Newton method. The working precision is determined from the
    # accuracy of *z*, but limited by *max_prec* bits.
    # This method slightly inflates the enclosure *z* to improve the chances
    # that the interval Newton step will succeed. Uniqueness on this larger
    # interval implies uniqueness of the original interval, but not
    # existence; when existence has not been ensured a priori,
    # :func:`_qqbar_validate_existence_uniqueness` should be used instead.

    int _qqbar_validate_existence_uniqueness(acb_t res, const fmpz_poly_t poly, const acb_t z, slong max_prec)
    # Given any complex interval *z*, certifies that the enclosure contains a
    # unique root of *poly*, and in that case sets *res* to a new (possibly
    # improved) enclosure for the same root, returning 1. Returns 0 if
    # existence and uniqueness cannot be certified.
    # The enclosure is validated by performing a single step with the
    # interval Newton method. The working precision is determined from the
    # accuracy of *z*, but limited by *max_prec* bits.

    void _qqbar_enclosure_raw(acb_t res, const fmpz_poly_t poly, const acb_t z, slong prec)
    void qqbar_enclosure_raw(acb_t res, const qqbar_t x, slong prec)
    # Sets *res* to an enclosure of *x* accurate to about *prec* bits
    # (the actual accuracy can be slightly lower, or higher).
    # This function uses repeated interval Newton steps to polish the initial
    # enclosure *z*, doubling the working precision each time. If any step
    # fails to improve the accuracy significantly, the root is recomputed
    # from scratch to higher precision.
    # If the initial enclosure is accurate enough, *res* is set to this value
    # without rounding and without further computation.

    int _qqbar_acb_lindep(fmpz * rel, acb_srcptr vec, slong len, int check, slong prec)
    # Attempts to find an integer vector *rel* giving a linear relation between
    # the elements of the real or complex vector *vec*, using the LLL algorithm.
    # The working precision is set to the minimum of *prec* and the relative
    # accuracy of *vec* (that is, the difference between the largest magnitude
    # and the largest error magnitude within *vec*). 95% of the bits within the
    # working precision are used for the LLL matrix, and the remaining 5% bits
    # are used to validate the linear relation by evaluating the linear
    # combination and checking that the resulting interval contains zero.
    # This validation does not prove the existence or nonexistence
    # of a linear relation, but it provides a quick heuristic way to eliminate
    # spurious relations.
    # If *check* is set, the return value indicates whether the validation
    # was successful; otherwise, the return value simply indicates whether
    # the algorithm was executed normally (failure may occur, for example,
    # if the input vector is non-finite).
    # In principle, this method can be used to produce a proof that no linear
    # relation exists with coefficients up to a specified bit size, but this has
    # not yet been implemented.
