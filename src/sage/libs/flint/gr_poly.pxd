# distutils: libraries = flint
# distutils: depends = flint/gr_poly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx)

    void gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx)

    void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx)

    gr_ptr gr_poly_entry_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)

    slong gr_poly_length(const gr_poly_t poly, gr_ctx_t ctx)

    void gr_poly_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)

    void gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

    void _gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx)

    void _gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)

    int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
    int gr_poly_get_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx)
    int gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx)
    int gr_poly_set_gr_poly_other(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_ctx, gr_ctx_t ctx)

    int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
    int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

    int gr_poly_truncate(gr_poly_t res, const gr_poly_t poly, slong newlen, gr_ctx_t ctx)

    int gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
    int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx)
    int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx)
    int gr_poly_gen(gr_poly_t poly, gr_ctx_t ctx)

    int gr_poly_write(gr_stream_t out, const gr_poly_t poly, const char * x, gr_ctx_t ctx)
    int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)

    int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)

    truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

    truth_t gr_poly_is_zero(const gr_poly_t poly, gr_ctx_t ctx)
    truth_t gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx)
    truth_t gr_poly_is_gen(const gr_poly_t poly, gr_ctx_t ctx)
    truth_t gr_poly_is_scalar(const gr_poly_t poly, gr_ctx_t ctx)

    int gr_poly_set_scalar(gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
    int gr_poly_set_si(gr_poly_t poly, slong c, gr_ctx_t ctx)
    int gr_poly_set_ui(gr_poly_t poly, ulong c, gr_ctx_t ctx)
    int gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
    int gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)

    int gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr c, gr_ctx_t ctx)
    int gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong c, gr_ctx_t ctx)
    int gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong c, gr_ctx_t ctx)
    int gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t c, gr_ctx_t ctx)
    int gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t c, gr_ctx_t ctx)

    int gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

    int gr_poly_neg(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

    int _gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

    int _gr_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

    int _gr_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)

    int _gr_poly_mullow_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
    int _gr_poly_mullow(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
    int gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong len, gr_ctx_t ctx)

    int gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)

    int _gr_poly_pow_series_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
    int gr_poly_pow_series_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

    int _gr_poly_pow_series_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
    int gr_poly_pow_series_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)

    int _gr_poly_pow_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
    int gr_poly_pow_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)

    int _gr_poly_pow_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
    int gr_poly_pow_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx)

    int gr_poly_pow_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx)

    int _gr_poly_pow_series_fmpq_recurrence(gr_ptr h, gr_srcptr f, slong flen, const fmpq_t exp, slong len, int precomp, gr_ctx_t ctx)
    int gr_poly_pow_series_fmpq_recurrence(gr_poly_t res, const gr_poly_t poly, const fmpq_t exp, slong len, gr_ctx_t ctx)

    int _gr_poly_shift_left(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
    int gr_poly_shift_left(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

    int _gr_poly_shift_right(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
    int gr_poly_shift_right(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)

    int _gr_poly_divrem_divconquer_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_divrem_divconquer_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
    int gr_poly_divrem_divconquer(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_divrem_basecase_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
    int _gr_poly_divrem_basecase_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int _gr_poly_divrem_basecase(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_divrem_basecase(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    # These functions implement Euclidean division with remainder:
    # given polynomials `A, B \in K[x]` where `K` is a field, with `B \ne 0`,
    # there is a unique quotient `Q` and remainder `R` such that `A = BQ + R`
    # and either `R = 0` or `\deg(R) < \deg(B)`.
    # If *B* is provably zero, ``GR_DOMAIN`` is returned.
    # When `K` is a commutative ring and `\operatorname{lc}(B)` is a unit in `K`,
    # the situation is the same as over fields. In particular, Euclidean division
    # with remainder always makes sense over commutative rings when `B` is monic.
    # If `\operatorname{lc}(B)` is not a unit, the division still makes sense if
    # the coefficient quotient `\operatorname{lc}(r)`  / `\operatorname{lc}(B)`
    # exists for each partial remainder `r`. Indeed,
    # the *basecase* and *divconquer* algorithms return ``GR_DOMAIN`` precisely when
    # encountering a leading quotient `\operatorname{lc}(r)`  / `\operatorname{lc}(B) \not \in K`.
    # However, the *newton* algorithm as currently implemented
    # returns ``GR_DOMAIN`` when `\operatorname{lc}(B)^{-1} \not \in K`.
    # The underscore methods make the following assumptions:
    # * *Q* has room for ``lenA - lenB + 1`` coefficients.
    # * *R* has room for ``lenB - 1`` coefficients.
    # * ``lenA >= lenB >= 1``.
    # * *Q* is not aliased with either *A* or *B*.
    # * *R* is not aliased with *B*.
    # * *R* may be aliased with *A*, in which case all ``lenA``
    # entries may be used as scratch space. Note that in this case,
    # only the low ``lenB - 1`` coefficients of *R* actually represent
    # valid coefficients on output: the higher scratch coefficients will not
    # necessarily be zeroed.
    # * The divisor *B* is normalized to have nonzero leading coefficient.
    # (The non-underscore methods check for leading coefficients that
    # are not provably nonzero and return ``GR_UNABLE``.)
    # The *preinv1* functions take a precomputed inverse of the
    # leading coefficient as input.
    # The *noinv* versions perform repeated checked divisions
    # by the leading coefficient.

    int _gr_poly_div_divconquer_preinv1(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_div_divconquer_noinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_div_divconquer(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
    int gr_poly_div_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_div_basecase_preinv1(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
    int _gr_poly_div_basecase_noinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int _gr_poly_div_basecase(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_div_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_div_newton(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_div_newton(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_div(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_div(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    # Versions of the *divrem* functions which output only the quotient.
    # These are generally faster.

    int _gr_poly_rem(gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_rem(gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    # Versions of the *divrem* functions which output only the remainder.

    int _gr_poly_inv_series_newton(gr_ptr res, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
    int gr_poly_inv_series_newton(gr_poly_t res, const gr_poly_t A, slong len, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_inv_series_basecase_preinv1(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr Ainv, slong len, gr_ctx_t ctx)
    int _gr_poly_inv_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
    int gr_poly_inv_series_basecase(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)
    int _gr_poly_inv_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
    int gr_poly_inv_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx)

    int _gr_poly_div_series_newton(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, slong cutoff, gr_ctx_t ctx)
    int gr_poly_div_series_newton(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_div_series_divconquer(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
    int gr_poly_div_series_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_div_series_invmul(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
    int gr_poly_div_series_invmul(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
    int _gr_poly_div_series_basecase_preinv1(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_srcptr Binv, slong len, gr_ctx_t ctx)
    int _gr_poly_div_series_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
    int _gr_poly_div_series_basecase(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
    int gr_poly_div_series_basecase(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
    int _gr_poly_div_series(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
    int gr_poly_div_series(gr_poly_t res, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)

    int _gr_poly_divexact_basecase_bidirectional(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
    int gr_poly_divexact_basecase_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_divexact_bidirectional(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
    int gr_poly_divexact_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_divexact_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
    int _gr_poly_divexact_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, gr_ctx_t ctx)
    int gr_poly_divexact_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    int _gr_poly_divexact_series_basecase_noinv(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
    int _gr_poly_divexact_series_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
    int gr_poly_divexact_series_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)

    int _gr_poly_sqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
    int gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_sqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_sqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_sqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_sqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_sqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_sqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

    int _gr_poly_rsqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
    int gr_poly_rsqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_rsqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_rsqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_rsqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_rsqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_rsqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_rsqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

    int _gr_poly_evaluate_rectangular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
    int gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

    int _gr_poly_evaluate_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
    int gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)

    int _gr_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx)
    int gr_poly_evaluate(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)
    # Set *res* to *poly* evaluated at *x*.

    int _gr_poly_evaluate_other_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    int gr_poly_evaluate_other_horner(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    int _gr_poly_evaluate_other_rectangular(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    int gr_poly_evaluate_other_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    int _gr_poly_evaluate_other(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    int gr_poly_evaluate_other(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
    # Set *res* to *poly* evaluated at *x*, where the coefficients of *f*
    # belong to *ctx* while both *x* and *res* belong to *x_ctx*.

    gr_ptr * _gr_poly_tree_alloc(slong len, gr_ctx_t ctx)

    void _gr_poly_tree_free(gr_ptr * tree, slong len, gr_ctx_t ctx)

    int _gr_poly_tree_build(gr_ptr * tree, gr_srcptr roots, slong len, gr_ctx_t ctx)

    int _gr_poly_evaluate_vec_fast_precomp(gr_ptr vs, gr_srcptr poly, slong plen, gr_ptr * tree, slong len, gr_ctx_t ctx)

    int _gr_poly_evaluate_vec_fast(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)

    int gr_poly_evaluate_vec_fast(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

    int _gr_poly_evaluate_vec_iter(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)

    int gr_poly_evaluate_vec_iter(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)

    int _gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
    int gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
    int _gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
    int gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
    int _gr_poly_taylor_shift_convolution(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
    int gr_poly_taylor_shift_convolution(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
    int _gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
    int gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
    # Sets *res* to the Taylor shift `f(x+c)`, where *f* is given by
    # *poly*, computed respectively using
    # an optimized form of Horner's rule, divide-and-conquer, a single
    # convolution, and an automatic choice between the three algorithms.
    # The underscore methods support aliasing.

    int _gr_poly_compose_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_compose_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
    int _gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_compose_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
    int _gr_poly_compose(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_compose(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
    # Sets *res* to the composition `f(g(x))` where *f* is given by *poly1*
    # and *g* is given by *poly2*, respectively using Horner's rule,
    # divide-and-conquer, and an automatic choice between the two algorithms.
    # The default algorithm also handles special-form input `g = ax^n + c`
    # efficiently by performing a Taylor shift followed by a rescaling.
    # The underscore methods do not support aliasing of the output
    # with either input polynomial.

    int _gr_poly_compose_series_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
    int gr_poly_compose_series_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
    int _gr_poly_compose_series_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
    int gr_poly_compose_series_brent_kung(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
    int _gr_poly_compose_series_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
    int gr_poly_compose_series_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
    int _gr_poly_compose_series(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
    int gr_poly_compose_series(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx)
    # Sets *res* to the power series composition `h(x) = f(g(x))` truncated
    # to order `O(x^n)` where `f` is given by *poly1* and `g` is given by *poly2*,
    # respectively using Horner's rule, the Brent-Kung baby step-giant step
    # algorithm [BrentKung1978]_, divide-and-conquer, and an automatic choice between the algorithms.
    # The default algorithm also handles short input and
    # special-form input `g = ax^n` efficiently.
    # We require that the constant term in `g(x)` is exactly zero.
    # The underscore methods do not support aliasing of the output
    # with either input polynomial, and do not zero-pad the result.

    int _gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
    int gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

    int _gr_poly_nth_derivative(gr_ptr res, gr_srcptr poly, ulong n, slong len, gr_ctx_t ctx)
    int gr_poly_nth_derivative(gr_poly_t res, const gr_poly_t poly, ulong n, gr_ctx_t ctx)

    int _gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
    int gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)

    int _gr_poly_make_monic(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
    int gr_poly_make_monic(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)

    truth_t _gr_poly_is_monic(gr_srcptr poly, slong len, gr_ctx_t ctx)
    truth_t gr_poly_is_monic(const gr_poly_t res, gr_ctx_t ctx)

    int _gr_poly_hgcd(gr_ptr r, slong * sgn, gr_ptr * M, slong * lenM, gr_ptr A, slong * lenA, gr_ptr B, slong * lenB, gr_srcptr a, slong lena, gr_srcptr b, slong lenb, slong cutoff, gr_ctx_t ctx)
    # Computes the HGCD of `a` and `b`, that is, a matrix `M`, a sign `\sigma`
    # and two polynomials `A` and `B` such that
    # .. math ::
    # (A,B)^t = \sigma M^{-1} (a,b)^t.
    # Assumes that `\operatorname{len}(a) > \operatorname{len}(b) > 0`.
    # Assumes that `A` and `B` have space of size at least `\operatorname{len}(a)`
    # and `\operatorname{len}(b)`, respectively.  On exit, ``*lenA`` and ``*lenB``
    # will contain the correct lengths of `A` and `B`.
    # Assumes that ``M[0]``, ``M[1]``, ``M[2]``, and ``M[3]``
    # each point to a vector of size at least `\operatorname{len}(a)`.
    # If `r` is not ``NULL``, writes to that variable the corresponding value
    # for computing resultants using the HGCD algorithm.

    int _gr_poly_gcd_hgcd(gr_ptr G, slong * _lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
    int gr_poly_gcd_hgcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_gcd_euclidean(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_gcd_euclidean(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    int _gr_poly_gcd(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_gcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
    # Polynomial GCD. Currently only useful over fields.
    # The underscore methods assume ``lenA >= lenB >= 1`` and that both
    # *A* and *B* have nonzero leading coefficient.
    # The underscore methods do not attempt to make the result monic.
    # The time complexity of the half-GCD algorithm is `\mathcal{O}(n \log^2 n)`
    # ring operations. For further details, see [ThullYap1990]_.

    int _gr_poly_xgcd_euclidean(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
    int gr_poly_xgcd_euclidean(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)

    int _gr_poly_xgcd_hgcd(slong * Glen, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)
    int gr_poly_xgcd_hgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)

    int _gr_poly_resultant_euclidean(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_resultant_euclidean(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
    int _gr_poly_resultant_hgcd(gr_ptr res, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
    int gr_poly_resultant_hgcd(gr_ptr res, const gr_poly_t f, const gr_poly_t g, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_resultant_sylvester(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_resultant_sylvester(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
    int _gr_poly_resultant_small(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_resultant_small(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
    int _gr_poly_resultant(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
    int gr_poly_resultant(gr_ptr res, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
    # Sets *res* to the resultant of *poly1* and *poly2*.
    # The underscore methods assume that `len1 \ge len2 \ge 1`
    # and that the leading coefficients are nonzero.
    # The *euclidean* algorithm is the ordinary Euclidean algorithm.
    # The *hgcd* version uses the quasilinear half-GCD algorithm.
    # It requires two extra tuning parameters ``inner_cutoff``
    # (recursion threshold passed forward to the HGCD algorithm)
    # and ``cutoff``. Both algorithms can fail when run over
    # non-fields; they will return ``GR_DOMAIN``
    # when encountering an impossible inverse.
    # The *small* version uses division-free straight-line programs
    # optimized for short polynomials.
    # It returns ``GR_UNABLE`` if the polynomials are too large.
    # Currently this function handles the cases where `len1 \le 2`
    # or `len2 \le 3`.
    # The *sylvester* version constructs the Sylvester matrix
    # and computes its determinant. This is useful over inexact rings
    # and as a fallback for rings without division.
    # The default version attempts to choose an appropriate
    # algorithm automatically.
    # Currently no algorithm has been implemented that is appropriate for
    # integral domains.

    int gr_poly_factor_squarefree(gr_ptr c, gr_vec_t fac, gr_vec_t exp, const gr_poly_t poly, gr_ctx_t ctx)
    # Computes a squarefree factorization of *poly*.
    # The constant *c* is set to an element of the scalar ring.
    # The factors in *fac* are set to polynomials; the user must thus
    # initialize it to a vector of polynomials of the same type as
    # *poly* (and *not* to the parent *ctx*).
    # The exponent vector *exp* must be initialized to the *fmpz* type.

    int gr_poly_squarefree_part(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
    # Sets *res* to the squarefreepart of *poly*.

    int gr_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
    int gr_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t poly_ctx, int flags, gr_ctx_t ctx)
    # Finds all roots of the given polynomial in the ring defined by *ctx*,
    # storing the roots without duplication in *roots* (a vector with
    # elements of type ``ctx``) and the corresponding multiplicities in
    # *mult* (a vector with elements of type ``fmpz``).
    # If the target ring is not an algebraically closed field, then
    # the sum of multiplicities can be smaller than the degree of the
    # polynomial. For example, with ``fmpz`` coefficients, we only
    # find integer roots.
    # The *other* version of this function takes as input a polynomial
    # with entries in a different ring ``poly_ctx``. For example,
    # we can compute ``qqbar`` or ``arb`` roots for a polynomial
    # with ``fmpz`` coefficients.
    # Whether the roots are sorted in any particular order is
    # ring-dependent.
    # We consider roots of the zero polynomial to be ill-defined and return
    # ``GR_DOMAIN`` in that case.

    int _gr_poly_asin_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_asin_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_asinh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_asinh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_acos_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_acos_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_acosh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_acosh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_atan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_atan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_atanh_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_atanh_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

    int _gr_poly_log_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_log_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
    int _gr_poly_log1p_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_log1p_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)

    int _gr_poly_exp_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
    int gr_poly_exp_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
    int _gr_poly_exp_series_basecase_mul(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
    int gr_poly_exp_series_basecase_mul(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
    int _gr_poly_exp_series_newton(gr_ptr f, gr_ptr g, gr_srcptr h, slong hlen, slong n, slong cutoff, gr_ctx_t ctx)
    int gr_poly_exp_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_exp_series_generic(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
    int _gr_poly_exp_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
    int gr_poly_exp_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)

    int _gr_poly_sin_cos_series_basecase(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, int times_pi, gr_ctx_t ctx)
    int gr_poly_sin_cos_series_basecase(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
    int _gr_poly_sin_cos_series_tangent(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, int times_pi, gr_ctx_t ctx)
    int gr_poly_sin_cos_series_tangent(gr_poly_t s, gr_poly_t c, const gr_poly_t h, slong n, int times_pi, gr_ctx_t ctx)
    # The *basecase* version uses a simple recurrence for the coefficients,
    # requiring `O(nm)` operations where `m` is the length of `h`.
    # The *tangent* version uses the tangent half-angle formulas to compute
    # the sine and cosine via :func:`_acb_poly_tan_series`. This
    # requires `O(M(n))` operations.
    # When `h = h_0 + h_1` where the constant term `h_0` is nonzero,
    # the evaluation is done as
    # `\sin(h_0 + h_1) = \cos(h_0) \sin(h_1) + \sin(h_0) \cos(h_1)`,
    # `\cos(h_0 + h_1) = \cos(h_0) \cos(h_1) - \sin(h_0) \sin(h_1)`.
    # The *basecase* and *tangent* versions take a flag *times_pi*
    # specifying that the input is to be multiplied by `\pi`.

    int _gr_poly_tan_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
    int gr_poly_tan_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
    int _gr_poly_tan_series_newton(gr_ptr f, gr_srcptr h, slong hlen, slong n, slong cutoff, gr_ctx_t ctx)
    int gr_poly_tan_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx)
    int _gr_poly_tan_series(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
    int gr_poly_tan_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
