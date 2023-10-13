# distutils: libraries = flint
# distutils: depends = flint/gr_mpoly.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void gr_mpoly_init(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Initializes and sets *A* to the zero polynomial.

    void gr_mpoly_init3(gr_mpoly_t A, long alloc, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_init2(gr_mpoly_t A, long alloc, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Initializes *A* with space allocated for the given number
    # of coefficients and exponents with the given number of bits.

    void gr_mpoly_clear(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Clears *A*, freeing all allocated data.

    void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Swaps *A* and *B* efficiently.

    int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to *B*.

    int gr_mpoly_zero(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to the zero polynomial.

    truth_t gr_mpoly_is_zero(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Returns whether *A* is the zero polynomial.

    int gr_mpoly_gen(gr_mpoly_t A, long var, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to the generator with index *var* (indexed from zero).

    truth_t gr_mpoly_is_gen(const gr_mpoly_t A, long var, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Returns whether *A* is the generator with index *var* (indexed from zero).

    truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Returns whether *A* and *B* are equal.

    int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, long length, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to a random polynomial with up to *length* terms
    # and up to *exp_bits* bits in the exponents.

    int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, const char ** x, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_print_pretty(const gr_mpoly_t A, const char ** x, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Prints *A* using the strings in *x* for the variables.
    # If *x* is *NULL*, defaults are used.

    int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *c* to the coefficient in *A* with exponents *exp*.

    int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, unsigned long c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, long c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, unsigned long c, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, long c, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets the coefficient with exponents *exp* in *A* to the scalar *c*
    # which must be an element of or coercible to the coefficient ring.

    int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to the negation of *B*.

    int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to the difference of *B* and *C*.

    int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to the difference of *B* and *C*.

    int gr_mpoly_mul(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_johnson(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to the product of *B* and *C*.
    # The *monomial* version assumes that *C* is a monomial.

    int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, long c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, unsigned long c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Sets *A* to *B* multiplied by the scalar *c* which must be
    # an element of or coercible to the coefficient ring.

    void _gr_mpoly_fit_length(gr_ptr * coeffs, long * coeffs_alloc, unsigned long ** exps, long * exps_alloc, long N, long length, gr_ctx_t cctx)

    void gr_mpoly_fit_length(gr_mpoly_t A, long len, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    # Ensures that *A* has space for *len* coefficients and exponents.

    void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void gr_mpoly_fit_length_fit_bits(gr_mpoly_t A, long len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void gr_mpoly_fit_length_reset_bits(gr_mpoly_t A, long len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void _gr_mpoly_set_length(gr_mpoly_t A, long newlen, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const unsigned long * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void gr_mpoly_sort_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    int gr_mpoly_combine_like_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)

    void gr_mpoly_assert_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
