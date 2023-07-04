# distutils: libraries = flint
# distutils: depends = flint/nmod_poly.h

from sage.libs.gmp.types cimport *
from sage.libs.flint.types cimport *

# flint/nmod_poly.h
cdef extern from "flint_wrap.h":
    # Memory management
    cdef void nmod_poly_init(nmod_poly_t poly, mp_limb_t n)
    cdef void nmod_poly_init_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv)
    cdef void nmod_poly_init2(nmod_poly_t poly, mp_limb_t n, long alloc)
    cdef void nmod_poly_init2_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv, long alloc)

    cdef void nmod_poly_clear(nmod_poly_t poly)

    cdef void nmod_poly_realloc(nmod_poly_t poly, long alloc)

    cdef void nmod_poly_fit_length(nmod_poly_t poly, long alloc)
    cdef void _nmod_poly_normalise(nmod_poly_t poly)

    # Getting and setting coefficients
    cdef unsigned long nmod_poly_get_coeff_ui(nmod_poly_t poly, unsigned long j)
    cdef void nmod_poly_set_coeff_ui(nmod_poly_t poly, unsigned long j, unsigned long c)

    # Input and output
    cdef char * nmod_poly_get_str(nmod_poly_t poly)
    cdef int nmod_poly_set_str(char * s, nmod_poly_t poly)
    cdef int nmod_poly_print(nmod_poly_t a)
    cdef int nmod_poly_fread(FILE * f, nmod_poly_t poly)
    cdef int nmod_poly_fprint(FILE * f, nmod_poly_t poly)
    cdef int nmod_poly_read(nmod_poly_t poly)

    # Polynomial parameters
    cdef long nmod_poly_length(nmod_poly_t poly)
    cdef long nmod_poly_degree(nmod_poly_t poly)
    cdef mp_limb_t nmod_poly_modulus(nmod_poly_t poly)
    cdef mp_bitcnt_t nmod_poly_max_bits(nmod_poly_t poly)

    # Assignment and basic manipulation
    cdef void nmod_poly_set(nmod_poly_t a, nmod_poly_t b)
    cdef void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_zero(nmod_poly_t res)
    cdef void nmod_poly_one(nmod_poly_t res)
    cdef void nmod_poly_truncate(nmod_poly_t poly, long len)
    cdef void nmod_poly_reverse(nmod_poly_t output, nmod_poly_t input, long m)
    cdef void nmod_poly_revert_series(nmod_poly_t output, nmod_poly_t intput, long m)

    cdef int nmod_poly_equal(nmod_poly_t a, nmod_poly_t b)

    # Powering
    cdef void nmod_poly_pow(nmod_poly_t res, nmod_poly_t poly, unsigned long e)
    cdef void nmod_poly_pow_trunc(nmod_poly_t res, nmod_poly_t poly, unsigned long e, long trunc)
    cdef void nmod_poly_powmod_ui_binexp(nmod_poly_t res, const nmod_poly_t poly, ulong e, const nmod_poly_t f)

    # Inflation and deflation
    cdef unsigned long nmod_poly_deflation(nmod_poly_t input)
    cdef void nmod_poly_deflate(nmod_poly_t result, nmod_poly_t input, unsigned long deflation)
    cdef void nmod_poly_inflate(nmod_poly_t result, nmod_poly_t input, unsigned long inflation)

    # Comparison
    cdef int nmod_poly_is_zero(nmod_poly_t poly)
    cdef int nmod_poly_is_one(nmod_poly_t poly)

    # Addition and subtraction
    cdef void nmod_poly_add(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_sub(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_neg(nmod_poly_t res, nmod_poly_t poly1)

    # Shifting
    cdef void nmod_poly_shift_left(nmod_poly_t res, nmod_poly_t poly, long k)
    cdef void nmod_poly_shift_right(nmod_poly_t res, nmod_poly_t poly, long k)

    # Multiplication
    cdef void nmod_poly_mul(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    cdef void nmod_poly_mullow(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, long trunc)
    cdef void nmod_poly_mulhigh(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, long n)
    cdef void nmod_poly_mulmod(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, nmod_poly_t f)

    # Square roots
    cdef void nmod_poly_invsqrt_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_sqrt_series(nmod_poly_t g, nmod_poly_t h, long n)
    int nmod_poly_sqrt(nmod_poly_t b, nmod_poly_t a)

    # Scalar multiplication and division
    cdef void nmod_poly_scalar_mul_nmod(nmod_poly_t res, nmod_poly_t poly1, mp_limb_t c)
    cdef void nmod_poly_make_monic(nmod_poly_t output, nmod_poly_t input)

    # Division
    cdef void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_div(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_inv_series(nmod_poly_t Qinv, nmod_poly_t Q, long n)
    cdef void nmod_poly_div_series(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B, long n)

    # GCD
    cdef void nmod_poly_gcd(nmod_poly_t G, nmod_poly_t A, nmod_poly_t B)
    cdef void nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, nmod_poly_t A, nmod_poly_t B)
    mp_limb_t nmod_poly_resultant(nmod_poly_t f, nmod_poly_t g)


    # Evaluation
    cdef mp_limb_t nmod_poly_evaluate_nmod(nmod_poly_t poly, mp_limb_t c)
    cdef void nmod_poly_evaluate_nmod_vec(mp_ptr ys, nmod_poly_t poly, mp_srcptr xs, long n)

    # Interpolation
    cdef void nmod_poly_interpolate_nmod_vec(nmod_poly_t poly, mp_srcptr xs, mp_srcptr ys, long n)

    # Composition
    cdef void nmod_poly_compose(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)

    # Power series composition and reversion
    cdef void nmod_poly_compose_series(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, long n)
    cdef void nmod_poly_reverse_series(nmod_poly_t Qinv, nmod_poly_t Q, long n)

    # Factoring
    cdef void nmod_poly_factor_clear(nmod_poly_factor_t fac)
    cdef void nmod_poly_factor_init(nmod_poly_factor_t fac)
    cdef void nmod_poly_factor_insert(nmod_poly_factor_t fac, nmod_poly_t poly, unsigned long exp)
    cdef void nmod_poly_factor_print(nmod_poly_factor_t fac)
    cdef void nmod_poly_factor_concat(nmod_poly_factor_t res, nmod_poly_factor_t fac)
    cdef void nmod_poly_factor_pow(nmod_poly_factor_t fac, unsigned long exp)
    cdef unsigned long nmod_poly_remove(nmod_poly_t f, nmod_poly_t p)
    cdef int nmod_poly_is_irreducible(nmod_poly_t f)
    cdef int nmod_poly_is_squarefree(nmod_poly_t f)
    cdef void nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, nmod_poly_t f)
    cdef void nmod_poly_factor_berlekamp(nmod_poly_factor_t factors, nmod_poly_t f)
    cdef void nmod_poly_factor_squarefree(nmod_poly_factor_t res, nmod_poly_t f)
    cdef mp_limb_t nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result, nmod_poly_t input)
    cdef mp_limb_t nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t result, nmod_poly_t input)
    cdef mp_limb_t nmod_poly_factor(nmod_poly_factor_t result, nmod_poly_t input)

    # Derivative
    cdef void nmod_poly_derivative(nmod_poly_t x_prime, nmod_poly_t x)
    cdef void nmod_poly_integral(nmod_poly_t x_int, nmod_poly_t x)

    # Transcendental functions
    cdef void nmod_poly_atan_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_tan_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_asin_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_sin_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_cos_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_asinh_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_atanh_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_sinh_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_cosh_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_tanh_series(nmod_poly_t g, nmod_poly_t h, long n)
    cdef void nmod_poly_log_series(nmod_poly_t res, nmod_poly_t f, long n)
    cdef void nmod_poly_exp_series(nmod_poly_t f, nmod_poly_t h, long)
