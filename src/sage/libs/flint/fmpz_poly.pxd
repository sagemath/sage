from libc.stdio cimport FILE
from sage.libs.gmp.types cimport mpz_t
from sage.libs.flint.types cimport *

cdef extern from "flint/fmpz_poly.h":
    # Memory management
    void fmpz_poly_init(fmpz_poly_t)
    void fmpz_poly_init2(fmpz_poly_t, slong)

    void fmpz_poly_realloc(fmpz_poly_t, slong)

    void fmpz_poly_fit_length(fmpz_poly_t, slong)

    void fmpz_poly_clear(fmpz_poly_t)

    # Polynomial parameters
    slong fmpz_poly_length(const fmpz_poly_t)
    slong fmpz_poly_degree(const fmpz_poly_t)

    # Assignment and basic manipulation
    void fmpz_poly_set(fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_set_ui(fmpz_poly_t, ulong)
    void fmpz_poly_set_si(fmpz_poly_t, slong)
    void fmpz_poly_set_fmpz(fmpz_poly_t, const fmpz_t)
    void fmpz_poly_set_mpz(fmpz_poly_t, const mpz_t)
    int fmpz_poly_set_str(fmpz_poly_t, const char *)

    char *fmpz_poly_get_str(const fmpz_poly_t)
    char *fmpz_poly_get_str_pretty(const fmpz_poly_t, const char *)

    void fmpz_poly_zero(fmpz_poly_t)
    void fmpz_poly_one(fmpz_poly_t)

    void fmpz_poly_zero_coeffs(fmpz_poly_t, slong, slong)

    void fmpz_poly_reverse(fmpz_poly_t, const fmpz_poly_t, slong)

    void fmpz_poly_truncate(fmpz_poly_t, slong)

    # Getting and setting coefficients
    void fmpz_poly_get_coeff_fmpz(fmpz_t, const fmpz_poly_t, slong)
    slong fmpz_poly_get_coeff_si(const fmpz_poly_t, slong)
    ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t, slong)

    fmpz *fmpz_poly_get_coeff_ptr(const fmpz_poly_t, slong)
    fmpz *fmpz_poly_lead(const fmpz_poly_t)

    void fmpz_poly_set_coeff_fmpz(fmpz_poly_t, slong, const fmpz_t)
    void fmpz_poly_set_coeff_si(fmpz_poly_t, slong, slong)
    void fmpz_poly_set_coeff_ui(fmpz_poly_t, slong, ulong)

    # Comparison
    int fmpz_poly_equal(const fmpz_poly_t, const fmpz_poly_t)
    int fmpz_poly_is_zero(const fmpz_poly_t)
    int fmpz_poly_is_one(const fmpz_poly_t)
    int fmpz_poly_is_unit(const fmpz_poly_t)
    int fmpz_poly_equal_fmpz(const fmpz_poly_t, const fmpz_t)

    # Addition and subtraction
    void fmpz_poly_add(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_sub(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_neg(fmpz_poly_t, const fmpz_poly_t)

    # Scalar multiplication and division
    void fmpz_poly_scalar_mul_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_scalar_mul_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t)
    void fmpz_poly_scalar_mul_si(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_scalar_mul_ui(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_scalar_mul_2exp(fmpz_poly_t, const fmpz_poly_t, ulong)

    void fmpz_poly_scalar_addmul_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_scalar_submul_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)

    void fmpz_poly_scalar_fdiv_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_scalar_fdiv_si(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_scalar_fdiv_ui(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_scalar_fdiv_2exp(fmpz_poly_t, const fmpz_poly_t, ulong)

    void fmpz_poly_scalar_tdiv_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_scalar_tdiv_si(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_scalar_tdiv_2exp(fmpz_poly_t, const fmpz_poly_t, ulong)

    void fmpz_poly_scalar_divexact_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_scalar_divexact_si(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_scalar_divexact_ui(fmpz_poly_t, const fmpz_poly_t, ulong)

    void fmpz_poly_scalar_mod_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_scalar_smod_fmpz(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)

    # Multiplication
    void fmpz_poly_mul_classical(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_mullow_classical(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_mulhigh_classical(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_mulmid_classical(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_mul_karatsuba(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_mullow_karatsuba_n(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_mulhigh_karatsuba_n(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)

    void fmpz_poly_mul_KS(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_mullow_KS(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)

    void fmpz_poly_mul_SS(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_mullow_SS(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)

    void fmpz_poly_mul(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_mullow(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_mulhigh(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)

    # Squaring
    void fmpz_poly_sqr_classical(fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_sqrlow_classical(fmpz_poly_t, const fmpz_poly_t, slong)

    void fmpz_poly_sqr_karatsuba(fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_sqrlow_karatsuba_n(fmpz_poly_t, const fmpz_poly_t, slong)

    void fmpz_poly_sqr_KS(fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_sqrlow_KS(fmpz_poly_t, const fmpz_poly_t, slong)

    # Powering
    void fmpz_poly_pow_multinomial(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_pow_binomial(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_pow_addchains(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_pow_binexp(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_pow_small(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_pow(fmpz_poly_t, const fmpz_poly_t, ulong)
    void fmpz_poly_pow_trunc(fmpz_poly_t, const fmpz_poly_t, ulong, slong)

    # Shifting
    void fmpz_poly_shift_left(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_shift_right(fmpz_poly_t, const fmpz_poly_t, slong)

    # Bit sizes and norms
    ulong fmpz_poly_max_limbs(const fmpz_poly_t)
    slong fmpz_poly_max_bits(const fmpz_poly_t)
    void fmpz_poly_height(fmpz_t, const fmpz_poly_t)
    void fmpz_poly_2norm(fmpz_t, const fmpz_poly_t)

    # Greatest common divisor
    void fmpz_poly_gcd_subresultant(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    int fmpz_poly_gcd_heuristic(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_gcd_modular(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_gcd(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_xgcd_modular(
            fmpz_t, fmpz_poly_t, fmpz_poly_t,
                    const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_xgcd(
            fmpz_t, fmpz_poly_t, fmpz_poly_t,
                    const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_lcm(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_resultant(fmpz_t, const fmpz_poly_t, const fmpz_poly_t)

    # Gaussian content
    void fmpz_poly_content(fmpz_t, const fmpz_poly_t)
    void fmpz_poly_primitive_part(fmpz_poly_t, const fmpz_poly_t)

    # Square-free
    int fmpz_poly_is_squarefree(const fmpz_poly_t)

    # Euclidean division
    void fmpz_poly_divrem_basecase(
            fmpz_poly_t, fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_divrem_divconquer(
            fmpz_poly_t, fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_divrem(
            fmpz_poly_t, fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_div_basecase(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_div_divconquer(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_div(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_rem_basecase(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_rem(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_div_root(fmpz_poly_t, const fmpz_poly_t, const fmpz_t)

    # Division with precomputed inverse
    void fmpz_poly_preinvert(fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_div_preinv(
            fmpz_poly_t,
                    const fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_divrem_preinv(
            fmpz_poly_t, fmpz_poly_t,
                    const fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    # Power series division
    void fmpz_poly_inv_series_newton(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_inv_series(fmpz_poly_t, const fmpz_poly_t, slong)
    void fmpz_poly_div_series(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t, slong)

    # Pseudo division
    void fmpz_poly_pseudo_divrem_basecase(
            fmpz_poly_t, fmpz_poly_t,
                    ulong *, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_pseudo_divrem_divconquer(
            fmpz_poly_t, fmpz_poly_t,
                    ulong *, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_pseudo_divrem(
            fmpz_poly_t, fmpz_poly_t,
                    ulong *, const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_pseudo_div(
            fmpz_poly_t, ulong *, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_pseudo_rem(
            fmpz_poly_t, ulong *, const fmpz_poly_t, const fmpz_poly_t)

    void fmpz_poly_pseudo_divrem_cohen(
            fmpz_poly_t, fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_pseudo_rem_cohen(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    # Derivative
    void fmpz_poly_derivative(fmpz_poly_t, const fmpz_poly_t)

    # Evaluation
    void fmpz_poly_evaluate_divconquer_fmpz(
            fmpz_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_evaluate_horner_fmpz(
            fmpz_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_evaluate_fmpz(fmpz_t, const fmpz_poly_t, const fmpz_t)

    void fmpz_poly_evaluate_horner_fmpq(
            fmpq_t, const fmpz_poly_t, const fmpq_t)
    void fmpz_poly_evaluate_fmpq(fmpq_t, const fmpz_poly_t, const fmpq_t)

    mp_limb_t fmpz_poly_evaluate_mod(fmpz_poly_t, mp_limb_t, mp_limb_t)

    # Composition
    void fmpz_poly_compose_horner(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_compose_divconquer(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)
    void fmpz_poly_compose(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t)

    # Revert
    void fmpz_poly_revert_series_lagrange(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)
    void fmpz_poly_revert_series_lagrange_fast(fmpz_poly_t Qinv,
            fmpz_poly_t Q, long n)
    void fmpz_poly_revert_series_newton(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)
    void fmpz_poly_revert_series(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)

    #  Square root
    int fmpz_poly_sqrt_classical(fmpz_poly_t b, fmpz_poly_t a)
    int fmpz_poly_sqrt(fmpz_poly_t b, fmpz_poly_t a)

    #  Signature
    void fmpz_poly_signature(long * r1, long * r2, fmpz_poly_t poly)

    # Taylor shift
    void fmpz_poly_taylor_shift_horner(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_taylor_shift_divconquer(
            fmpz_poly_t, const fmpz_poly_t, const fmpz_t)
    void fmpz_poly_taylor_shift(fmpz_poly_t, const fmpz_poly_t, const fmpz_t)

    # Input and output
    int fmpz_poly_fprint(FILE *, fmpz_poly_t)
    int fmpz_poly_fprint_prett(FILE *, fmpz_poly_t, char *)
    int fmpz_poly_print(const fmpz_poly_t)
    int fmpz_poly_print_pretty(const fmpz_poly_t, const char *)

    int fmpz_poly_fread(FILE *, fmpz_poly_t)
    int fmpz_poly_fread_pretty(FILE *, fmpz_poly_t, char **)
    int fmpz_poly_read(fmpz_poly_t)
    int fmpz_poly_read_pretty(fmpz_poly_t, char **)

    # CRT
    void fmpz_poly_get_nmod_poly(nmod_poly_t, const fmpz_poly_t)

    void fmpz_poly_set_nmod_poly(fmpz_poly_t, const nmod_poly_t)
    void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t, const nmod_poly_t)

    void fmpz_poly_CRT_ui(
            fmpz_poly_t,
                    const fmpz_poly_t, const fmpz_t, const nmod_poly_t, int)

    # Some functions for backwards compatibility
    void fmpz_poly_scalar_mul_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t)
    void fmpz_poly_scalar_divexact_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t)
    void fmpz_poly_scalar_fdiv_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t)
    void fmpz_poly_set_coeff_mpz(fmpz_poly_t, slong, const mpz_t)
    void fmpz_poly_get_coeff_mpz(mpz_t, const fmpz_poly_t, slong)


# Wrapper Cython class
from sage.structure.sage_object cimport SageObject
cdef class Fmpz_poly(SageObject):
    cdef fmpz_poly_t poly
