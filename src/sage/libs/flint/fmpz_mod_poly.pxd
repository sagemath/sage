# distutils: libraries = flint
# distutils: depends = flint/fmpz_mod_poly.h

from sage.libs.flint.types cimport *
from sage.libs.flint.thread_pool cimport *

# flint/fmpz_mod_poly.h
cdef extern from "flint_wrap.h":
    void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_t p, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_clear(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, slong alloc,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_truncate(fmpz_mod_poly_t poly, slong len, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_set_trunc(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx)

    slong fmpz_mod_poly_degree(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
    slong fmpz_mod_poly_length(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
    fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_poly_is_one(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_poly_is_gen(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set(fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1,
                               fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_reverse(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_zero(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_one(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly,
                                   slong i, slong j, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_ui(fmpz_mod_poly_t f, ulong x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f,
                                const fmpz_poly_t g, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f,
                            const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_equal_trunc(const fmpz_mod_poly_t poly1,
                const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n,
                                     const fmpz_t x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_coeff_ui(fmpz_mod_poly_t poly, slong n,
                                            ulong x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_coeff_si(fmpz_mod_poly_t poly, slong n,
                                            slong x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly,
                                             slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_set_coeff_mpz(fmpz_mod_poly_t poly,
                             slong n, const mpz_t x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_get_coeff_mpz(mpz_t x,
                 const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx)


    void _fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly,
                                                           slong len, slong n)

    void fmpz_mod_poly_shift_left(fmpz_mod_poly_t f,
                   const fmpz_mod_poly_t g, slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_shift_right(fmpz_mod_poly_t f,
                   const fmpz_mod_poly_t g, slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_add(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_add_series(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                            slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_sub(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_sub_series(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                            slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_neg(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res,
         const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_scalar_mul_ui(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, ulong x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_scalar_div_fmpz(fmpz_mod_poly_t res,
         const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_mul(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                        const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_mullow(fmpz_mod_poly_t res,
                    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                            slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_sqr(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_mulmod(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_mulmod_preinv(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                     const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op,
                                            ulong e, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res,
                            const fmpz_mod_poly_t poly, ulong e, slong trunc,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_pow_trunc_binexp(fmpz_mod_poly_t res,
                              const fmpz_mod_poly_t poly, ulong e, slong trunc,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_powmod_ui_binexp(fmpz_mod_poly_t res,
                            const fmpz_mod_poly_t poly, ulong e,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, ulong e,
                         const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_powmod_fmpz_binexp(fmpz_mod_poly_t res,
                            const fmpz_mod_poly_t poly, const fmpz_t e,
                            const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz_mod_poly_t res,
                          const fmpz_mod_poly_t poly, const fmpz_t e,
                          const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res,
         const fmpz_t e, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                     const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_powers_mod_naive(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_powers_mod_bsgs(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_frobenius_powers_2exp_precomp(
            fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_poly_t f,
                const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_frobenius_powers_2exp_clear(
          fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_frobenius_power(fmpz_mod_poly_t res,
                            fmpz_mod_poly_frobenius_powers_2exp_t pow,
                   const fmpz_mod_poly_t f, ulong m, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_frobenius_powers_precomp(
                fmpz_mod_poly_frobenius_powers_t pow, const fmpz_mod_poly_t f,
                const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_frobenius_powers_clear(
               fmpz_mod_poly_frobenius_powers_t pow, const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_div_basecase(fmpz_mod_poly_t Q,
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_div_newton_n_preinv(fmpz_mod_poly_t Q,
                         const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                         const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q,
          fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                         const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)

    ulong fmpz_mod_poly_remove(fmpz_mod_poly_t f,
                            const fmpz_mod_poly_t p, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R,
                            const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_divrem_divconquer(fmpz_mod_poly_t Q,
        fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_divrem_f(fmpz_t f, fmpz_mod_poly_t Q,
          fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_rem(fmpz_mod_poly_t R, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_rem_f(fmpz_t f, fmpz_mod_poly_t R, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_inv_series_newton(fmpz_mod_poly_t Qinv,
                   const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_inv_series_newton_f(fmpz_t f, fmpz_mod_poly_t Qinv,
                   const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_inv_series(fmpz_mod_poly_t Qinv, const fmpz_mod_poly_t Q,
                                             slong n, const fmpz_mod_ctx_t ctx)

    void  fmpz_mod_poly_inv_series_f(fmpz_t f, fmpz_mod_poly_t Qinv,
                    const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_div_series(fmpz_mod_poly_t Q,
                    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, slong n,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_make_monic_f(fmpz_t f, fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcd_euclidean(fmpz_mod_poly_t G,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcd_f(fmpz_t f, fmpz_mod_poly_t G, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcd_hgcd(fmpz_mod_poly_t G,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_gcd(fmpz_mod_poly_t G, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)


    void fmpz_mod_poly_xgcd_euclidean(fmpz_mod_poly_t G,
                             fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G,
                             fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_xgcd_hgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
          fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_xgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T,
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                      const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S,
          fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                      const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcdinv_euclidean_f(fmpz_t f, fmpz_mod_poly_t G,
          fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcdinv_euclidean(fmpz_mod_poly_t G,
          fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
                          const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz_mod_poly_t G,
          fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                                                     const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_invmod(fmpz_mod_poly_t A, const fmpz_mod_poly_t B,
                            const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_invmod_f(fmpz_t f, fmpz_mod_poly_t A,
                         const fmpz_mod_poly_t B, const fmpz_mod_poly_t P,
                                                     const fmpz_mod_ctx_t ctx)

    void  fmpz_mod_poly_minpoly_bm(fmpz_mod_poly_t poly, const fmpz* seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)
    void  fmpz_mod_poly_minpoly_hgcd(fmpz_mod_poly_t poly, const fmpz* seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)

    void  fmpz_mod_poly_minpoly(fmpz_mod_poly_t poly, const fmpz* seq, slong len,
                                                      const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_resultant_euclidean(fmpz_t r,
                            const fmpz_mod_poly_t f, const fmpz_mod_poly_t g,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_resultant_hgcd(fmpz_t res, const fmpz_mod_poly_t A,
                            const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

    void  fmpz_mod_poly_resultant(fmpz_t res, const fmpz_mod_poly_t f,
                             const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_discriminant(fmpz_t d, const fmpz_mod_poly_t f,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_derivative(fmpz_mod_poly_t res,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_evaluate_fmpz(fmpz_t res,
                                 const fmpz_mod_poly_t poly, const fmpz_t a,
                                                     const fmpz_mod_ctx_t ctx)

    fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(slong len)

    void fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys,
                        const fmpz_mod_poly_t poly, const fmpz * xs, slong n,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_horner(fmpz_mod_poly_t res,
                    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_divconquer(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                        const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                        const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_precompute_matrix(fmpz_mat_t A,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                     const fmpz_mod_poly_t poly2inv, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mat_t A,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res,
                   const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                   const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res,
                     const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                        const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(
      fmpz_mod_poly_struct * res, const fmpz_mod_poly_struct * polys,
      slong len1,slong n, const fmpz_mod_poly_t g, const fmpz_mod_poly_t poly,
                      const fmpz_mod_poly_t polyinv, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded_pool(fmpz_mod_poly_struct * res,
            const fmpz_mod_poly_struct * polys, slong len1, slong n,
            const fmpz_mod_poly_t g, const fmpz_mod_poly_t poly,
            const fmpz_mod_poly_t polyinv, const fmpz_mod_ctx_t ctx,
                              thread_pool_handle * threads, slong num_threads)

    void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv_threaded(fmpz_mod_poly_struct * res,
                    const fmpz_mod_poly_struct * polys, slong len1, slong n,
                    const fmpz_mod_poly_t g, const fmpz_mod_poly_t poly,
                      const fmpz_mod_poly_t polyinv, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_radix_init(fmpz_mod_poly_radix_t D,
                const fmpz_mod_poly_t R, slong degF, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_radix_clear(fmpz_mod_poly_radix_t D)

    void fmpz_mod_poly_radix(fmpz_mod_poly_struct **B, const fmpz_mod_poly_t F,
                      const fmpz_mod_poly_radix_t D, const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly,
                                                     const fmpz_mod_ctx_t ctx)

    int fmpz_mod_poly_fread(FILE * file, fmpz_mod_poly_t poly,
                                                           fmpz_mod_ctx_t ctx)


    int fmpz_mod_poly_fprint_pretty(FILE * file, const fmpz_mod_poly_t poly,
                                      const char * x, const fmpz_mod_ctx_t ctx)



    int fmpz_mod_poly_print(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)


    int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_product_roots_fmpz_vec(fmpz_poly_t poly, const fmpz * xs,
                                                    slong n, const fmpz_t mod)

    int fmpz_mod_poly_find_distinct_nonzero_roots(fmpz * roots,
                            const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_init(
                      fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_start_over(
                      fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_clear(fmpz_mod_berlekamp_massey_t B,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_print(
                const fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_add_points(
                   fmpz_mod_berlekamp_massey_t B, const fmpz * a, slong count,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_add_zeros(
                                   fmpz_mod_berlekamp_massey_t B, slong count,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_add_point(
                                fmpz_mod_berlekamp_massey_t B, const fmpz_t a,
                                                     const fmpz_mod_ctx_t ctx)

    void fmpz_mod_berlekamp_massey_add_point_ui(
                                      fmpz_mod_berlekamp_massey_t B, ulong a,
                                                     const fmpz_mod_ctx_t ctx)

    int fmpz_mod_berlekamp_massey_reduce(
                      fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

    const fmpz * fmpz_mod_berlekamp_massey_points(
                        const fmpz_mod_berlekamp_massey_t B)

    slong fmpz_mod_berlekamp_massey_point_count(
                        const fmpz_mod_berlekamp_massey_t B)

    const fmpz_mod_poly_struct * fmpz_mod_berlekamp_massey_V_poly(
                        const fmpz_mod_berlekamp_massey_t B)

    const fmpz_mod_poly_struct * fmpz_mod_berlekamp_massey_R_poly(
                        const fmpz_mod_berlekamp_massey_t B)

    void fmpz_mod_poly_add_si(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong c, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_sub_si(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong c, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_si_sub(fmpz_mod_poly_t res, slong c,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_add_fmpz(fmpz_mod_poly_t res,
               const fmpz_mod_poly_t poly, fmpz_t c, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_sub_fmpz(fmpz_mod_poly_t res,
               const fmpz_mod_poly_t poly, fmpz_t c, const fmpz_mod_ctx_t ctx)

    void fmpz_mod_poly_fmpz_sub(fmpz_mod_poly_t res, fmpz_t c,
                         const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
