# distutils: libraries = flint
# distutils: depends = flint/arb_fpwrap.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    int arb_fpwrap_double_exp(double * res, double x, int flags)
    int arb_fpwrap_cdouble_exp(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_expm1(double * res, double x, int flags)
    int arb_fpwrap_cdouble_expm1(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_log(double * res, double x, int flags)
    int arb_fpwrap_cdouble_log(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_log1p(double * res, double x, int flags)
    int arb_fpwrap_cdouble_log1p(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_pow(double * res, double x, double y, int flags)
    int arb_fpwrap_cdouble_pow(complex_double * res, complex_double x, complex_double y, int flags)

    int arb_fpwrap_double_sqrt(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sqrt(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_rsqrt(double * res, double x, int flags)
    int arb_fpwrap_cdouble_rsqrt(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cbrt(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cbrt(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sin(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sin(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cos(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cos(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_tan(double * res, double x, int flags)
    int arb_fpwrap_cdouble_tan(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cot(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cot(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sec(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sec(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_csc(double * res, double x, int flags)
    int arb_fpwrap_cdouble_csc(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sinc(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sinc(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sin_pi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sin_pi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cos_pi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cos_pi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_tan_pi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_tan_pi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cot_pi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cot_pi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sinc_pi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sinc_pi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_asin(double * res, double x, int flags)
    int arb_fpwrap_cdouble_asin(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_acos(double * res, double x, int flags)
    int arb_fpwrap_cdouble_acos(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_atan(double * res, double x, int flags)
    int arb_fpwrap_cdouble_atan(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_atan2(double * res, double x1, double x2, int flags)

    int arb_fpwrap_double_asinh(double * res, double x, int flags)
    int arb_fpwrap_cdouble_asinh(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_acosh(double * res, double x, int flags)
    int arb_fpwrap_cdouble_acosh(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_atanh(double * res, double x, int flags)
    int arb_fpwrap_cdouble_atanh(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_lambertw(double * res, double x, long branch, int flags)
    int arb_fpwrap_cdouble_lambertw(complex_double * res, complex_double x, long branch, int flags)

    int arb_fpwrap_double_rising(double * res, double x, double n, int flags)
    int arb_fpwrap_cdouble_rising(complex_double * res, complex_double x, complex_double n, int flags)
    # Rising factorial.

    int arb_fpwrap_double_gamma(double * res, double x, int flags)
    int arb_fpwrap_cdouble_gamma(complex_double * res, complex_double x, int flags)
    # Gamma function.

    int arb_fpwrap_double_rgamma(double * res, double x, int flags)
    int arb_fpwrap_cdouble_rgamma(complex_double * res, complex_double x, int flags)
    # Reciprocal gamma function.

    int arb_fpwrap_double_lgamma(double * res, double x, int flags)
    int arb_fpwrap_cdouble_lgamma(complex_double * res, complex_double x, int flags)
    # Log-gamma function.

    int arb_fpwrap_double_digamma(double * res, double x, int flags)
    int arb_fpwrap_cdouble_digamma(complex_double * res, complex_double x, int flags)
    # Digamma function.

    int arb_fpwrap_double_zeta(double * res, double x, int flags)
    int arb_fpwrap_cdouble_zeta(complex_double * res, complex_double x, int flags)
    # Riemann zeta function.

    int arb_fpwrap_double_hurwitz_zeta(double * res, double s, double z, int flags)
    int arb_fpwrap_cdouble_hurwitz_zeta(complex_double * res, complex_double s, complex_double z, int flags)
    # Hurwitz zeta function.

    int arb_fpwrap_double_lerch_phi(double * res, double z, double s, double a, int flags)
    int arb_fpwrap_cdouble_lerch_phi(complex_double * res, complex_double z, complex_double s, complex_double a, int flags)
    # Lerch transcendent.

    int arb_fpwrap_double_barnes_g(double * res, double x, int flags)
    int arb_fpwrap_cdouble_barnes_g(complex_double * res, complex_double x, int flags)
    # Barnes G-function.

    int arb_fpwrap_double_log_barnes_g(double * res, double x, int flags)
    int arb_fpwrap_cdouble_log_barnes_g(complex_double * res, complex_double x, int flags)
    # Logarithmic Barnes G-function.

    int arb_fpwrap_double_polygamma(double * res, double s, double z, int flags)
    int arb_fpwrap_cdouble_polygamma(complex_double * res, complex_double s, complex_double z, int flags)
    # Polygamma function.

    int arb_fpwrap_double_polylog(double * res, double s, double z, int flags)
    int arb_fpwrap_cdouble_polylog(complex_double * res, complex_double s, complex_double z, int flags)
    # Polylogarithm.

    int arb_fpwrap_cdouble_dirichlet_eta(complex_double * res, complex_double s, int flags)

    int arb_fpwrap_cdouble_riemann_xi(complex_double * res, complex_double s, int flags)

    int arb_fpwrap_cdouble_hardy_theta(complex_double * res, complex_double z, int flags)

    int arb_fpwrap_cdouble_hardy_z(complex_double * res, complex_double z, int flags)

    int arb_fpwrap_cdouble_zeta_zero(complex_double * res, unsigned long n, int flags)

    int arb_fpwrap_double_erf(double * res, double x, int flags)
    int arb_fpwrap_cdouble_erf(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_erfc(double * res, double x, int flags)
    int arb_fpwrap_cdouble_erfc(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_erfi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_erfi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_erfinv(double * res, double x, int flags)

    int arb_fpwrap_double_erfcinv(double * res, double x, int flags)

    int arb_fpwrap_double_fresnel_s(double * res, double x, int normalized, int flags)
    int arb_fpwrap_cdouble_fresnel_s(complex_double * res, complex_double x, int normalized, int flags)

    int arb_fpwrap_double_fresnel_c(double * res, double x, int normalized, int flags)
    int arb_fpwrap_cdouble_fresnel_c(complex_double * res, complex_double x, int normalized, int flags)

    int arb_fpwrap_double_gamma_upper(double * res, double s, double z, int regularized, int flags)
    int arb_fpwrap_cdouble_gamma_upper(complex_double * res, complex_double s, complex_double z, int regularized, int flags)

    int arb_fpwrap_double_gamma_lower(double * res, double s, double z, int regularized, int flags)
    int arb_fpwrap_cdouble_gamma_lower(complex_double * res, complex_double s, complex_double z, int regularized, int flags)

    int arb_fpwrap_double_beta_lower(double * res, double a, double b, double z, int regularized, int flags)
    int arb_fpwrap_cdouble_beta_lower(complex_double * res, complex_double a, complex_double b, complex_double z, int regularized, int flags)

    int arb_fpwrap_double_exp_integral_e(double * res, double s, double z, int flags)
    int arb_fpwrap_cdouble_exp_integral_e(complex_double * res, complex_double s, complex_double z, int flags)

    int arb_fpwrap_double_exp_integral_ei(double * res, double x, int flags)
    int arb_fpwrap_cdouble_exp_integral_ei(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sin_integral(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sin_integral(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cos_integral(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cos_integral(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_sinh_integral(double * res, double x, int flags)
    int arb_fpwrap_cdouble_sinh_integral(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_cosh_integral(double * res, double x, int flags)
    int arb_fpwrap_cdouble_cosh_integral(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_log_integral(double * res, double x, int offset, int flags)
    int arb_fpwrap_cdouble_log_integral(complex_double * res, complex_double x, int offset, int flags)

    int arb_fpwrap_double_dilog(double * res, double x, int flags)
    int arb_fpwrap_cdouble_dilog(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_bessel_j(double * res, double nu, double x, int flags)
    int arb_fpwrap_cdouble_bessel_j(complex_double * res, complex_double nu, complex_double x, int flags)

    int arb_fpwrap_double_bessel_y(double * res, double nu, double x, int flags)
    int arb_fpwrap_cdouble_bessel_y(complex_double * res, complex_double nu, complex_double x, int flags)

    int arb_fpwrap_double_bessel_i(double * res, double nu, double x, int flags)
    int arb_fpwrap_cdouble_bessel_i(complex_double * res, complex_double nu, complex_double x, int flags)

    int arb_fpwrap_double_bessel_k(double * res, double nu, double x, int flags)
    int arb_fpwrap_cdouble_bessel_k(complex_double * res, complex_double nu, complex_double x, int flags)

    int arb_fpwrap_double_bessel_k_scaled(double * res, double nu, double x, int flags)
    int arb_fpwrap_cdouble_bessel_k_scaled(complex_double * res, complex_double nu, complex_double x, int flags)

    int arb_fpwrap_double_airy_ai(double * res, double x, int flags)
    int arb_fpwrap_cdouble_airy_ai(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_airy_ai_prime(double * res, double x, int flags)
    int arb_fpwrap_cdouble_airy_ai_prime(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_airy_bi(double * res, double x, int flags)
    int arb_fpwrap_cdouble_airy_bi(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_airy_bi_prime(double * res, double x, int flags)
    int arb_fpwrap_cdouble_airy_bi_prime(complex_double * res, complex_double x, int flags)

    int arb_fpwrap_double_airy_ai_zero(double * res, unsigned long n, int flags)

    int arb_fpwrap_double_airy_ai_prime_zero(double * res, unsigned long n, int flags)

    int arb_fpwrap_double_airy_bi_zero(double * res, unsigned long n, int flags)

    int arb_fpwrap_double_airy_bi_prime_zero(double * res, unsigned long n, int flags)

    int arb_fpwrap_double_coulomb_f(double * res, double l, double eta, double x, int flags)
    int arb_fpwrap_cdouble_coulomb_f(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

    int arb_fpwrap_double_coulomb_g(double * res, double l, double eta, double x, int flags)
    int arb_fpwrap_cdouble_coulomb_g(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

    int arb_fpwrap_cdouble_coulomb_hpos(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)
    int arb_fpwrap_cdouble_coulomb_hneg(complex_double * res, complex_double l, complex_double eta, complex_double x, int flags)

    int arb_fpwrap_double_chebyshev_t(double * res, double n, double x, int flags)
    int arb_fpwrap_cdouble_chebyshev_t(complex_double * res, complex_double n, complex_double x, int flags)

    int arb_fpwrap_double_chebyshev_u(double * res, double n, double x, int flags)
    int arb_fpwrap_cdouble_chebyshev_u(complex_double * res, complex_double n, complex_double x, int flags)

    int arb_fpwrap_double_jacobi_p(double * res, double n, double a, double b, double x, int flags)
    int arb_fpwrap_cdouble_jacobi_p(complex_double * res, complex_double n, complex_double a, complex_double b, complex_double x, int flags)

    int arb_fpwrap_double_gegenbauer_c(double * res, double n, double m, double x, int flags)
    int arb_fpwrap_cdouble_gegenbauer_c(complex_double * res, complex_double n, complex_double m, complex_double x, int flags)

    int arb_fpwrap_double_laguerre_l(double * res, double n, double m, double x, int flags)
    int arb_fpwrap_cdouble_laguerre_l(complex_double * res, complex_double n, complex_double m, complex_double x, int flags)

    int arb_fpwrap_double_hermite_h(double * res, double n, double x, int flags)
    int arb_fpwrap_cdouble_hermite_h(complex_double * res, complex_double n, complex_double x, int flags)

    int arb_fpwrap_double_legendre_p(double * res, double n, double m, double x, int type, int flags)
    int arb_fpwrap_cdouble_legendre_p(complex_double * res, complex_double n, complex_double m, complex_double x, int type, int flags)

    int arb_fpwrap_double_legendre_q(double * res, double n, double m, double x, int type, int flags)
    int arb_fpwrap_cdouble_legendre_q(complex_double * res, complex_double n, complex_double m, complex_double x, int type, int flags)

    int arb_fpwrap_double_legendre_root(double * res1, double * res2, unsigned long n, unsigned long k, int flags)
    # Sets *res1* to the index *k* root of the Legendre polynomial `P_n(x)`,
    # and simultaneously sets *res2* to the corresponding weight for
    # Gauss-Legendre quadrature.

    int arb_fpwrap_cdouble_spherical_y(complex_double * res, long n, long m, complex_double x1, complex_double x2, int flags)

    int arb_fpwrap_double_hypgeom_0f1(double * res, double a, double x, int regularized, int flags)
    int arb_fpwrap_cdouble_hypgeom_0f1(complex_double * res, complex_double a, complex_double x, int regularized, int flags)

    int arb_fpwrap_double_hypgeom_1f1(double * res, double a, double b, double x, int regularized, int flags)
    int arb_fpwrap_cdouble_hypgeom_1f1(complex_double * res, complex_double a, complex_double b, complex_double x, int regularized, int flags)

    int arb_fpwrap_double_hypgeom_u(double * res, double a, double b, double x, int flags)
    int arb_fpwrap_cdouble_hypgeom_u(complex_double * res, complex_double a, complex_double b, complex_double x, int flags)

    int arb_fpwrap_double_hypgeom_2f1(double * res, double a, double b, double c, double x, int regularized, int flags)
    int arb_fpwrap_cdouble_hypgeom_2f1(complex_double * res, complex_double a, complex_double b, complex_double c, complex_double x, int regularized, int flags)

    int arb_fpwrap_double_hypgeom_pfq(double * res, const double * a, long p, const double * b, long q, double z, int regularized, int flags)
    int arb_fpwrap_cdouble_hypgeom_pfq(complex_double * res, const complex_double * a, long p, const complex_double * b, long q, complex_double z, int regularized, int flags)

    int arb_fpwrap_double_agm(double * res, double x, double y, int flags)
    int arb_fpwrap_cdouble_agm(complex_double * res, complex_double x, complex_double y, int flags)
    # Arithmetic-geometric mean.

    int arb_fpwrap_cdouble_elliptic_k(complex_double * res, complex_double m, int flags)

    int arb_fpwrap_cdouble_elliptic_e(complex_double * res, complex_double m, int flags)

    int arb_fpwrap_cdouble_elliptic_pi(complex_double * res, complex_double n, complex_double m, int flags)

    int arb_fpwrap_cdouble_elliptic_f(complex_double * res, complex_double phi, complex_double m, int pi, int flags)

    int arb_fpwrap_cdouble_elliptic_e_inc(complex_double * res, complex_double phi, complex_double m, int pi, int flags)

    int arb_fpwrap_cdouble_elliptic_pi_inc(complex_double * res, complex_double n, complex_double phi, complex_double m, int pi, int flags)
    # Complete and incomplete elliptic integrals.

    int arb_fpwrap_cdouble_elliptic_rf(complex_double * res, complex_double x, complex_double y, complex_double z, int option, int flags)

    int arb_fpwrap_cdouble_elliptic_rg(complex_double * res, complex_double x, complex_double y, complex_double z, int option, int flags)

    int arb_fpwrap_cdouble_elliptic_rj(complex_double * res, complex_double x, complex_double y, complex_double z, complex_double w, int option, int flags)
    # Carlson symmetric elliptic integrals.

    int arb_fpwrap_cdouble_elliptic_p(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_elliptic_p_prime(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_elliptic_inv_p(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_elliptic_zeta(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_elliptic_sigma(complex_double * res, complex_double z, complex_double tau, int flags)
    # Weierstrass elliptic functions.

    int arb_fpwrap_cdouble_jacobi_theta_1(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_jacobi_theta_2(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_jacobi_theta_3(complex_double * res, complex_double z, complex_double tau, int flags)

    int arb_fpwrap_cdouble_jacobi_theta_4(complex_double * res, complex_double z, complex_double tau, int flags)
    # Jacobi theta functions.

    int arb_fpwrap_cdouble_dedekind_eta(complex_double * res, complex_double tau, int flags)

    int arb_fpwrap_cdouble_modular_j(complex_double * res, complex_double tau, int flags)

    int arb_fpwrap_cdouble_modular_lambda(complex_double * res, complex_double tau, int flags)

    int arb_fpwrap_cdouble_modular_delta(complex_double * res, complex_double tau, int flags)
