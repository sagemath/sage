# distutils: libraries = flint
# distutils: depends = flint/arb_hypgeom.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void _arb_hypgeom_rising_coeffs_1(unsigned long * c, unsigned long k, long n)
    void _arb_hypgeom_rising_coeffs_2(unsigned long * c, unsigned long k, long n)
    void _arb_hypgeom_rising_coeffs_fmpz(fmpz * c, unsigned long k, long n)
    # Sets *c* to the coefficients of the rising factorial polynomial
    # `(X+k)_n`. The *1* and *2* versions respectively
    # compute single-word and double-word coefficients, without checking for
    # overflow, while the *fmpz* version allows arbitrarily large coefficients.
    # These functions are mostly intended for internal use; the *fmpz* version
    # does not use an asymptotically fast algorithm.
    # The degree *n* must be at least 2.

    void arb_hypgeom_rising_ui_forward(arb_t res, const arb_t x, unsigned long n, long prec)
    void arb_hypgeom_rising_ui_bs(arb_t res, const arb_t x, unsigned long n, long prec)
    void arb_hypgeom_rising_ui_rs(arb_t res, const arb_t x, unsigned long n, unsigned long m, long prec)
    void arb_hypgeom_rising_ui_rec(arb_t res, const arb_t x, unsigned long n, long prec)
    void arb_hypgeom_rising_ui(arb_t res, const arb_t x, unsigned long n, long prec)
    void arb_hypgeom_rising(arb_t res, const arb_t x, const arb_t n, long prec)
    # Computes the rising factorial `(x)_n`.
    # The *forward* version uses the forward recurrence.
    # The *bs* version uses binary splitting.
    # The *rs* version uses rectangular splitting. It takes an extra tuning
    # parameter *m* which can be set to zero to choose automatically.
    # The *rec* version chooses an algorithm automatically, avoiding
    # use of the gamma function (so that it can be used in the computation
    # of the gamma function).
    # The default versions (*rising_ui* and *rising_ui*) choose an algorithm
    # automatically and may additionally fall back on the gamma function.

    void arb_hypgeom_rising_ui_jet_powsum(arb_ptr res, const arb_t x, unsigned long n, long len, long prec)
    void arb_hypgeom_rising_ui_jet_bs(arb_ptr res, const arb_t x, unsigned long n, long len, long prec)
    void arb_hypgeom_rising_ui_jet_rs(arb_ptr res, const arb_t x, unsigned long n, unsigned long m, long len, long prec)
    void arb_hypgeom_rising_ui_jet(arb_ptr res, const arb_t x, unsigned long n, long len, long prec)
    # Computes the jet of the rising factorial `(x)_n`, truncated to length *len*.
    # In other words, constructs the polynomial `(X + x)_n \in \mathbb{R}[X]`,
    # truncated if `\operatorname{len} < n + 1` (and zero-extended
    # if `\operatorname{len} > n + 1`).
    # The *powsum* version computes the sequence of powers of *x* and forms integral
    # linear combinations of these.
    # The *bs* version uses binary splitting.
    # The *rs* version uses rectangular splitting. It takes an extra tuning
    # parameter *m* which can be set to zero to choose automatically.
    # The default version chooses an algorithm automatically.

    void _arb_hypgeom_gamma_stirling_term_bounds(long * bound, const mag_t zinv, long N)
    # For `1 \le n < N`, sets *bound* to an exponent bounding the *n*-th term
    # in the Stirling series for the gamma function, given a precomputed upper
    # bound for `|z|^{-1}`. This function is intended for internal use and
    # does not check for underflow or underflow in the exponents.

    void arb_hypgeom_gamma_stirling_sum_horner(arb_t res, const arb_t z, long N, long prec)
    void arb_hypgeom_gamma_stirling_sum_improved(arb_t res, const arb_t z, long N, long K, long prec)
    # Sets *res* to the final sum in the Stirling series for the gamma function
    # truncated before the term with index *N*, i.e. computes
    # `\sum_{n=1}^{N-1} B_{2n} / (2n(2n-1) z^{2n-1})`.
    # The *horner* version uses Horner scheme with gradual precision adjustments.
    # The *improved* version uses rectangular splitting for the low-index
    # terms and reexpands the high-index terms as hypergeometric polynomials,
    # using a splitting parameter *K* (which can be set to 0 to use a default
    # value).

    void arb_hypgeom_gamma_stirling(arb_t res, const arb_t x, int reciprocal, long prec)
    # Sets *res* to the gamma function of *x* computed using the Stirling
    # series together with argument reduction. If *reciprocal* is set,
    # the reciprocal gamma function is computed instead.

    int arb_hypgeom_gamma_taylor(arb_t res, const arb_t x, int reciprocal, long prec)
    # Attempts to compute the gamma function of *x* using Taylor series
    # together with argument reduction. This is only supported if *x* and *prec*
    # are both small enough. If successful, returns 1; otherwise, does nothing
    # and returns 0. If *reciprocal* is set, the reciprocal gamma function is
    # computed instead.

    void arb_hypgeom_gamma(arb_t res, const arb_t x, long prec)
    void arb_hypgeom_gamma_fmpq(arb_t res, const fmpq_t x, long prec)
    void arb_hypgeom_gamma_fmpz(arb_t res, const fmpz_t x, long prec)
    # Sets *res* to the gamma function of *x* computed using a default
    # algorithm choice.

    void arb_hypgeom_rgamma(arb_t res, const arb_t x, long prec)
    # Sets *res* to the reciprocal gamma function of *x* computed using a default
    # algorithm choice.

    void arb_hypgeom_lgamma(arb_t res, const arb_t x, long prec)
    # Sets *res* to the log-gamma function of *x* computed using a default
    # algorithm choice.

    void arb_hypgeom_central_bin_ui(arb_t res, unsigned long n, long prec)
    # Computes the central binomial coefficient `{2n \choose n}`.

    void arb_hypgeom_pfq(arb_t res, arb_srcptr a, long p, arb_srcptr b, long q, const arb_t z, int regularized, long prec)
    # Computes the generalized hypergeometric function `{}_pF_{q}(z)`,
    # or the regularized version if *regularized* is set.

    void arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, long prec)
    # Computes the confluent hypergeometric limit function
    # `{}_0F_1(a,z)`, or `\frac{1}{\Gamma(a)} {}_0F_1(a,z)` if *regularized*
    # is set.

    void arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    # Computes the confluent hypergeometric function
    # `M(a,b,z) = {}_1F_1(a,b,z)`, or
    # `\mathbf{M}(a,b,z) = \frac{1}{\Gamma(b)} {}_1F_1(a,b,z)` if *regularized* is set.

    void arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    # Alias for :func:`arb_hypgeom_m`.

    void arb_hypgeom_1f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    # Computes the confluent hypergeometric function using numerical integration
    # of the representation
    # .. math ::
    # {}_1F_1(a,b,z) = \frac{\Gamma(b)}{\Gamma(a) \Gamma(b-a)} \int_0^1 e^{zt} t^{a-1} (1-t)^{b-a-1} dt.
    # This algorithm can be useful if the parameters are large. This will currently
    # only return a finite enclosure if `a \ge 1` and `b - a \ge 1`.

    void arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, long prec)
    # Computes the confluent hypergeometric function `U(a,b,z)`.

    void arb_hypgeom_u_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, long prec)
    # Computes the confluent hypergeometric function `U(a,b,z)` using numerical integration
    # of the representation
    # .. math ::
    # U(a,b,z) = \frac{1}{\Gamma(a)} \int_0^{\infty} e^{-zt} t^{a-1} (1+t)^{b-a-1} dt.
    # This algorithm can be useful if the parameters are large. This will currently
    # only return a finite enclosure if `a \ge 1` and `z > 0`.

    void arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, long prec)
    # Computes the Gauss hypergeometric function
    # `{}_2F_1(a,b,c,z)`, or
    # `\mathbf{F}(a,b,c,z) = \frac{1}{\Gamma(c)} {}_2F_1(a,b,c,z)`
    # if *regularized* is set.
    # Additional evaluation flags can be passed via the *regularized*
    # argument; see :func:`acb_hypgeom_2f1` for documentation.

    void arb_hypgeom_2f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, long prec)
    # Computes the Gauss hypergeometric function using numerical integration
    # of the representation
    # .. math ::
    # {}_2F_1(a,b,c,z) = \frac{\Gamma(a)}{\Gamma(b) \Gamma(c-b)} \int_0^1 t^{b-1} (1-t)^{c-b-1} (1-zt)^{-a} dt.
    # This algorithm can be useful if the parameters are large. This will currently
    # only return a finite enclosure if `b \ge 1` and `c - b \ge 1` and
    # `z < 1`, possibly with *a* and *b* exchanged.

    void arb_hypgeom_erf(arb_t res, const arb_t z, long prec)
    # Computes the error function `\operatorname{erf}(z)`.

    void _arb_hypgeom_erf_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_erf_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the error function of the power series *z*,
    # truncated to length *len*.

    void arb_hypgeom_erfc(arb_t res, const arb_t z, long prec)
    # Computes the complementary error function
    # `\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)`.
    # This function avoids catastrophic cancellation for large positive *z*.

    void _arb_hypgeom_erfc_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_erfc_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the complementary error function of the power series *z*,
    # truncated to length *len*.

    void arb_hypgeom_erfi(arb_t res, const arb_t z, long prec)
    # Computes the imaginary error function
    # `\operatorname{erfi}(z) = -i\operatorname{erf}(iz)`.

    void _arb_hypgeom_erfi_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_erfi_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the imaginary error function of the power series *z*,
    # truncated to length *len*.

    void arb_hypgeom_erfinv(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erfcinv(arb_t res, const arb_t z, long prec)
    # Computes the inverse error function `\operatorname{erf}^{-1}(z)`
    # or inverse complementary error function `\operatorname{erfc}^{-1}(z)`.

    void arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, long prec)
    # Sets *res1* to the Fresnel sine integral `S(z)` and *res2* to
    # the Fresnel cosine integral `C(z)`. Optionally, just a single function
    # can be computed by passing *NULL* as the other output variable.
    # The definition `S(z) = \int_0^z \sin(t^2) dt` is used if *normalized* is 0,
    # and `S(z) = \int_0^z \sin(\tfrac{1}{2} \pi t^2) dt` is used if
    # *normalized* is 1 (the latter is the Abramowitz & Stegun convention).
    # `C(z)` is defined analogously.

    void _arb_hypgeom_fresnel_series(arb_ptr res1, arb_ptr res2, arb_srcptr z, long zlen, int normalized, long len, long prec)
    void arb_hypgeom_fresnel_series(arb_poly_t res1, arb_poly_t res2, const arb_poly_t z, int normalized, long len, long prec)
    # Sets *res1* to the Fresnel sine integral and *res2* to the Fresnel
    # cosine integral of the power series *z*, truncated to length *len*.
    # Optionally, just a single function can be computed by passing *NULL*
    # as the other output variable.

    void arb_hypgeom_gamma_upper(arb_t res, const arb_t s, const arb_t z, int regularized, long prec)
    # If *regularized* is 0, computes the upper incomplete gamma function
    # `\Gamma(s,z)`.
    # If *regularized* is 1, computes the regularized upper incomplete
    # gamma function `Q(s,z) = \Gamma(s,z) / \Gamma(s)`.
    # If *regularized* is 2, computes the generalized exponential integral
    # `z^{-s} \Gamma(s,z) = E_{1-s}(z)` instead (this option is mainly
    # intended for internal use; :func:`arb_hypgeom_expint` is the intended
    # interface for computing the exponential integral).

    void arb_hypgeom_gamma_upper_integration(arb_t res, const arb_t s, const arb_t z, int regularized, long prec)
    # Computes the upper incomplete gamma function using numerical
    # integration.

    void _arb_hypgeom_gamma_upper_series(arb_ptr res, const arb_t s, arb_srcptr z, long zlen, int regularized, long n, long prec)
    void arb_hypgeom_gamma_upper_series(arb_poly_t res, const arb_t s, const arb_poly_t z, int regularized, long n, long prec)
    # Sets *res* to an upper incomplete gamma function where *s* is
    # a constant and *z* is a power series, truncated to length *n*.
    # The *regularized* argument has the same interpretation as in
    # :func:`arb_hypgeom_gamma_upper`.

    void arb_hypgeom_gamma_lower(arb_t res, const arb_t s, const arb_t z, int regularized, long prec)
    # If *regularized* is 0, computes the lower incomplete gamma function
    # `\gamma(s,z) = \frac{z^s}{s} {}_1F_1(s, s+1, -z)`.
    # If *regularized* is 1, computes the regularized lower incomplete
    # gamma function `P(s,z) = \gamma(s,z) / \Gamma(s)`.
    # If *regularized* is 2, computes a further regularized lower incomplete
    # gamma function `\gamma^{*}(s,z) = z^{-s} P(s,z)`.

    void _arb_hypgeom_gamma_lower_series(arb_ptr res, const arb_t s, arb_srcptr z, long zlen, int regularized, long n, long prec)
    void arb_hypgeom_gamma_lower_series(arb_poly_t res, const arb_t s, const arb_poly_t z, int regularized, long n, long prec)
    # Sets *res* to an lower incomplete gamma function where *s* is
    # a constant and *z* is a power series, truncated to length *n*.
    # The *regularized* argument has the same interpretation as in
    # :func:`arb_hypgeom_gamma_lower`.

    void arb_hypgeom_beta_lower(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    # Computes the (lower) incomplete beta function, defined by
    # `B(a,b;z) = \int_0^z t^{a-1} (1-t)^{b-1}`,
    # optionally the regularized incomplete beta function
    # `I(a,b;z) = B(a,b;z) / B(a,b;1)`.

    void _arb_hypgeom_beta_lower_series(arb_ptr res, const arb_t a, const arb_t b, arb_srcptr z, long zlen, int regularized, long n, long prec)
    void arb_hypgeom_beta_lower_series(arb_poly_t res, const arb_t a, const arb_t b, const arb_poly_t z, int regularized, long n, long prec)
    # Sets *res* to the lower incomplete beta function `B(a,b;z)` (optionally
    # the regularized version `I(a,b;z)`) where *a* and *b* are constants
    # and *z* is a power series, truncating the result to length *n*.
    # The underscore method requires positive lengths and does not support
    # aliasing.

    void _arb_hypgeom_gamma_lower_sum_rs_1(arb_t res, unsigned long p, unsigned long q, const arb_t z, long N, long prec)
    # Computes `\sum_{k=0}^{N-1} z^k / (a)_k` where `a = p/q` using
    # rectangular splitting. It is assumed that `p + qN` fits in a limb.

    void _arb_hypgeom_gamma_upper_sum_rs_1(arb_t res, unsigned long p, unsigned long q, const arb_t z, long N, long prec)
    # Computes `\sum_{k=0}^{N-1} (a)_k / z^k` where `a = p/q` using
    # rectangular splitting. It is assumed that `p + qN` fits in a limb.

    long _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol)
    # Returns number of terms *N* and sets *err* to the truncation error for evaluating
    # `\Gamma(a,z)` using the asymptotic series at infinity, targeting an absolute
    # tolerance of *abs_tol*. The error may be set to *err* if the tolerance
    # cannot be achieved. Assumes that *z* is positive.

    void _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(arb_t res, const fmpq_t a, const arb_t z, long N, long prec)
    # Sets *res* to the approximation of `\Gamma(a,z)` obtained by truncating
    # the asymptotic series at infinity before term *N*.
    # The truncation error bound has to be added separately.

    long _arb_hypgeom_gamma_lower_fmpq_0_choose_N(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol)
    # Returns number of terms *N* and sets *err* to the truncation error for evaluating
    # `\gamma(a,z)` using the Taylor series at zero, targeting an absolute
    # tolerance of *abs_tol*. Assumes that *z* is positive.

    void _arb_hypgeom_gamma_lower_fmpq_0_bsplit(arb_t res, const fmpq_t a, const arb_t z, long N, long prec)
    # Sets *res* to the approximation of `\gamma(a,z)` obtained by truncating
    # the Taylor series at zero before term *N*.
    # The truncation error bound has to be added separately.

    long _arb_hypgeom_gamma_upper_singular_si_choose_N(mag_t err, long n, const arb_t z, const mag_t abs_tol)
    # Returns number of terms *N* and sets *err* to the truncation error for evaluating
    # `\Gamma(-n,z)` using the Taylor series at zero, targeting an absolute
    # tolerance of *abs_tol*.

    void _arb_hypgeom_gamma_upper_singular_si_bsplit(arb_t res, long n, const arb_t z, long N, long prec)
    # Sets *res* to the approximation of `\Gamma(-n,z)` obtained by truncating
    # the Taylor series at zero before term *N*.
    # The truncation error bound has to be added separately.

    void _arb_gamma_upper_fmpq_step_bsplit(arb_t Gz1, const fmpq_t a, const arb_t z0, const arb_t z1, const arb_t Gz0, const arb_t expmz0, const mag_t abs_tol, long prec)
    # Given *Gz0* and *expmz0* representing the values `\Gamma(a,z_0)` and `\exp(-z_0)`,
    # computes `\Gamma(a,z_1)` using the Taylor series at `z_0` evaluated
    # using binary splitting,
    # targeting an absolute error of *abs_tol*.
    # Assumes that `z_0` and `z_1` are positive.

    void arb_hypgeom_expint(arb_t res, const arb_t s, const arb_t z, long prec)
    # Computes the generalized exponential integral `E_s(z)`.

    void arb_hypgeom_ei(arb_t res, const arb_t z, long prec)
    # Computes the exponential integral `\operatorname{Ei}(z)`.

    void _arb_hypgeom_ei_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_ei_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the exponential integral of the power series *z*,
    # truncated to length *len*.

    void _arb_hypgeom_si_asymp(arb_t res, const arb_t z, long N, long prec)
    void _arb_hypgeom_si_1f2(arb_t res, const arb_t z, long N, long wp, long prec)
    void arb_hypgeom_si(arb_t res, const arb_t z, long prec)
    # Computes the sine integral `\operatorname{Si}(z)`.

    void _arb_hypgeom_si_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_si_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the sine integral of the power series *z*,
    # truncated to length *len*.

    void _arb_hypgeom_ci_asymp(arb_t res, const arb_t z, long N, long prec)
    void _arb_hypgeom_ci_2f3(arb_t res, const arb_t z, long N, long wp, long prec)
    void arb_hypgeom_ci(arb_t res, const arb_t z, long prec)
    # Computes the cosine integral `\operatorname{Ci}(z)`.
    # The result is indeterminate if `z < 0` since the value of the function would be complex.

    void _arb_hypgeom_ci_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_ci_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the cosine integral of the power series *z*,
    # truncated to length *len*.

    void arb_hypgeom_shi(arb_t res, const arb_t z, long prec)
    # Computes the hyperbolic sine integral `\operatorname{Shi}(z) = -i \operatorname{Si}(iz)`.

    void _arb_hypgeom_shi_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_shi_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the hyperbolic sine integral of the power series *z*,
    # truncated to length *len*.

    void arb_hypgeom_chi(arb_t res, const arb_t z, long prec)
    # Computes the hyperbolic cosine integral `\operatorname{Chi}(z)`.
    # The result is indeterminate if `z < 0` since the value of the function would be complex.

    void _arb_hypgeom_chi_series(arb_ptr res, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_chi_series(arb_poly_t res, const arb_poly_t z, long len, long prec)
    # Computes the hyperbolic cosine integral of the power series *z*,
    # truncated to length *len*.

    void arb_hypgeom_li(arb_t res, const arb_t z, int offset, long prec)
    # If *offset* is zero, computes the logarithmic integral
    # `\operatorname{li}(z) = \operatorname{Ei}(\log(z))`.
    # If *offset* is nonzero, computes the offset logarithmic integral
    # `\operatorname{Li}(z) = \operatorname{li}(z) - \operatorname{li}(2)`.
    # The result is indeterminate if `z < 0` since the value of the function would be complex.

    void _arb_hypgeom_li_series(arb_ptr res, arb_srcptr z, long zlen, int offset, long len, long prec)
    void arb_hypgeom_li_series(arb_poly_t res, const arb_poly_t z, int offset, long len, long prec)
    # Computes the logarithmic integral (optionally the offset version)
    # of the power series *z*, truncated to length *len*.

    void arb_hypgeom_bessel_j(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes the Bessel function of the first kind `J_{\nu}(z)`.

    void arb_hypgeom_bessel_y(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes the Bessel function of the second kind `Y_{\nu}(z)`.

    void arb_hypgeom_bessel_jy(arb_t res1, arb_t res2, const arb_t nu, const arb_t z, long prec)
    # Sets *res1* to `J_{\nu}(z)` and *res2* to `Y_{\nu}(z)`, computed
    # simultaneously.

    void arb_hypgeom_bessel_i(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes the modified Bessel function of the first kind
    # `I_{\nu}(z) = z^{\nu} (iz)^{-\nu} J_{\nu}(iz)`.

    void arb_hypgeom_bessel_i_scaled(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes the function `e^{-z} I_{\nu}(z)`.

    void arb_hypgeom_bessel_k(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes the modified Bessel function of the second kind `K_{\nu}(z)`.

    void arb_hypgeom_bessel_k_scaled(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes the function `e^{z} K_{\nu}(z)`.

    void arb_hypgeom_bessel_i_integration(arb_t res, const arb_t nu, const arb_t z, int scaled, long prec)
    void arb_hypgeom_bessel_k_integration(arb_t res, const arb_t nu, const arb_t z, int scaled, long prec)
    # Computes the modified Bessel functions using numerical integration.

    void arb_hypgeom_airy(arb_t ai, arb_t ai_prime, arb_t bi, arb_t bi_prime, const arb_t z, long prec)
    # Computes the Airy functions `(\operatorname{Ai}(z), \operatorname{Ai}'(z), \operatorname{Bi}(z), \operatorname{Bi}'(z))`
    # simultaneously. Any of the four function values can be omitted by passing
    # *NULL* for the unwanted output variables, speeding up the evaluation.

    void arb_hypgeom_airy_jet(arb_ptr ai, arb_ptr bi, const arb_t z, long len, long prec)
    # Writes to *ai* and *bi* the respective Taylor expansions of the Airy functions
    # at the point *z*, truncated to length *len*.
    # Either of the outputs can be *NULL* to avoid computing that function.
    # The variable *z* is not allowed to be aliased with the outputs.
    # To simplify the implementation, this method does not compute the
    # series expansions of the primed versions directly; these are
    # easily obtained by computing one extra coefficient and differentiating
    # the output with :func:`_arb_poly_derivative`.

    void _arb_hypgeom_airy_series(arb_ptr ai, arb_ptr ai_prime, arb_ptr bi, arb_ptr bi_prime, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_airy_series(arb_poly_t ai, arb_poly_t ai_prime, arb_poly_t bi, arb_poly_t bi_prime, const arb_poly_t z, long len, long prec)
    # Computes the Airy functions evaluated at the power series *z*,
    # truncated to length *len*. As with the other Airy methods, any of the
    # outputs can be *NULL*.

    void arb_hypgeom_airy_zero(arb_t a, arb_t a_prime, arb_t b, arb_t b_prime, const fmpz_t n, long prec)
    # Computes the *n*-th real zero `a_n`, `a'_n`, `b_n`, or `b'_n`
    # for the respective Airy function or Airy function derivative.
    # Any combination of the four output variables can be *NULL*.
    # The zeros are indexed by increasing magnitude, starting with
    # `n = 1` to follow the convention in the literature.
    # An index *n* that is not positive is invalid input.
    # The implementation uses asymptotic expansions for the zeros
    # [PS1991]_ together with the interval Newton method for refinement.

    void arb_hypgeom_coulomb(arb_t F, arb_t G, const arb_t l, const arb_t eta, const arb_t z, long prec)
    # Writes to *F*, *G* the values of the respective
    # Coulomb wave functions `F_{\ell}(\eta,z)` and `G_{\ell}(\eta,z)`.
    # Either of the outputs can be *NULL*.

    void arb_hypgeom_coulomb_jet(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, const arb_t z, long len, long prec)
    # Writes to *F*, *G* the respective Taylor expansions of the
    # Coulomb wave functions at the point *z*, truncated to length *len*.
    # Either of the outputs can be *NULL*.

    void _arb_hypgeom_coulomb_series(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, arb_srcptr z, long zlen, long len, long prec)
    void arb_hypgeom_coulomb_series(arb_poly_t F, arb_poly_t G, const arb_t l, const arb_t eta, const arb_poly_t z, long len, long prec)
    # Computes the Coulomb wave functions evaluated at the power series *z*,
    # truncated to length *len*. Either of the outputs can be *NULL*.

    void arb_hypgeom_chebyshev_t(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_chebyshev_u(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t z, long prec)
    void arb_hypgeom_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t z, long prec)
    void arb_hypgeom_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t z, long prec)
    void arb_hypgeom_hermite_h(arb_t res, const arb_t nu, const arb_t z, long prec)
    # Computes Chebyshev, Jacobi, Gegenbauer, Laguerre or Hermite polynomials,
    # or their extensions to non-integer orders.

    void arb_hypgeom_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, long prec)
    void arb_hypgeom_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, long prec)
    # Computes Legendre functions of the first and second kind.
    # See :func:`acb_hypgeom_legendre_p` and :func:`acb_hypgeom_legendre_q`
    # for definitions.

    void arb_hypgeom_legendre_p_ui_deriv_bound(mag_t dp, mag_t dp2, unsigned long n, const arb_t x, const arb_t x2sub1)
    # Sets *dp* to an upper bound for `P'_n(x)` and *dp2* to an upper
    # bound for `P''_n(x)` given *x* assumed to represent a real
    # number with `|x| \le 1`. The variable *x2sub1* must contain
    # the precomputed value `1-x^2` (or `x^2-1`). This method is used
    # internally to bound the propagated error for Legendre polynomials.

    void arb_hypgeom_legendre_p_ui_zero(arb_t res, arb_t res_prime, unsigned long n, const arb_t x, long K, long prec)
    void arb_hypgeom_legendre_p_ui_one(arb_t res, arb_t res_prime, unsigned long n, const arb_t x, long K, long prec)
    void arb_hypgeom_legendre_p_ui_asymp(arb_t res, arb_t res_prime, unsigned long n, const arb_t x, long K, long prec)
    void arb_hypgeom_legendre_p_ui_rec(arb_t res, arb_t res_prime, unsigned long n, const arb_t x, long prec)
    void arb_hypgeom_legendre_p_ui(arb_t res, arb_t res_prime, unsigned long n, const arb_t x, long prec)
    # Evaluates the ordinary Legendre polynomial `P_n(x)`. If *res_prime* is
    # non-NULL, simultaneously evaluates the derivative `P'_n(x)`.
    # The overall algorithm is described in [JM2018]_.
    # The versions *zero*, *one* respectively use the hypergeometric series
    # expansions at `x = 0` and `x = 1` while the *asymp* version uses an
    # asymptotic series on `(-1,1)` intended for large *n*. The parameter *K*
    # specifies the exact number of expansion terms to use (if the series
    # expansion truncated at this point does not give the exact polynomial,
    # an error bound is computed automatically).
    # The asymptotic expansion with error bounds is given in [Bog2012]_.
    # The *rec* version uses the forward recurrence implemented using
    # fixed-point arithmetic; it is only intended for the interval `(-1,1)`,
    # moderate *n* and modest precision.
    # The default version attempts to choose the best algorithm automatically.
    # It also estimates the amount of cancellation in the hypergeometric series
    # and increases the working precision to compensate, bounding the
    # propagated error using derivative bounds.

    void arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, unsigned long n, unsigned long k, long prec)
    # Sets *res* to the *k*-th root of the Legendre polynomial `P_n(x)`.
    # We index the roots in decreasing order
    # .. math ::
    # 1 > x_0 > x_1 > \ldots > x_{n-1} > -1
    # (which corresponds to ordering the roots of `P_n(\cos(\theta))`
    # in order of increasing `\theta`).
    # If *weight* is non-NULL, it is set to the weight corresponding
    # to the node `x_k` for Gaussian quadrature on `[-1,1]`.
    # Note that only `\lceil n / 2 \rceil` roots need to be computed,
    # since the remaining roots are given by `x_k = -x_{n-1-k}`.
    # We compute an enclosing interval using an asymptotic approximation followed
    # by some number of Newton iterations, using the error bounds given
    # in [Pet1999]_. If very high precision is requested, the root is
    # subsequently refined using interval Newton steps with doubling working
    # precision.

    void arb_hypgeom_dilog(arb_t res, const arb_t z, long prec)
    # Computes the dilogarithm `\operatorname{Li}_2(z)`.

    void arb_hypgeom_sum_fmpq_arb_forward(arb_t res, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    void arb_hypgeom_sum_fmpq_arb_rs(arb_t res, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    void arb_hypgeom_sum_fmpq_arb(arb_t res, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    # Sets *res* to the finite hypergeometric sum
    # `\sum_{n=0}^{N-1} (\textbf{a})_n z^n / (\textbf{b})_n`
    # where `\textbf{x}_n = (x_1)_n (x_2)_n \cdots`,
    # given vectors of rational parameters *a* (of length *alen*)
    # and *b* (of length *blen*).
    # If *reciprocal* is set, replace `z` by `1 / z`.
    # The *forward* version uses the forward recurrence, optimized by
    # delaying divisions, the *rs* version
    # uses rectangular splitting, and the default version uses
    # an automatic algorithm choice.

    void arb_hypgeom_sum_fmpq_imag_arb_forward(arb_t res1, arb_t res2, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    void arb_hypgeom_sum_fmpq_imag_arb_rs(arb_t res1, arb_t res2, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    void arb_hypgeom_sum_fmpq_imag_arb_bs(arb_t res1, arb_t res2, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    void arb_hypgeom_sum_fmpq_imag_arb(arb_t res1, arb_t res2, const fmpq * a, long alen, const fmpq * b, long blen, const arb_t z, int reciprocal, long N, long prec)
    # Sets *res1* and *res2* to the real and imaginary part of the
    # finite hypergeometric sum
    # `\sum_{n=0}^{N-1} (\textbf{a})_n (i z)^n / (\textbf{b})_n`.
    # If *reciprocal* is set, replace `z` by `1 / z`.
