# distutils: libraries = flint
# distutils: depends = flint/acb_dirichlet.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void acb_dirichlet_roots_init(acb_dirichlet_roots_t roots, ulong n, slong num, slong prec)
    # Initializes *roots* with precomputed data for fast evaluation of roots of
    # unity `e^{2\pi i k/n}` of a fixed order *n*. The precomputation is
    # optimized for *num* evaluations.
    # For very small *num*, only the single root `e^{2\pi i/n}` will be
    # precomputed, which can then be raised to a power. For small *prec*
    # and large *n*, this method might even skip precomputing this single root
    # if it estimates that evaluating roots of unity from scratch will be faster
    # than powering.
    # If *num* is large enough, the whole set of roots in the first quadrant
    # will be precomputed at once. However, this is automatically avoided for
    # large *n* if too much memory would be used. For intermediate *num*,
    # baby-step giant-step tables are computed.

    void acb_dirichlet_roots_clear(acb_dirichlet_roots_t roots)
    # Clears the structure.

    void acb_dirichlet_root(acb_t res, const acb_dirichlet_roots_t roots, ulong k, slong prec)
    # Computes `e^{2\pi i k/n}`.

    void acb_dirichlet_powsum_term(acb_ptr res, arb_t log_prev, ulong * prev, const acb_t s, ulong k, int integer, int critical_line, slong len, slong prec)
    # Sets *res* to `k^{-(s+x)}` as a power series in *x* truncated to length *len*.
    # The flags *integer* and *critical_line* respectively specify optimizing
    # for *s* being an integer or having real part 1/2.
    # On input *log_prev* should contain the natural logarithm of the integer
    # at *prev*. If *prev* is close to *k*, this can be used to speed up
    # computations. If `\log(k)` is computed internally by this function, then
    # *log_prev* is overwritten by this value, and the integer at *prev* is
    # overwritten by *k*, allowing *log_prev* to be recycled for the next
    # term when evaluating a power sum.

    void acb_dirichlet_powsum_sieved(acb_ptr res, const acb_t s, ulong n, slong len, slong prec)
    # Sets *res* to `\sum_{k=1}^n k^{-(s+x)}`
    # as a power series in *x* truncated to length *len*.
    # This function stores a table of powers that have already been calculated,
    # computing `(ij)^r` as `i^r j^r` whenever `k = ij` is
    # composite. As a further optimization, it groups all even `k` and
    # evaluates the sum as a polynomial in `2^{-(s+x)}`.
    # This scheme requires about `n / \log n` powers, `n / 2` multiplications,
    # and temporary storage of `n / 6` power series. Due to the extra
    # power series multiplications, it is only faster than the naive
    # algorithm when *len* is small.

    void acb_dirichlet_powsum_smooth(acb_ptr res, const acb_t s, ulong n, slong len, slong prec)
    # Sets *res* to `\sum_{k=1}^n k^{-(s+x)}`
    # as a power series in *x* truncated to length *len*.
    # This function performs partial sieving by adding multiples of 5-smooth *k*
    # into separate buckets. Asymptotically, this requires computing 4/15
    # of the powers, which is slower than *sieved*, but only requires
    # logarithmic extra space. It is also faster for large *len*, since most
    # power series multiplications are traded for additions.
    # A slightly bigger gain for larger *n* could be achieved by using more
    # small prime factors, at the expense of space.

    void acb_dirichlet_zeta(acb_t res, const acb_t s, slong prec)
    # Computes `\zeta(s)` using an automatic choice of algorithm.

    void acb_dirichlet_zeta_jet(acb_t res, const acb_t s, int deflate, slong len, slong prec)
    # Computes the first *len* terms of the Taylor series of the Riemann zeta
    # function at *s*. If *deflate* is nonzero, computes the deflated
    # function `\zeta(s) - 1/(s-1)` instead.

    void acb_dirichlet_zeta_bound(mag_t res, const acb_t s)
    # Computes an upper bound for `|\zeta(s)|` quickly. On the critical strip (and
    # slightly outside of it), formula (43.3) in [Rad1973]_ is used.
    # To the right, evaluating at the real part of *s* gives a trivial bound.
    # To the left, the functional equation is used.

    void acb_dirichlet_zeta_deriv_bound(mag_t der1, mag_t der2, const acb_t s)
    # Sets *der1* to a bound for `|\zeta'(s)|` and *der2* to a bound for
    # `|\zeta''(s)|`. These bounds are mainly intended for use in the critical
    # strip and will not be tight.

    void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec)
    # Sets *res* to the Dirichlet eta function
    # `\eta(s) = \sum_{k=1}^{\infty} (-1)^{k+1} / k^s = (1-2^{1-s}) \zeta(s)`,
    # also known as the alternating zeta function.
    # Note that the alternating character `\{1,-1\}` is not itself
    # a Dirichlet character.

    void acb_dirichlet_xi(acb_t res, const acb_t s, slong prec)
    # Sets *res* to the Riemann xi function
    # `\xi(s) = \frac{1}{2} s (s-1) \pi^{-s/2} \Gamma(\frac{1}{2} s) \zeta(s)`.
    # The functional equation for xi is `\xi(1-s) = \xi(s)`.

    void acb_dirichlet_zeta_rs_f_coeffs(acb_ptr f, const arb_t p, slong n, slong prec)
    # Computes the coefficients `F^{(j)}(p)` for `0 \le j < n`.
    # Uses power series division. This method breaks down when `p = \pm 1/2`
    # (which is not problem if *s* is an exact floating-point number).

    void acb_dirichlet_zeta_rs_d_coeffs(arb_ptr d, const arb_t sigma, slong k, slong prec)
    # Computes the coefficients `d_j^{(k)}` for `0 \le j \le \lfloor 3k/2 \rfloor + 1`.
    # On input, the array *d* must contain the coefficients for `d_j^{(k-1)}`
    # unless `k = 0`, and these coefficients will be updated in-place.

    void acb_dirichlet_zeta_rs_bound(mag_t err, const acb_t s, slong K)
    # Bounds the error term `RS_K` following Theorem 4.2 in Arias de Reyna.

    void acb_dirichlet_zeta_rs_r(acb_t res, const acb_t s, slong K, slong prec)
    # Computes `\mathcal{R}(s)` in the upper half plane. Uses precisely *K*
    # asymptotic terms in the RS formula if this input parameter is positive;
    # otherwise chooses the number of terms automatically based on *s* and the
    # precision.

    void acb_dirichlet_zeta_rs(acb_t res, const acb_t s, slong K, slong prec)
    # Computes `\zeta(s)` using the Riemann-Siegel formula. Uses precisely
    # *K* asymptotic terms in the RS formula if this input parameter is positive;
    # otherwise chooses the number of terms automatically based on *s* and the
    # precision.

    void acb_dirichlet_zeta_jet_rs(acb_ptr res, const acb_t s, slong len, slong prec)
    # Computes the first *len* terms of the Taylor series of the Riemann zeta
    # function at *s* using the Riemann Siegel formula. This function currently
    # only supports *len* = 1 or *len* = 2. A finite difference is used
    # to compute the first derivative.

    void acb_dirichlet_hurwitz(acb_t res, const acb_t s, const acb_t a, slong prec)
    # Computes the Hurwitz zeta function `\zeta(s, a)`.
    # This function automatically delegates to the code for the Riemann zeta function
    # when `a = 1`. Some other special cases may also be handled by direct
    # formulas. In general, Euler-Maclaurin summation is used.

    void acb_dirichlet_hurwitz_precomp_init(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, slong A, slong K, slong N, slong prec)
    # Precomputes a grid of Taylor polynomials for fast evaluation of
    # `\zeta(s,a)` on `a \in (0,1]` with fixed *s*.
    # *A* is the initial shift to apply to *a*, *K* is the number of Taylor terms,
    # *N* is the number of grid points.  The precomputation requires *NK*
    # evaluations of the Hurwitz zeta function, and each subsequent evaluation
    # requires *2K* simple arithmetic operations (polynomial evaluation) plus
    # *A* powers. As *K* grows, the error is at most `O(1/(2AN)^K)`.
    # This function can be called with *A* set to zero, in which case
    # no Taylor series precomputation is performed. This means that evaluation
    # will be identical to calling :func:`acb_dirichlet_hurwitz` directly.
    # Otherwise, we require that *A*, *K* and *N* are all positive. For a finite
    # error bound, we require `K+\operatorname{re}(s) > 1`.
    # To avoid an initial "bump" that steals precision
    # and slows convergence, *AN* should be at least roughly as large as `|s|`,
    # e.g. it is a good idea to have at least `AN > 0.5 |s|`.
    # If *deflate* is set, the deflated Hurwitz zeta function is used,
    # removing the pole at `s = 1`.

    void acb_dirichlet_hurwitz_precomp_init_num(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, double num_eval, slong prec)
    # Initializes *pre*, choosing the parameters *A*, *K*, and *N*
    # automatically to minimize the cost of *num_eval* evaluations of the
    # Hurwitz zeta function at argument *s* to precision *prec*.

    void acb_dirichlet_hurwitz_precomp_clear(acb_dirichlet_hurwitz_precomp_t pre)
    # Clears the precomputed data.

    void acb_dirichlet_hurwitz_precomp_choose_param(ulong * A, ulong * K, ulong * N, const acb_t s, double num_eval, slong prec)
    # Chooses precomputation parameters *A*, *K* and *N* to minimize
    # the cost of *num_eval* evaluations of the Hurwitz zeta function
    # at argument *s* to precision *prec*.
    # If it is estimated that evaluating each Hurwitz zeta function from
    # scratch would be better than performing a precomputation, *A*, *K* and *N*
    # are all set to 0.

    void acb_dirichlet_hurwitz_precomp_bound(mag_t res, const acb_t s, slong A, slong K, slong N)
    # Computes an upper bound for the truncation error (not accounting for
    # roundoff error) when evaluating `\zeta(s,a)` with precomputation parameters
    # *A*, *K*, *N*, assuming that `0 < a \le 1`.
    # For details, see :ref:`algorithms_hurwitz`.

    void acb_dirichlet_hurwitz_precomp_eval(acb_t res, const acb_dirichlet_hurwitz_precomp_t pre, ulong p, ulong q, slong prec)
    # Evaluates `\zeta(s,p/q)` using precomputed data, assuming that `0 < p/q \le 1`.

    void acb_dirichlet_lerch_phi_integral(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
    void acb_dirichlet_lerch_phi_direct(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
    void acb_dirichlet_lerch_phi(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
    # Computes the Lerch transcendent
    # .. math ::
    # \Phi(z,s,a) = \sum_{k=0}^{\infty} \frac{z^k}{(k+a)^s}
    # which is analytically continued for `|z| \ge 1`.
    # The *direct* version evaluates a truncation of the defining series.
    # The *integral* version uses the Hankel contour integral
    # .. math ::
    # \Phi(z,s,a) = -\frac{\Gamma(1-s)}{2 \pi i} \int_C \frac{(-t)^{s-1} e^{-a t}}{1 - z e^{-t}} dt
    # where the path is deformed as needed to avoid poles and branch
    # cuts of the integrand.
    # The default method chooses an algorithm automatically and also
    # checks for some special cases where the function can be expressed
    # in terms of simpler functions (Hurwitz zeta, polylogarithms).

    void acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, slong prec)
    # Given a nonnegative integer *n*, sets *res* to the generalized Stieltjes constant
    # `\gamma_n(a)` which is the coefficient in the Laurent series of the
    # Hurwitz zeta function at the pole
    # .. math ::
    # \zeta(s,a) = \frac{1}{s-1} + \sum_{n=0}^\infty \frac{(-1)^n}{n!} \gamma_n(a) (s-1)^n.
    # With `a = 1`, this gives the ordinary Stieltjes constants for the
    # Riemann zeta function.
    # This function uses an integral representation to permit fast computation
    # for extremely large *n* [JB2018]_. If *n* is moderate and the precision
    # is high enough, it falls back to evaluating the Hurwitz zeta function
    # of a power series and reading off the last coefficient.
    # Note that for computing a range of values
    # `\gamma_0(a), \ldots, \gamma_n(a)`, it is
    # generally more efficient to evaluate the Hurwitz zeta function series
    # expansion once at `s = 1` than to call this function repeatedly,
    # unless *n* is extremely large (at least several hundred).

    void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, slong prec)
    # Sets *res* to `\chi(n)`, the value of the Dirichlet character *chi*
    # at the integer *n*.

    void acb_dirichlet_chi_vec(acb_ptr v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv, slong prec)
    # Compute the *nv* first Dirichlet values.

    void acb_dirichlet_pairing(acb_t res, const dirichlet_group_t G, ulong m, ulong n, slong prec)

    void acb_dirichlet_pairing_char(acb_t res, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b, slong prec)
    # Sets *res* to the value of the Dirichlet pairing `\chi(m,n)` at numbers `m` and `n`.
    # The second form takes two characters as input.

    void acb_dirichlet_gauss_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_gauss_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_gauss_sum_order2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_gauss_sum_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_gauss_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_jacobi_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)

    void acb_dirichlet_jacobi_sum_factor(acb_t res,  const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)

    void acb_dirichlet_jacobi_sum_gauss(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)

    void acb_dirichlet_jacobi_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1,  const dirichlet_char_t chi2, slong prec)

    void acb_dirichlet_jacobi_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, ulong b, slong prec)

    void acb_dirichlet_chi_theta_arb(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t t, slong prec)

    void acb_dirichlet_ui_theta_arb(acb_t res, const dirichlet_group_t G, ulong a, const arb_t t, slong prec)
    # Compute the theta series `\Theta_q(a,t)` for real argument `t>0`.
    # Beware that if `t<1` the functional equation
    # .. math::
    # t \theta(a,t) = \epsilon(\chi) \theta\left(\frac1a, \frac1t\right)
    # should be used, which is not done automatically (to avoid recomputing the
    # Gauss sum).
    # We call *theta series* of a Dirichlet character the quadratic series
    # .. math::
    # \Theta_q(a) = \sum_{n\geq 0} \chi_q(a, n) n^p x^{n^2}
    # where `p` is the parity of the character `\chi_q(a,\cdot)`.
    # For `\Re(t)>0` we write `x(t)=\exp(-\frac{\pi}{N}t^2)` and define
    # .. math::
    # \Theta_q(a,t) = \sum_{n\geq 0} \chi_q(a, n) x(t)^{n^2}.

    ulong acb_dirichlet_theta_length(ulong q, const arb_t t, slong prec)

    void acb_dirichlet_qseries_arb_powers_naive(acb_t res, const arb_t x, int p, const ulong * a, const acb_dirichlet_roots_t z, slong len, slong prec)

    void acb_dirichlet_qseries_arb_powers_smallorder(acb_t res, const arb_t x, int p, const ulong * a, const acb_dirichlet_roots_t z, slong len, slong prec)

    void acb_dirichlet_dft_conrey(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

    void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

    void acb_dirichlet_root_number_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_root_number(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
    # Computes `L(s,\chi)` using decomposition in terms of the Hurwitz zeta function
    # .. math::
    # L(s,\chi) = q^{-s}\sum_{k=1}^q \chi(k) \,\zeta\!\left(s,\frac kq\right).
    # If `s = 1` and `\chi` is non-principal, the deflated Hurwitz zeta function
    # is used to avoid poles.
    # If *precomp* is *NULL*, each Hurwitz zeta function value is computed
    # directly. If a pre-initialized *precomp* object is provided, this will be
    # used instead to evaluate the Hurwitz zeta function.

    void acb_dirichlet_l_euler_product(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

    void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s, const signed char * chi, int mod, int reciprocal, slong prec)
    # Computes `L(s,\chi)` directly using the Euler product. This is
    # efficient if *s* has large positive real part. As implemented, this
    # function only gives a finite result if `\operatorname{re}(s) \ge 2`.
    # An error bound is computed via :func:`mag_hurwitz_zeta_uiui`.
    # If *s* is complex, replace it with its real part. Since
    # .. math ::
    # \frac{1}{L(s,\chi)} = \prod_{p} \left(1 - \frac{\chi(p)}{p^s}\right)
    # = \sum_{k=1}^{\infty} \frac{\mu(k)\chi(k)}{k^s}
    # and the truncated product gives all smooth-index terms in the series, we have
    # .. math ::
    # \left|\prod_{p < N} \left(1 - \frac{\chi(p)}{p^s}\right) - \frac{1}{L(s,\chi)}\right|
    # \le \sum_{k=N}^{\infty} \frac{1}{k^s} = \zeta(s,N).
    # The underscore version specialized for integer *s* assumes that `\chi` is
    # a real Dirichlet character given by the explicit list *chi* of character
    # values at 0, 1, ..., *mod* - 1. If *reciprocal* is set, it computes
    # `1 / L(s,\chi)` (this is faster if the reciprocal can be used directly).

    void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
    # Computes `L(s,\chi)` using a default choice of algorithm.

    void acb_dirichlet_l_fmpq(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
    void acb_dirichlet_l_fmpq_afe(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
    # Computes `L(s,\chi)` where *s* is a rational number.
    # The *afe* version uses the approximate functional equation;
    # the default version chooses an algorithm automatically.

    void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp, const dirichlet_group_t G, slong prec)
    # Compute all values `L(s,\chi)` for `\chi` mod `q`, using the
    # Hurwitz zeta function and a discrete Fourier transform.
    # The output *res* is assumed to have length *G->phi_q* and values
    # are stored by lexicographically ordered
    # Conrey logs. See :func:`acb_dirichlet_dft_conrey`.
    # If *precomp* is *NULL*, each Hurwitz zeta function value is computed
    # directly. If a pre-initialized *precomp* object is provided, this will be
    # used instead to evaluate the Hurwitz zeta function.

    void acb_dirichlet_l_jet(acb_ptr res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)
    # Computes the Taylor expansion of `L(s,\chi)` to length *len*,
    # i.e. `L(s), L'(s), \ldots, L^{(len-1)}(s) / (len-1)!`.
    # If *deflate* is set, computes the expansion of
    # .. math ::
    # L(s,\chi) - \frac{\sum_{k=1}^q \chi(k)}{(s-1)q}
    # instead. If *chi* is a principal character, then this has the effect of
    # subtracting the pole with residue `\sum_{k=1}^q \chi(k) = \phi(q) / q`
    # that is located at `s = 1`. In particular, when evaluated at `s = 1`, this
    # gives the regular part of the Laurent expansion.
    # When *chi* is non-principal, *deflate* has no effect.

    void _acb_dirichlet_l_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)

    void acb_dirichlet_l_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)
    # Sets *res* to the power series `L(s,\chi)` where *s* is a given power series, truncating the result to length *len*.
    # See :func:`acb_dirichlet_l_jet` for the meaning of the *deflate* flag.

    void acb_dirichlet_hardy_theta(acb_ptr res, const acb_t t, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
    # Computes the phase function used to construct the Z-function.
    # We have
    # .. math ::
    # \theta(t) = -\frac{t}{2} \log(\pi/q) - \frac{i \log(\epsilon)}{2}
    # + \frac{\log \Gamma((s+\delta)/2) - \log \Gamma((1-s+\delta)/2)}{2i}
    # where `s = 1/2+it`, `\delta` is the parity of *chi*, and `\epsilon`
    # is the root number as computed by :func:`acb_dirichlet_root_number`.
    # The first *len* terms in the Taylor expansion are written to the output.

    void acb_dirichlet_hardy_z(acb_ptr res, const acb_t t, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
    # Computes the Hardy Z-function, also known as the Riemann-Siegel Z-function
    # `Z(t) = e^{i \theta(t)} L(1/2+it)`, which is real-valued for real *t*.
    # The first *len* terms in the Taylor expansion are written to the output.

    void _acb_dirichlet_hardy_theta_series(acb_ptr res, acb_srcptr t, slong tlen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)

    void acb_dirichlet_hardy_theta_series(acb_poly_t res, const acb_poly_t t, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
    # Sets *res* to the power series `\theta(t)` where *t* is a given power series, truncating the result to length *len*.

    void _acb_dirichlet_hardy_z_series(acb_ptr res, acb_srcptr t, slong tlen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)

    void acb_dirichlet_hardy_z_series(acb_poly_t res, const acb_poly_t t, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
    # Sets *res* to the power series `Z(t)` where *t* is a given power series, truncating the result to length *len*.

    void acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
    # Sets *res* to the *n*-th Gram point `g_n`, defined as the unique solution
    # in `[7, \infty)` of `\theta(g_n) = \pi n`. Currently only the Gram points
    # corresponding to the Riemann zeta function are supported and *G* and *chi*
    # must both be set to *NULL*. Requires `n \ge -1`.

    ulong acb_dirichlet_turing_method_bound(const fmpz_t p)
    # Computes an upper bound *B* for the minimum number of consecutive good
    # Gram blocks sufficient to count nontrivial zeros of the Riemann zeta
    # function using Turing's method [Tur1953]_ as updated by [Leh1970]_,
    # [Bre1979]_, and [Tru2011]_.
    # Let `N(T)` denote the number of zeros (counted according to their
    # multiplicities) of `\zeta(s)` in the region `0 < \operatorname{Im}(s) \le T`.
    # If at least *B* consecutive Gram blocks with union `[g_n, g_p)`
    # satisfy Rosser's rule, then `N(g_n) \le n + 1` and `N(g_p) \ge p + 1`.

    int _acb_dirichlet_definite_hardy_z(arb_t res, const arf_t t, slong * pprec)
    # Sets *res* to the Hardy Z-function `Z(t)`.
    # The initial precision (* *pprec*) is increased as necessary
    # to determine the sign of `Z(t)`. The sign is returned.

    void _acb_dirichlet_isolate_gram_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
    # Uses Gram's law to compute an interval `(a, b)` that
    # contains the *n*-th zero of the Hardy Z-function and no other zero.
    # Requires `1 \le n \le 126`.

    void _acb_dirichlet_isolate_rosser_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
    # Uses Rosser's rule to compute an interval `(a, b)` that
    # contains the *n*-th zero of the Hardy Z-function and no other zero.
    # Requires `1 \le n \le 13999526`.

    void _acb_dirichlet_isolate_turing_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
    # Computes an interval `(a, b)` that contains the *n*-th zero of the
    # Hardy Z-function and no other zero, following Turing's method.
    # Requires `n \ge 2`.

    void acb_dirichlet_isolate_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
    # Computes an interval `(a, b)` that contains the *n*-th zero of the
    # Hardy Z-function and contains no other zero, using the most appropriate
    # underscore version of this function. Requires `n \ge 1`.

    void _acb_dirichlet_refine_hardy_z_zero(arb_t res, const arf_t a, const arf_t b, slong prec)
    # Sets *res* to the unique zero of the Hardy Z-function in the
    # interval `(a, b)`.

    void acb_dirichlet_hardy_z_zero(arb_t res, const fmpz_t n, slong prec)
    # Sets *res* to the *n*-th zero of the Hardy Z-function, requiring `n \ge 1`.

    void acb_dirichlet_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)
    # Sets the entries of *res* to *len* consecutive zeros of the
    # Hardy Z-function, beginning with the *n*-th zero. Requires positive *n*.

    void acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, slong prec)
    # Sets *res* to the *n*-th nontrivial zero of `\zeta(s)`, requiring `n \ge 1`.

    void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec)
    # Sets the entries of *res* to *len* consecutive nontrivial zeros of `\zeta(s)`
    # beginning with the *n*-th zero. Requires positive *n*.

    void _acb_dirichlet_exact_zeta_nzeros(fmpz_t res, const arf_t t)

    void acb_dirichlet_zeta_nzeros(arb_t res, const arb_t t, slong prec)
    # Compute the number of zeros (counted according to their multiplicities)
    # of `\zeta(s)` in the region `0 < \operatorname{Im}(s) \le t`.

    void acb_dirichlet_backlund_s(arb_t res, const arb_t t, slong prec)
    # Compute `S(t) = \frac{1}{\pi}\operatorname{arg}\zeta(\frac{1}{2} + it)`
    # where the argument is defined by continuous variation of `s` in `\zeta(s)`
    # starting at `s = 2`, then vertically to `s = 2 + it`, then horizontally
    # to `s = \frac{1}{2} + it`. In particular `\operatorname{arg}` in this
    # context is not the principal value of the argument, and it cannot be
    # computed directly by :func:`acb_arg`. In practice `S(t)` is computed as
    # `S(t) = N(t) - \frac{1}{\pi}\theta(t) - 1` where `N(t)` is
    # :func:`acb_dirichlet_zeta_nzeros` and `\theta(t)` is
    # :func:`acb_dirichlet_hardy_theta`.

    void acb_dirichlet_backlund_s_bound(mag_t res, const arb_t t)
    # Compute an upper bound for `|S(t)|` quickly. Theorem 1
    # and the bounds in (1.2) in [Tru2014]_ are used.

    void acb_dirichlet_zeta_nzeros_gram(fmpz_t res, const fmpz_t n)
    # Compute `N(g_n)`. That is, compute the number of zeros (counted according
    # to their multiplicities) of `\zeta(s)` in the region
    # `0 < \operatorname{Im}(s) \le g_n` where `g_n` is the *n*-th Gram point.
    # Requires `n \ge -1`.

    slong acb_dirichlet_backlund_s_gram(const fmpz_t n)
    # Compute `S(g_n)` where `g_n` is the *n*-th Gram point. Requires `n \ge -1`.

    void acb_dirichlet_platt_scaled_lambda(arb_t res, const arb_t t, slong prec)
    # Compute `\Lambda(t) e^{\pi t/4}` where
    # .. math ::
    # \Lambda(t) = \pi^{-\frac{it}{2}}
    # \Gamma\left(\frac{\frac{1}{2}+it}{2}\right)
    # \zeta\left(\frac{1}{2} + it\right)
    # is defined in the beginning of section 3 of [Pla2017]_. As explained in
    # [Pla2011]_ this function has the same zeros as `\zeta(1/2 + it)` and is
    # real-valued by the functional equation, and the exponential factor is
    # designed to counteract the decay of the gamma factor as `t` increases.

    void acb_dirichlet_platt_scaled_lambda_vec(arb_ptr res, const fmpz_t T, slong A, slong B, slong prec)

    void acb_dirichlet_platt_multieval(arb_ptr res, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma, slong prec)

    void acb_dirichlet_platt_multieval_threaded(arb_ptr res, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma, slong prec)
    # Compute :func:`acb_dirichlet_platt_scaled_lambda` at `N=AB` points on a
    # grid, following the notation of [Pla2017]_. The first point on the grid
    # is `T - B/2` and the distance between grid points is `1/A`. The product
    # `N=AB` must be an even integer. The multieval versions evaluate the
    # function at all points on the grid simultaneously using discrete Fourier
    # transforms, and they require the four additional tuning parameters
    # *h*, *J*, *K*, and *sigma*. The *threaded* multieval version splits the
    # computation over the number of threads returned by
    # *flint_get_num_threads()*, while the default multieval version chooses
    # whether to use multithreading automatically.

    void acb_dirichlet_platt_ws_interpolation(arb_t res, arf_t deriv, const arb_t t0, arb_srcptr p, const fmpz_t T, slong A, slong B, slong Ns_max, const arb_t H, slong sigma, slong prec)
    # Compute :func:`acb_dirichlet_platt_scaled_lambda` at *t0* by
    # Gaussian-windowed Whittaker-Shannon interpolation of points evaluated by
    # :func:`acb_dirichlet_platt_scaled_lambda_vec`. The derivative is
    # also approximated if the output parameter *deriv* is not *NULL*.
    # *Ns_max* defines the maximum number of supporting points to be used in
    # the interpolation on either side of *t0*. *H* is the standard deviation
    # of the Gaussian window centered on *t0* to be applied before the
    # interpolation. *sigma* is an odd positive integer tuning parameter
    # `\sigma \in 2\mathbb{Z}_{>0}+1` used in computing error bounds.

    slong _acb_dirichlet_platt_local_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma_grid, slong Ns_max, const arb_t H, slong sigma_interp, slong prec)

    slong acb_dirichlet_platt_local_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)

    slong acb_dirichlet_platt_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)
    # Sets at most the first *len* entries of *res* to consecutive
    # zeros of the Hardy Z-function starting with the *n*-th zero.
    # The number of obtained consecutive zeros is returned. The first two
    # function variants each make a single call to Platt's grid evaluation
    # of the scaled Lambda function, whereas the third variant performs as many
    # evaluations as necessary to obtain *len* consecutive zeros.
    # The final several parameters of the underscored local variant have the same
    # meanings as in the functions :func:`acb_dirichlet_platt_multieval`
    # and :func:`acb_dirichlet_platt_ws_interpolation`. The non-underscored
    # variants currently expect `10^4 \leq n \leq 10^{23}`. The user has the
    # option of multi-threading through *flint_set_num_threads(numthreads)*.

    slong acb_dirichlet_platt_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec)
    # Sets at most the first *len* entries of *res* to consecutive
    # zeros of the Riemann zeta function starting with the *n*-th zero.
    # The number of obtained consecutive zeros is returned. It currently
    # expects `10^4 \leq n \leq 10^{23}`. The user has the option of
    # multi-threading through *flint_set_num_threads(numthreads)*.
