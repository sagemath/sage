# distutils: libraries = flint
# distutils: depends = flint/hypgeom.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void hypgeom_init(hypgeom_t hyp)

    void hypgeom_clear(hypgeom_t hyp)

    slong hypgeom_estimate_terms(const mag_t z, int r, slong d)
    # Computes an approximation of the largest `n` such
    # that `|z|^n/(n!)^r = 2^{-d}`, giving a first-order estimate of the
    # number of terms needed to approximate the sum of a hypergeometric
    # series of weight `r \ge 0` and argument `z` to an absolute
    # precision of `d \ge 0` bits. If `r = 0`, the direct solution of the
    # equation is given by `n = (\log(1-z) - d \log 2) / \log z`.
    # If `r > 0`, using `\log n! \approx n \log n - n` gives an equation
    # that can be solved in terms of the Lambert *W*-function as
    # `n = (d \log 2) / (r\,W\!(t))` where
    # `t = (d \log 2) / (e r z^{1/r})`.
    # The evaluation is done using double precision arithmetic.
    # The function aborts if the computed value of `n` is greater
    # than or equal to LONG_MAX / 2.

    slong hypgeom_bound(mag_t error, int r, slong C, slong D, slong K, const mag_t TK, const mag_t z, slong prec)
    # Computes a truncation parameter sufficient to achieve *prec* bits
    # of absolute accuracy, according to the strategy described above.
    # The input consists of `r`, `C`, `D`, `K`, precomputed bound for `T(K)`,
    # and `\tilde z = z (a_p / b_q)`, such that for `k > K`, the hypergeometric
    # term ratio is bounded by
    # .. math ::
    # \frac{\tilde z}{k^r} \frac{k(k-D)}{(k-C)(k-2D)}.
    # Given this information, we compute a `\varepsilon` and an
    # integer `n` such that
    # `\left| \sum_{k=n}^{\infty} T(k) \right| \le \varepsilon \le 2^{-\mathrm{prec}}`.
    # The output variable *error* is set to the value of `\varepsilon`,
    # and `n` is returned.

    void hypgeom_precompute(hypgeom_t hyp)
    # Precomputes the bounds data `C`, `D`, `K` and an upper bound for `T(K)`.

    void arb_hypgeom_sum(arb_t P, arb_t Q, const hypgeom_t hyp, slong n, slong prec)
    # Computes `P, Q` such that `P / Q = \sum_{k=0}^{n-1} T(k)` where `T(k)`
    # is defined by *hyp*,
    # using binary splitting and a working precision of *prec* bits.

    void arb_hypgeom_infsum(arb_t P, arb_t Q, hypgeom_t hyp, slong tol, slong prec)
    # Computes `P, Q` such that `P / Q = \sum_{k=0}^{\infty} T(k)` where `T(k)`
    # is defined by *hyp*, using binary splitting and
    # working precision of *prec* bits.
    # The number of terms is chosen automatically to bound the
    # truncation error by at most `2^{-\mathrm{tol}}`.
    # The bound for the truncation error is included in the output
    # as part of *P*.
