# distutils: libraries = flint
# distutils: depends = flint/acb_modular.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void psl2z_init(psl2z_t g)
    # Initializes *g* and set it to the identity element.

    void psl2z_clear(psl2z_t g)
    # Clears *g*.

    void psl2z_swap(psl2z_t f, psl2z_t g)
    # Swaps *f* and *g* efficiently.

    void psl2z_set(psl2z_t f, const psl2z_t g)
    # Sets *f* to a copy of *g*.

    void psl2z_one(psl2z_t g)
    # Sets *g* to the identity element.

    int psl2z_is_one(const psl2z_t g)
    # Returns nonzero iff *g* is the identity element.

    void psl2z_print(const psl2z_t g)
    # Prints *g* to standard output.

    void psl2z_fprint(FILE * file, const psl2z_t g)
    # Prints *g* to the stream *file*.

    int psl2z_equal(const psl2z_t f, const psl2z_t g)
    # Returns nonzero iff *f* and *g* are equal.

    void psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g)
    # Sets *h* to the product of *f* and *g*, namely the matrix product
    # with the signs canonicalized.

    void psl2z_inv(psl2z_t h, const psl2z_t g)
    # Sets *h* to the inverse of *g*.

    int psl2z_is_correct(const psl2z_t g)
    # Returns nonzero iff *g* contains correct data, i.e.
    # satisfying `ad-bc = 1`, `c \ge 0`, and `d > 0` if `c = 0`.

    void psl2z_randtest(psl2z_t g, flint_rand_t state, slong bits)
    # Sets *g* to a random element of `\text{PSL}(2, \mathbb{Z})`
    # with entries of bit length at most *bits*
    # (or 1, if *bits* is not positive). We first generate *a* and *d*, compute
    # their Bezout coefficients, divide by the GCD, and then correct the signs.

    void acb_modular_transform(acb_t w, const psl2z_t g, const acb_t z, slong prec)
    # Applies the modular transformation *g* to the complex number *z*,
    # evaluating
    # .. math ::
    # w = g z = \frac{az+b}{cz+d}.

    void acb_modular_fundamental_domain_approx_d(psl2z_t g, double x, double y, double one_minus_eps)

    void acb_modular_fundamental_domain_approx_arf(psl2z_t g, const arf_t x, const arf_t y, const arf_t one_minus_eps, slong prec)
    # Attempts to determine a modular transformation *g* that maps the
    # complex number `x+yi` to the fundamental domain or just
    # slightly outside the fundamental domain, where the target tolerance
    # (not a strict bound) is specified by *one_minus_eps*.
    # The inputs are assumed to be finite numbers, with *y* positive.
    # Uses floating-point iteration, repeatedly applying either
    # the transformation `z \gets z + b` or `z \gets -1/z`. The iteration is
    # terminated if `|x| \le 1/2` and `x^2 + y^2 \ge 1 - \varepsilon` where
    # `1 - \varepsilon` is passed as *one_minus_eps*. It is also terminated
    # if too many steps have been taken without convergence, or if the numbers
    # end up too large or too small for the working precision.
    # The algorithm can fail to produce a satisfactory transformation.
    # The output *g* is always set to *some* correct modular transformation,
    # but it is up to the user to verify a posteriori that *g* maps `x+yi`
    # close enough to the fundamental domain.

    void acb_modular_fundamental_domain_approx(acb_t w, psl2z_t g, const acb_t z, const arf_t one_minus_eps, slong prec)
    # Attempts to determine a modular transformation *g* that maps the
    # complex number `z` to the fundamental domain or just
    # slightly outside the fundamental domain, where the target tolerance
    # (not a strict bound) is specified by *one_minus_eps*. It also computes
    # the transformed value `w = gz`.
    # This function first tries to use
    # :func:`acb_modular_fundamental_domain_approx_d` and checks if the
    # result is acceptable. If this fails, it calls
    # :func:`acb_modular_fundamental_domain_approx_arf` with higher precision.
    # Finally, `w = gz` is evaluated by a single application of *g*.
    # The algorithm can fail to produce a satisfactory transformation.
    # The output *g* is always set to *some* correct modular transformation,
    # but it is up to the user to verify a posteriori that `w` is close enough
    # to the fundamental domain.

    int acb_modular_is_in_fundamental_domain(const acb_t z, const arf_t tol, slong prec)
    # Returns nonzero if it is certainly true that `|z| \ge 1 - \varepsilon` and
    # `|\operatorname{Re}(z)| \le 1/2 + \varepsilon` where `\varepsilon` is
    # specified by *tol*. Returns zero if this is false or cannot be determined.

    void acb_modular_fill_addseq(slong * tab, slong len)
    # Builds a near-optimal addition sequence for a sequence of integers
    # which is assumed to be reasonably dense.
    # As input, the caller should set each entry in *tab* to `-1` if
    # that index is to be part of the addition sequence, and to 0 otherwise.
    # On output, entry *i* in *tab* will either be zero (if the number is
    # not part of the sequence), or a value *j* such that both
    # *j* and `i - j` are also marked.
    # The first two entries in *tab* are ignored (the number 1 is always
    # assumed to be part of the sequence).

    void acb_modular_theta_transform(int * R, int * S, int * C, const psl2z_t g)
    # We wish to write a theta function with quasiperiod `\tau` in terms
    # of a theta function with quasiperiod `\tau' = g \tau`, given
    # some `g = (a, b; c, d) \in \text{PSL}(2, \mathbb{Z})`.
    # For `i = 0, 1, 2, 3`, this function computes integers `R_i` and `S_i`
    # (*R* and *S* should be arrays of length 4)
    # and `C \in \{0, 1\}` such that
    # .. math ::
    # \theta_{1+i}(z,\tau) = \exp(\pi i R_i / 4) \cdot A \cdot B \cdot \theta_{1+S_i}(z',\tau')
    # where `z' = z, A = B = 1` if `C = 0`, and
    # .. math ::
    # z' = \frac{-z}{c \tau + d}, \quad
    # A = \sqrt{\frac{i}{c \tau + d}}, \quad
    # B = \exp\left(-\pi i c \frac{z^2}{c \tau + d}\right)
    # if `C = 1`. Note that `A` is well-defined with the principal branch
    # of the square root since `A^2 = i/(c \tau + d)` lies in the right half-plane.
    # Firstly, if `c = 0`, we have
    # `\theta_i(z, \tau) = \exp(-\pi i b / 4) \theta_i(z, \tau+b)`
    # for `i = 1, 2`, whereas
    # `\theta_3` and `\theta_4` remain unchanged when `b` is even
    # and swap places with each other when `b` is odd.
    # In this case we set `C = 0`.
    # For an arbitrary `g` with `c > 0`, we set `C = 1`. The general
    # transformations are given by Rademacher [Rad1973]_.
    # We need the function `\theta_{m,n}(z,\tau)` defined for `m, n \in \mathbb{Z}` by
    # (beware of the typos in [Rad1973]_)
    # .. math ::
    # \theta_{0,0}(z,\tau) = \theta_3(z,\tau), \quad
    # \theta_{0,1}(z,\tau) = \theta_4(z,\tau)
    # .. math ::
    # \theta_{1,0}(z,\tau) = \theta_2(z,\tau), \quad
    # \theta_{1,1}(z,\tau) = i \theta_1(z,\tau)
    # .. math ::
    # \theta_{m+2,n}(z,\tau) = (-1)^n \theta_{m,n}(z,\tau)
    # .. math ::
    # \theta_{m,n+2}(z,\tau) = \theta_{m,n}(z,\tau).
    # Then we may write
    # .. math ::
    # \theta_1(z,\tau) &= \varepsilon_1 A B \theta_1(z', \tau')
    # \theta_2(z,\tau) &= \varepsilon_2 A B \theta_{1-c,1+a}(z', \tau')
    # \theta_3(z,\tau) &= \varepsilon_3 A B \theta_{1+d-c,1-b+a}(z', \tau')
    # \theta_4(z,\tau) &= \varepsilon_4 A B \theta_{1+d,1-b}(z', \tau')
    # where `\varepsilon_i` is an 8th root of unity.
    # Specifically, if we denote the 24th root of unity
    # in the transformation formula of the Dedekind eta
    # function by `\varepsilon(a,b,c,d) = \exp(\pi i R(a,b,c,d) / 12)`
    # (see :func:`acb_modular_epsilon_arg`), then:
    # .. math ::
    # \varepsilon_1(a,b,c,d) &= \exp(\pi i [R(-d,b,c,-a) + 1] / 4)
    # \varepsilon_2(a,b,c,d) &= \exp(\pi i [-R(a,b,c,d) + (5+(2-c)a)] / 4)
    # \varepsilon_3(a,b,c,d) &= \exp(\pi i [-R(a,b,c,d) + (4+(c-d-2)(b-a))] / 4)
    # \varepsilon_4(a,b,c,d) &= \exp(\pi i [-R(a,b,c,d) + (3-(2+d)b)] / 4)
    # These formulas are easily derived from the formulas in [Rad1973]_
    # (Rademacher has the transformed/untransformed variables exchanged,
    # and his "`\varepsilon`" differs from ours by a constant
    # offset in the phase).

    void acb_modular_addseq_theta(slong * exponents, slong * aindex, slong * bindex, slong num)
    # Constructs an addition sequence for the first *num* squares and triangular
    # numbers interleaved (excluding zero), i.e. 1, 2, 4, 6, 9, 12, 16, 20, 25, 30 etc.

    void acb_modular_theta_sum(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t w, int w_is_unit, const acb_t q, slong len, slong prec)
    # Simultaneously computes the first *len* coefficients of each of the
    # formal power series
    # .. math ::
    # \theta_1(z+x,\tau) / q_{1/4} \in \mathbb{C}[[x]]
    # \theta_2(z+x,\tau) / q_{1/4} \in \mathbb{C}[[x]]
    # \theta_3(z+x,\tau) \in \mathbb{C}[[x]]
    # \theta_4(z+x,\tau) \in \mathbb{C}[[x]]
    # given `w = \exp(\pi i z)` and `q = \exp(\pi i \tau)`, by summing
    # a finite truncation of the respective theta function series.
    # In particular, with *len* equal to 1, computes the respective
    # value of the theta function at the point *z*.
    # We require *len* to be positive.
    # If *w_is_unit* is nonzero, *w* is assumed to lie on the unit circle,
    # i.e. *z* is assumed to be real.
    # Note that the factor `q_{1/4}` is removed from `\theta_1` and `\theta_2`.
    # To get the true theta function values, the user has to multiply
    # this factor back. This convention avoids unnecessary computations,
    # since the user can compute `q_{1/4} = \exp(\pi i \tau / 4)` followed by
    # `q = (q_{1/4})^4`, and in many cases when computing products or quotients
    # of theta functions, the factor `q_{1/4}` can be eliminated entirely.
    # This function is intended for `|q| \ll 1`. It can be called with any
    # `q`, but will return useless intervals if convergence is not rapid.
    # For general evaluation of theta functions, the user should only call
    # this function after applying a suitable modular transformation.
    # We consider the sums together, alternatingly updating `(\theta_1, \theta_2)`
    # or `(\theta_3, \theta_4)`. For `k = 0, 1, 2, \ldots`, the powers of `q`
    # are `\lfloor (k+2)^2 / 4 \rfloor = 1, 2, 4, 6, 9` etc. and the powers of `w` are
    # `\pm (k+2) = \pm 2, \pm 3, \pm 4, \ldots` etc. The scheme
    # is illustrated by the following table:
    # .. math ::
    # \begin{array}{llll}
    # & \theta_1, \theta_2 & q^0 & (w^1 \pm w^{-1}) \\
    # k = 0  & \theta_3, \theta_4 & q^1 & (w^2 \pm w^{-2}) \\
    # k = 1  & \theta_1, \theta_2 & q^2 & (w^3 \pm w^{-3}) \\
    # k = 2  & \theta_3, \theta_4 & q^4 & (w^4 \pm w^{-4}) \\
    # k = 3  & \theta_1, \theta_2 & q^6 & (w^5 \pm w^{-5}) \\
    # k = 4  & \theta_3, \theta_4 & q^9 & (w^6 \pm w^{-6}) \\
    # k = 5  & \theta_1, \theta_2 & q^{12} & (w^7 \pm w^{-7}) \\
    # \end{array}
    # For some integer `N \ge 1`, the summation is stopped just before term
    # `k = N`. Let `Q = |q|`, `W = \max(|w|,|w^{-1}|)`,
    # `E = \lfloor (N+2)^2 / 4 \rfloor` and
    # `F = \lfloor (N+1)/2 \rfloor + 1`. The error of the
    # zeroth derivative can be bounded as
    # .. math ::
    # 2 Q^E W^{N+2} \left[ 1 + Q^F W + Q^{2F} W^2 + \ldots \right]
    # = \frac{2 Q^E W^{N+2}}{1 - Q^F W}
    # provided that the denominator is positive (otherwise we set
    # the error bound to infinity).
    # When *len* is greater than 1, consider the derivative of order *r*.
    # The term of index *k* and order *r* picks up a factor of magnitude
    # `(k+2)^r` from differentiation of `w^{k+2}` (it also picks up a factor
    # `\pi^r`, but we omit this until we rescale the coefficients
    # at the end of the computation). Thus we have the error bound
    # .. math ::
    # 2 Q^E W^{N+2} (N+2)^r \left[ 1 + Q^F W \frac{(N+3)^r}{(N+2)^r} + Q^{2F} W^2 \frac{(N+4)^r}{(N+2)^r} + \ldots \right]
    # which by the inequality `(1 + m/(N+2))^r \le \exp(mr/(N+2))`
    # can be bounded as
    # .. math ::
    # \frac{2 Q^E W^{N+2} (N+2)^r}{1 - Q^F W \exp(r/(N+2))},
    # again valid when the denominator is positive.
    # To actually evaluate the series, we write the even
    # cosine terms as `w^{2n} + w^{-2n}`, the odd cosine terms as
    # `w (w^{2n} + w^{-2n-2})`, and the sine terms as `w (w^{2n} - w^{-2n-2})`.
    # This way we only need even powers of `w` and `w^{-1}`.
    # The implementation is not yet optimized for real `z`, in which case
    # further work can be saved.
    # This function does not permit aliasing between input and output
    # arguments.

    void acb_modular_theta_const_sum_basecase(acb_t theta2, acb_t theta3, acb_t theta4, const acb_t q, slong N, slong prec)

    void acb_modular_theta_const_sum_rs(acb_t theta2, acb_t theta3, acb_t theta4, const acb_t q, slong N, slong prec)
    # Computes the truncated theta constant sums
    # `\theta_2 = \sum_{k(k+1) < N} q^{k(k+1)}`,
    # `\theta_3 = \sum_{k^2 < N} q^{k^2}`,
    # `\theta_4 = \sum_{k^2 < N} (-1)^k q^{k^2}`.
    # The *basecase* version uses a short addition sequence.
    # The *rs* version uses rectangular splitting.
    # The algorithms are described in [EHJ2016]_.

    void acb_modular_theta_const_sum(acb_t theta2, acb_t theta3, acb_t theta4, const acb_t q, slong prec)
    # Computes the respective theta constants by direct summation
    # (without applying modular transformations). This function
    # selects an appropriate *N*, calls either
    # :func:`acb_modular_theta_const_sum_basecase` or
    # :func:`acb_modular_theta_const_sum_rs` or depending on *N*,
    # and adds a bound for the truncation error.

    void acb_modular_theta_notransform(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, slong prec)
    # Evaluates the Jacobi theta functions `\theta_i(z,\tau)`, `i = 1, 2, 3, 4`
    # simultaneously. This function does not move `\tau` to the fundamental domain.
    # This is generally worse than :func:`acb_modular_theta`, but can
    # be slightly better for moderate input.

    void acb_modular_theta(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, slong prec)
    # Evaluates the Jacobi theta functions `\theta_i(z,\tau)`, `i = 1, 2, 3, 4`
    # simultaneously. This function moves `\tau` to the fundamental domain
    # and then also reduces `z` modulo `\tau`
    # before calling :func:`acb_modular_theta_sum`.

    void acb_modular_theta_jet_notransform(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau, slong len, slong prec)

    void acb_modular_theta_jet(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau, slong len, slong prec)
    # Evaluates the Jacobi theta functions along with their derivatives
    # with respect to *z*, writing the first *len* coefficients in the power
    # series `\theta_i(z+x,\tau) \in \mathbb{C}[[x]]` to
    # each respective output variable. The *notransform* version does not
    # move `\tau` to the fundamental domain or reduce `z` during the computation.

    void _acb_modular_theta_series(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec)

    void acb_modular_theta_series(acb_poly_t theta1, acb_poly_t theta2, acb_poly_t theta3, acb_poly_t theta4, const acb_poly_t z, const acb_t tau, slong len, slong prec)
    # Evaluates the respective Jacobi theta functions of the power series *z*,
    # truncated to length *len*. Either of the output variables can be *NULL*.

    void acb_modular_addseq_eta(slong * exponents, slong * aindex, slong * bindex, slong num)
    # Constructs an addition sequence for the first *num* generalized pentagonal
    # numbers (excluding zero), i.e. 1, 2, 5, 7, 12, 15, 22, 26, 35, 40 etc.

    void acb_modular_eta_sum(acb_t eta, const acb_t q, slong prec)
    # Evaluates the Dedekind eta function
    # without the leading 24th root, i.e.
    # .. math :: \exp(-\pi i \tau/12) \eta(\tau) = \sum_{n=-\infty}^{\infty} (-1)^n q^{(3n^2-n)/2}
    # given `q = \exp(2 \pi i \tau)`, by summing the defining series.
    # This function is intended for `|q| \ll 1`. It can be called with any
    # `q`, but will return useless intervals if convergence is not rapid.
    # For general evaluation of the eta function, the user should only call
    # this function after applying a suitable modular transformation.
    # The series is evaluated using either a short addition sequence or
    # rectangular splitting, depending on the number of terms.
    # The algorithms are described in [EHJ2016]_.

    int acb_modular_epsilon_arg(const psl2z_t g)
    # Given `g = (a, b; c, d)`, computes an integer `R` such that
    # `\varepsilon(a,b,c,d) = \exp(\pi i R / 12)` is the 24th root of unity in
    # the transformation formula for the Dedekind eta function,
    # .. math ::
    # \eta\left(\frac{a\tau+b}{c\tau+d}\right) = \varepsilon (a,b,c,d)
    # \sqrt{c\tau+d} \eta(\tau).

    void acb_modular_eta(acb_t r, const acb_t tau, slong prec)
    # Computes the Dedekind eta function `\eta(\tau)` given `\tau` in the upper
    # half-plane. This function applies the functional equation to move
    # `\tau` to the fundamental domain before calling
    # :func:`acb_modular_eta_sum`.

    void acb_modular_j(acb_t r, const acb_t tau, slong prec)
    # Computes Klein's j-invariant `j(\tau)` given `\tau` in the upper
    # half-plane. The function is normalized so that `j(i) = 1728`.
    # We first move `\tau` to the fundamental domain, which does not change
    # the value of the function. Then we use the formula
    # `j(\tau) = 32 (\theta_2^8+\theta_3^8+\theta_4^8)^3 / (\theta_2 \theta_3 \theta_4)^8` where
    # `\theta_i = \theta_i(0,\tau)`.

    void acb_modular_lambda(acb_t r, const acb_t tau, slong prec)
    # Computes the lambda function
    # `\lambda(\tau) = \theta_2^4(0,\tau) / \theta_3^4(0,\tau)`, which
    # is invariant under modular transformations `(a, b; c, d)`
    # where `a, d` are odd and `b, c` are even.

    void acb_modular_delta(acb_t r, const acb_t tau, slong prec)
    # Computes the modular discriminant `\Delta(\tau) = \eta(\tau)^{24}`,
    # which transforms as
    # .. math ::
    # \Delta\left(\frac{a\tau+b}{c\tau+d}\right) = (c\tau+d)^{12} \Delta(\tau).
    # The modular discriminant is sometimes defined with an extra factor
    # `(2\pi)^{12}`, which we omit in this implementation.

    void acb_modular_eisenstein(acb_ptr r, const acb_t tau, slong len, slong prec)
    # Computes simultaneously the first *len* entries in the sequence
    # of Eisenstein series `G_4(\tau), G_6(\tau), G_8(\tau), \ldots`,
    # defined by
    # .. math ::
    # G_{2k}(\tau) = \sum_{m^2 + n^2 \ne 0} \frac{1}{(m+n\tau )^{2k}}
    # and satisfying
    # .. math ::
    # G_{2k} \left(\frac{a\tau+b}{c\tau+d}\right) = (c\tau+d)^{2k} G_{2k}(\tau).
    # We first evaluate `G_4(\tau)` and `G_6(\tau)` on the fundamental
    # domain using theta functions, and then compute the Eisenstein series
    # of higher index using a recurrence relation.

    void acb_modular_elliptic_k(acb_t w, const acb_t m, slong prec)

    void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, slong len, slong prec)

    void acb_modular_elliptic_e(acb_t w, const acb_t m, slong prec)

    void acb_modular_elliptic_p(acb_t wp, const acb_t z, const acb_t tau, slong prec)

    void acb_modular_elliptic_p_zpx(acb_ptr wp, const acb_t z, const acb_t tau, slong len, slong prec)

    void acb_modular_hilbert_class_poly(fmpz_poly_t res, slong D)
    # Sets *res* to the Hilbert class polynomial of discriminant *D*,
    # defined as
    # .. math ::
    # H_D(x) = \prod_{(a,b,c)} \left(x - j\left(\frac{-b+\sqrt{D}}{2a}\right)\right)
    # where `(a,b,c)` ranges over the primitive reduced positive
    # definite binary quadratic forms of discriminant `b^2 - 4ac = D`.
    # The Hilbert class polynomial is only defined if `D < 0` and *D*
    # is congruent to 0 or 1 mod 4. If some other value of *D* is passed as
    # input, *res* is set to the zero polynomial.
