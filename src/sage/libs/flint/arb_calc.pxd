# distutils: libraries = flint
# distutils: depends = flint/arb_calc.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void arf_interval_init(arf_interval_t v)

    void arf_interval_clear(arf_interval_t v)

    arf_interval_ptr _arf_interval_vec_init(long n)

    void _arf_interval_vec_clear(arf_interval_ptr v, long n)

    void arf_interval_set(arf_interval_t v, const arf_interval_t u)

    void arf_interval_swap(arf_interval_t v, arf_interval_t u)

    void arf_interval_get_arb(arb_t x, const arf_interval_t v, long prec)

    void arf_interval_printd(const arf_interval_t v, long n)
    # Helper functions for endpoint-based intervals.

    void arf_interval_fprintd(FILE * file, const arf_interval_t v, long n)
    # Helper functions for endpoint-based intervals.

    long arb_calc_isolate_roots(arf_interval_ptr * found, int ** flags, arb_calc_func_t func, void * param, const arf_interval_t interval, long maxdepth, long maxeval, long maxfound, long prec)
    # Rigorously isolates single roots of a real analytic function
    # on the interior of an interval.
    # This routine writes an array of *n* interesting subintervals of
    # *interval* to *found* and corresponding flags to *flags*, returning the integer *n*.
    # The output has the following properties:
    # * The function has no roots on *interval* outside of the output
    # subintervals.
    # * Subintervals are sorted in increasing order (with no overlap except
    # possibly starting and ending with the same point).
    # * Subintervals with a flag of 1 contain exactly one (single) root.
    # * Subintervals with any other flag may or may not contain roots.
    # If no flags other than 1 occur, all roots of the function on *interval*
    # have been isolated. If there are output subintervals on which the
    # existence or nonexistence of roots could not be determined,
    # the user may attempt further searches on those subintervals
    # (possibly with increased precision and/or increased
    # bounds for the breaking criteria). Note that roots of multiplicity
    # higher than one and roots located exactly at endpoints cannot be isolated
    # by the algorithm.
    # The following breaking criteria are implemented:
    # * At most *maxdepth* recursive subdivisions are attempted. The smallest
    # details that can be distinguished are therefore about
    # `2^{-\text{maxdepth}}` times the width of *interval*.
    # A typical, reasonable value might be between 20 and 50.
    # * If the total number of tested subintervals exceeds *maxeval*, the
    # algorithm is terminated and any untested subintervals are added
    # to the output. The total number of calls to *func* is thereby restricted
    # to a small multiple of *maxeval* (the actual count can be slightly
    # higher depending on implementation details).
    # A typical, reasonable value might be between 100 and 100000.
    # * The algorithm terminates if *maxfound* roots have been isolated.
    # In particular, setting *maxfound* to 1 can be used to locate
    # just one root of the function even if there are numerous roots.
    # To try to find all roots, *LONG_MAX* may be passed.
    # The argument *prec* denotes the precision used to evaluate the
    # function. It is possibly also used for some other arithmetic operations
    # performed internally by the algorithm. Note that it probably does not
    # make sense for *maxdepth* to exceed *prec*.
    # Warning: it is assumed that subdivision points of *interval* can be
    # represented exactly as floating-point numbers in memory.
    # Do not pass `1 \pm 2^{-10^{100}}` as input.

    int arb_calc_refine_root_bisect(arf_interval_t r, arb_calc_func_t func, void * param, const arf_interval_t start, long it, long prec)
    # Given an interval *start* known to contain a single root of *func*,
    # refines it using *iter* bisection steps. The algorithm can
    # return a failure code if the sign of the function at an evaluation
    # point is ambiguous. The output *r* is set to a valid isolating interval
    # (possibly just *start*) even if the algorithm fails.

    void arb_calc_newton_conv_factor(arf_t conv_factor, arb_calc_func_t func, void * param, const arb_t conv_region, long prec)
    # Given an interval `I` specified by *conv_region*, evaluates a bound
    # for `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    # where `f` is the function specified by *func* and *param*.
    # The bound is obtained by evaluating `f'(I)` and `f''(I)` directly.
    # If `f` is ill-conditioned, `I` may need to be extremely precise in
    # order to get an effective, finite bound for *C*.

    int arb_calc_newton_step(arb_t xnew, arb_calc_func_t func, void * param, const arb_t x, const arb_t conv_region, const arf_t conv_factor, long prec)
    # Performs a single step with an interval version of Newton's method.
    # The input consists of the function `f` specified
    # by *func* and *param*, a ball `x = [m-r, m+r]` known
    # to contain a single root of `f`, a ball `I` (*conv_region*)
    # containing `x` with an associated bound (*conv_factor*) for
    # `C = \sup_{t,u \in I} \frac{1}{2} |f''(t)| / |f'(u)|`,
    # and a working precision *prec*.
    # The Newton update consists of setting
    # `x' = [m'-r', m'+r']` where `m' = m - f(m) / f'(m)`
    # and `r' = C r^2`. The expression `m - f(m) / f'(m)` is evaluated
    # using ball arithmetic at a working precision of *prec* bits, and the
    # rounding error during this evaluation is accounted for in the output.
    # We now check that `x' \in I` and `r' < r`. If both conditions are
    # satisfied, we set *xnew* to `x'` and return *ARB_CALC_SUCCESS*.
    # If either condition fails, we set *xnew* to `x` and return
    # *ARB_CALC_NO_CONVERGENCE*, indicating that no progress
    # is made.

    int arb_calc_refine_root_newton(arb_t r, arb_calc_func_t func, void * param, const arb_t start, const arb_t conv_region, const arf_t conv_factor, long eval_extra_prec, long prec)
    # Refines a precise estimate of a single root of a function
    # to high precision by performing several Newton steps, using
    # nearly optimally chosen doubling precision steps.
    # The inputs are defined as for *arb_calc_newton_step*, except for
    # the precision parameters: *prec* is the target accuracy and
    # *eval_extra_prec* is the estimated number of guard bits that need
    # to be added to evaluate the function accurately close to the root
    # (for example, if the function is a polynomial with large coefficients
    # of alternating signs and Horner's rule is used to evaluate it,
    # the extra precision should typically be approximately
    # the bit size of the coefficients).
    # This function returns *ARB_CALC_SUCCESS* if all attempted
    # Newton steps are successful (note that this does not guarantee
    # that the computed root is accurate to *prec* bits, which has
    # to be verified by the user), only that it is more accurate
    # than the starting ball.
    # On failure, *ARB_CALC_IMPRECISE_INPUT*
    # or *ARB_CALC_NO_CONVERGENCE* may be returned. In this case, *r*
    # is set to a ball for the root which is valid but likely
    # does have full accuracy (it can possibly just be equal
    # to the starting ball).
