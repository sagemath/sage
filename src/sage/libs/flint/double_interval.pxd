# distutils: libraries = flint
# distutils: depends = flint/double_interval.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    di_t di_interval(double a, double b)
    # Returns the interval `[a, b]`. We require that the endpoints
    # are ordered and not NaN.

    di_t arb_get_di(const arb_t x)
    # Returns the ball *x* converted to a double-precision interval.

    void arb_set_di(arb_t res, di_t x, long prec)
    # Sets the ball *res* to the double-precision interval *x*,
    # rounded to *prec* bits.

    void di_print(di_t x)
    # Prints *x* to standard output. This simply prints decimal
    # representations of the floating-point endpoints; the
    # decimals are not guaranteed to be rounded outward.

    double d_randtest2(flint_rand_t state)
    # Returns a random non-NaN ``double`` with any exponent.
    # The value can be infinite or subnormal.

    di_t di_randtest(flint_rand_t state)
    # Returns an interval with random endpoints.

    di_t di_neg(di_t x)
    # Returns the exact negation of *x*.

    di_t di_fast_add(di_t x, di_t y)
    di_t di_fast_sub(di_t x, di_t y)
    di_t di_fast_mul(di_t x, di_t y)
    di_t di_fast_div(di_t x, di_t y)
    # Returns the sum, difference, product or quotient of *x* and *y*.
    # Division by zero is currently defined to return `[-\infty, +\infty]`.

    di_t di_fast_sqr(di_t x)
    # Returns the square of *x*. The output is clamped to
    # be nonnegative.

    di_t di_fast_add_d(di_t x, double y)
    di_t di_fast_sub_d(di_t x, double y)
    di_t di_fast_mul_d(di_t x, double y)
    di_t di_fast_div_d(di_t x, double y)
    # Arithmetic with an exact ``double`` operand.

    di_t di_fast_log_nonnegative(di_t x)
    # Returns an enclosure of `\log(x)`. The lower endpoint of *x*
    # is rounded up to 0 if it is negative.

    di_t di_fast_mid(di_t x)
    # Returns an enclosure of the midpoint of *x*.

    double di_fast_ubound_radius(di_t x)
    # Returns an upper bound for the radius of *x*.
