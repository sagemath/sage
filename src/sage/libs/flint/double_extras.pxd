# distutils: libraries = flint
# distutils: depends = flint/double_extras.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    double d_randtest(flint_rand_t state)
    # Returns a random number in the interval `[0.5, 1)`.

    double d_randtest_signed(flint_rand_t state, long minexp, long maxexp)
    # Returns a random signed number with exponent between ``minexp`` and
    # ``maxexp`` or zero.

    double d_randtest_special(flint_rand_t state, long minexp, long maxexp)
    # Returns a random signed number with exponent between ``minexp`` and
    # ``maxexp``, zero, ``D_NAN`` or `\pm`\ ``D_INF``.

    double d_polyval(const double * poly, int len, double x)
    # Uses Horner's rule to evaluate the polynomial defined by the given
    # ``len`` coefficients. Requires that ``len`` is nonzero.

    double d_lambertw(double x)
    # Computes the principal branch of the Lambert W function, solving
    # the equation `x = W(x) \exp(W(x))`. If `x < -1/e`, the solution is
    # complex, and NaN is returned.
    # Depending on the magnitude of `x`, we start from a piecewise rational
    # approximation or a zeroth-order truncation of the asymptotic expansion
    # at infinity, and perform 0, 1 or 2 iterations with Halley's
    # method to obtain full accuracy.
    # A test of `10^7` random inputs showed a maximum relative error smaller
    # than 0.95 times ``DBL_EPSILON`` (`2^{-52}`) for positive `x`.
    # Accuracy for negative `x` is slightly worse, and can grow to
    # about 10 times ``DBL_EPSILON`` close to `-1/e`.
    # However, accuracy may be worse depending on compiler flags and
    # the accuracy of the system libm functions.

    int d_is_nan(double x)
    # Returns a nonzero integral value if ``x`` is ``D_NAN``, and otherwise
    # returns 0.

    double d_log2(double x)
    # Returns the base 2 logarithm of ``x`` provided ``x`` is positive. If
    # a domain or pole error occurs, the appropriate error value is returned.
