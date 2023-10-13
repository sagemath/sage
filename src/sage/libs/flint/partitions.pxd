# distutils: libraries = flint
# distutils: depends = flint/partitions.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void partitions_rademacher_bound(arf_t b, const fmpz_t n, unsigned long N)
    # Sets `b` to an upper bound for
    # .. math ::
    # M(n,N) = \frac{44 \pi^2}{225 \sqrt 3} N^{-1/2}
    # + \frac{\pi \sqrt{2}}{75} \left( \frac{N}{n-1} \right)^{1/2}
    # \sinh\left(\frac{\pi}{N} \sqrt{\frac{2n}{3}}\right).
    # This formula gives an upper bound for the truncation error in the
    # Hardy-Ramanujan-Rademacher formula when the series is taken up
    # to the term `t(n,N)` inclusive.

    void partitions_hrr_sum_arb(arb_t x, const fmpz_t n, long N0, long N, int use_doubles)
    # Evaluates the partial sum `\sum_{k=N_0}^N t(n,k)` of the
    # Hardy-Ramanujan-Rademacher series.
    # If *use_doubles* is nonzero, doubles and the system's standard library math
    # functions are used to evaluate the smallest terms. This significantly
    # speeds up evaluation for small `n` (e.g. `n < 10^6`), and gives a small speed
    # improvement for larger `n`, but the result is not guaranteed to be correct.
    # In practice, the error is estimated very conservatively, and unless
    # the system's standard library is broken, use of doubles can be considered
    # safe. Setting *use_doubles* to zero gives a fully guaranteed
    # bound.

    void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)
    # Computes the partition function `p(n)` using the Hardy-Ramanujan-Rademacher
    # formula. This function computes a numerical ball containing `p(n)`
    # and verifies that the ball contains a unique integer.
    # If *n* is sufficiently large and a number of threads greater than 1
    # has been selected with :func:`flint_set_num_threads()`, the computation
    # time will be reduced by using two threads.
    # See :func:`partitions_hrr_sum_arb` for an explanation of the
    # *use_doubles* option.

    void partitions_fmpz_ui(fmpz_t p, unsigned long n)
    # Computes the partition function `p(n)` using the Hardy-Ramanujan-Rademacher
    # formula. This function computes a numerical ball containing `p(n)`
    # and verifies that the ball contains a unique integer.

    void partitions_fmpz_ui_using_doubles(fmpz_t p, unsigned long n)
    # Computes the partition function `p(n)`, enabling the use of doubles
    # internally. This significantly speeds up evaluation for small `n`
    # (e.g. `n < 10^6`), but the error bounds are not certified
    # (see remarks for :func:`partitions_hrr_sum_arb`).

    void partitions_leading_fmpz(arb_t res, const fmpz_t n, long prec)
    # Sets *res* to the leading term in the Hardy-Ramanujan series
    # for `p(n)` (without Rademacher's correction of this term, which is
    # vanishingly small when `n` is large), that is,
    # `\sqrt{12} (1-1/t) e^t / (24n-1)` where `t = \pi \sqrt{24n-1} / 6`.
