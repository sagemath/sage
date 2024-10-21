# sage_setup: distribution = sagemath-flint
# distutils: libraries = gmp
"""
Enumeration of totally real fields: data

AUTHORS:

- John Voight (2007-09-01): Initial version
- John Voight (2007-09-19): various optimization tweaks
- John Voight (2007-10-09): improvements: Smyth bound, Lagrange multipliers for b
- Craig Citro and John Voight (2007-11-04): type checking and other polishing
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein and John Voight
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.math cimport sqrt
from cysignals.memory cimport sig_malloc, sig_free

from sage.arith.misc import binomial
from sage.arith.misc import GCD as gcd
from sage.libs.gmp.mpz cimport *
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.real_mpfr import RealField
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer

# Other global variables
ZZx = PolynomialRing(ZZ, 'x')

from libc.math cimport lrint, floor, ceil, fabs, round


#*********************************************************************
# Auxiliary routine
# Hermite constant, naive Newton-Raphson, and a specialized Lagrange
# multiplier solver.
#*********************************************************************

def hermite_constant(n):
    r"""
    Return the `n`-th Hermite constant.

    The `n`-th Hermite constant (typically denoted `\gamma_n`), is defined
    to be

    .. MATH::

        \max_L \min_{0 \neq x \in L} ||x||^2

    where `L` runs over all lattices of dimension `n` and determinant `1`.

    For `n \leq 8` it returns the exact value of `\gamma_n`, and for
    `n > 9` it returns an upper bound on `\gamma_n`.

    INPUT:

    - ``n`` -- integer

    OUTPUT:

    (an upper bound for) the Hermite constant `\gamma_n`

    EXAMPLES::

        sage: hermite_constant(1) # trivial one-dimensional lattice
        1.0
        sage: hermite_constant(2) # Eisenstein lattice
        1.1547005383792515
        sage: 2/sqrt(3.)
        1.15470053837925
        sage: hermite_constant(8) # E_8
        2.0

    .. NOTE::

        The upper bounds used can be found in [CS1999]_ and [CE2003]_.

    AUTHORS:

    - John Voight (2007-09-03)
    """

    if n <= 8:
        # Exact values are known for gamma_n.
        gamman = [1, 1, 4./3, 2, 4, 8, 64./3, 64, 256][n]
        gamma = gamman**(1./n)
    elif n <= 36:
        gamma = [2.13263235569928, 2.26363016185702, 2.39334691240146,
                 2.52178702088414, 2.64929462619823, 2.77580405570023,
                 2.90147761892077, 3.02639364467182, 3.15067928476872,
                 3.27433066745617, 3.39744386110070, 3.52006195697466,
                 3.64224310140724, 3.76403701226104, 3.88547626036618,
                 4.00659977840648, 4.12744375027069, 4.24804458298177,
                 4.36843113799634, 4.48863097933934, 4.60866759008263,
                 4.72856660611662, 4.84834821242630, 4.96803435811402,
                 5.08764086822471, 5.20718687262715, 5.32668836123079,
                 5.44615801810606][n-9]
    else:
        # Mordell's inequality.
        gamma = 5.44615801810606**((n-1.)/35)

    return gamma


cdef double eval_seq_as_poly(int *f, int n, double x) noexcept:
    r"""
    Evaluates the sequence a, thought of as a polynomial with

    .. MATH::

        f[n]*x^n + f[n-1]*x^(n-1) + ... + f[0].
    """
    cdef double s

    # Horner's method: With polynomials of small degree, we shouldn't
    # expect asymptotic methods to be any faster.
    s = f[n]
    for i from n > i >= 0:
        s = s * x + f[i]
    return s

cdef double newton(int *f, int *df, int n, double x0, double eps) noexcept:
    r"""
    Find the real root x of f (with derivative df) near x0
    with provable precision eps, i.e. |x-z| < eps where z is the actual
    root.
    The sequence a corresponds to the polynomial f with

    .. MATH::

        f(x) = x^n + a[n-1]*x^(n-1) + ... + a[0].

    This function (for speed reasons) has no error checking and no
    guarantees are made as to the convergence; a naive Newton-Raphson
    method is used.
    """
    cdef double x, rdx, dx, fx

    x = x0
    dx = eval_seq_as_poly(f,n,x)/eval_seq_as_poly(df,n-1,x)
    x -= dx
    while fabs(dx) > eps:
        # In truly optimized code, one could tune by automatically
        # iterating a certain number of times based on the size of dx to
        # save on a few comparisons.
        # This savings should be almost negligible...?
        dx = eval_seq_as_poly(f,n,x)/eval_seq_as_poly(df,n-1,x)
        x -= dx

    # Small hack for improved performance elsewhere: if it is close to an
    # integer, give it full precision as an integer.
    rdx = round(x)
    if fabs(rdx-x) < eps:
        x = rdx

    # Now ensure that either f(x-eps) or f(x+eps) has opposite sign
    # as f(x), which implies that |x-z| < eps.
    fx = eval_seq_as_poly(f,n,x)
    while not (fx == 0 or fx*eval_seq_as_poly(f,n,x+eps) < 0 or
                          fx*eval_seq_as_poly(f,n,x-eps) < 0):
        dx = eval_seq_as_poly(f,n,x)/eval_seq_as_poly(df,n-1,x)
        x -= dx
        fx = eval_seq_as_poly(f,n,x)
    return x

cdef void newton_in_intervals(int *f, int *df, int n, double *beta,
                              double eps, double *rts) noexcept:
    r"""
    Find the real roots of f in the intervals specified by beta:

        (beta[0],beta[1]), (beta[1],beta[2]), ..., (beta[n-1], beta[n])

    Calls newton_kernel, so same provisos apply---in particular,
    each interval should contain a unique simple root.
    Note the derivative *df is passed but is recomputed--this is
    just a way to save a malloc and free for each call.
    """
    cdef int i

    for i from 0 <= i < n:
        df[i] = f[i+1]*(i+1)
    for i from 0 <= i < n:
        rts[i] = newton(f, df, n, (beta[i]+beta[i+1])/2, eps)

cpdef lagrange_degree_3(int n, int an1, int an2, int an3):
    r"""
    Private function.  Solves the equations which arise in the Lagrange multiplier
    for degree 3: for each `1 \leq r \leq n-2`, we solve

    .. MATH::

        r*x^i + (n-1-r)\cdot y^i + z^i = s_i \quad (i = 1,2,3)

    where the `s_i` are the power sums determined by the coefficients `a`.
    We output the largest value of `z` which occurs.
    We use a precomputed elimination ideal.

    EXAMPLES::

        sage: ls = sage.rings.number_field.totallyreal_data.lagrange_degree_3(3,0,1,2)
        sage: [RealField(10)(x) for x in ls]
        [-1.0, -1.0]
        sage: sage.rings.number_field.totallyreal_data.lagrange_degree_3(3,6,1,2) # random
        [-5.8878, -5.8878]

    TESTS:

    Check that :issue:`13101` is solved::

        sage: sage.rings.number_field.totallyreal_data.lagrange_degree_3(4,12,19,42)
        [0.0, -1.0]
    """
    cdef long coeffs[7]
    cdef int r, rsq, rcu
    cdef int nr, nrsq, nrcu
    cdef int s1, s1sq, s1cu, s1fo, s2, s2sq, s2cu, s3, s3sq

    RRx = PolynomialRing(RealField(20),'x')

    # Newton's relations.
    s1 = -an1
    s2 = -(an1*s1 + 2*an2)
    s3 = -(an1*s2 + an2*s1 + 3*an3)

    s1sq = s1*s1
    s1cu = s1*s1sq
    s1fo = s1*s1cu
    s2sq = s2*s2
    s2cu = s2*s2sq
    s3sq = s3*s3

    z4minmax = []

    for r from 1 <= r <= n-2:
        nr = n-1-r
        # common subexpressions
        rsq = r*r
        rcu = r*rsq
        nrsq = nr*nr
        nrcu = nr*nrsq

        ## x^6
        coeffs[6] = rcu*nr + rcu + 2*rsq*nrsq + 5*rsq*nr + 3*rsq + \
                    r*nrcu + 5*r*nrsq + 7*r*nr + 3*r + nrcu + \
                    3*nrsq + 3*nr + 1

        ## x^5
        coeffs[5] = -6*rsq*nr*s1 - 6*rsq*s1 - 6*r*nrsq*s1 - 18*r*nr*s1 - \
                    12*r*s1 - 6*nrsq*s1 - 12*nr*s1 - 6*s1

        ## x^4
        coeffs[4] = -3*rcu*s2 - 3*rsq*nr*s2 + 3*rsq*s1sq - 6*rsq*s2 - \
                    3*r*nrsq*s2 + 15*r*nr*s1sq - 6*r*nr*s2 + 18*r*s1sq - \
                    3*r*s2 - 3*nrcu*s2 + 3*nrsq*s1sq - 6*nrsq*s2 + \
                    18*nr*s1sq - 3*nr*s2 + 15*s1sq

        ## x^3
        coeffs[3] = -2*rcu*nr*s3 - 4*rsq*nrsq*s3 + 6*rsq*nr*s1*s2 - \
                    6*rsq*nr*s3 + 12*rsq*s1*s2 - 2*r*nrcu*s3 + \
                    6*r*nrsq*s1*s2 - 6*r*nrsq*s3 - 4*r*nr*s1cu + \
                    12*r*nr*s1*s2 - 4*r*nr*s3 - 12*r*s1cu + 12*r*s1*s2 + \
                    12*nrsq*s1*s2 - 12*nr*s1cu + 12*nr*s1*s2 - \
                    20*s1cu

        ## x^2
        coeffs[2] = 3*rcu*s2sq + 6*rsq*nr*s1*s3 - 3*rsq*nr*s2sq - \
                    6*rsq*s1sq*s2 + 3*rsq*s2sq + 6*r*nrsq*s1*s3 - \
                    3*r*nrsq*s2sq - 6*r*nr*s1sq*s2 + 12*r*nr*s1*s3 + \
                    3*r*nr*s2sq + 3*r*s1fo - 18*r*s1sq*s2 + \
                    3*nrcu*s2sq - 6*nrsq*s1sq*s2 + 3*nrsq*s2sq + \
                    3*nr*s1fo - 18*nr*s1sq*s2 + 15*s1fo

        ## x^1
        coeffs[1] = 6*rsq*nr*s2*s3 - 6*rsq*s1*s2sq + 6*r*nrsq*s2*s3 - \
                    12*r*nr*s1sq*s3 - 6*r*nr*s1*s2sq + 12*r*s1cu*s2 - \
                    6*nrsq*s1*s2sq + 12*nr*s1cu*s2 - 6*s1*s1fo

        ## x^0
        coeffs[0] = rcu*nr*s3sq - rcu*s2cu + 2*rsq*nrsq*s3sq - \
                    6*rsq*nr*s1*s2*s3 + rsq*nr*s2cu + 3*rsq*s1sq*s2sq + \
                    r*nrcu*s3sq - 6*r*nrsq*s1*s2*s3 + r*nrsq*s2cu + \
                    4*r*nr*s1cu*s3 + 3*r*nr*s1sq*s2sq - \
                    3*r*s1fo*s2 - nrcu*s2cu + \
                    3*nrsq*s1sq*s2sq - 3*nr*s1fo*s2 + \
                    s1sq*s1fo

        fcoeff = [ int(coeffs[i]) for i in range(7) ]
        f = ZZx(fcoeff)
        df = ZZx([i*coeffs[i] for i in range(1,7)])
        f = f//gcd(f, df)
        fcoeff = [int(c) for c in f.list()]

        rts = RRx(fcoeff).roots()

        if rts:
            rts = [rtsi[0] for rtsi in rts]
            z4minmax = [min(rts + z4minmax), max(rts + z4minmax)]

    if not z4minmax:
        return [0.0, -1.0]

    return z4minmax

cdef int __len_primes = 46
cdef long primessq[46]
primessq_py = [4, 9, 25, 49, 121, 169, 289, 361, 529, 841, 961, 1369, 1681, 1849, 2209, 2809, 3481, 3721, 4489, 5041, 5329, 6241, 6889, 7921, 9409, 10201, 10609, 11449, 11881, 12769, 16129, 17161, 18769, 19321, 22201, 22801, 24649, 26569, 27889, 29929, 32041, 32761, 36481, 37249, 38809, 39601]
for i from 0 <= i < 46:
    primessq[i] = primessq_py[i]


def int_has_small_square_divisor(sage.rings.integer.Integer d):
    r"""
    Return the largest `a` such that `a^2` divides `d` and `a` has prime divisors `< 200`.

    EXAMPLES::

        sage: from sage.rings.number_field.totallyreal_data import int_has_small_square_divisor
        sage: int_has_small_square_divisor(500)
        100
        sage: is_prime(691)
        True
        sage: int_has_small_square_divisor(691)
        1
        sage: int_has_small_square_divisor(691^2)
        1
    """

    cdef int i
    cdef Integer asq

    asq = ZZ(1)
    for i from 0 <= i < __len_primes:
        while mpz_divisible_ui_p(d.value, primessq[i]):
            asq *= primessq[i]
            mpz_divexact_ui(d.value, d.value, primessq[i])

    return asq


cdef int eval_seq_as_poly_int(int *f, int n, int x) noexcept:
    r"""
    Evaluates the sequence a, thought of as a polynomial with

    .. MATH::

        f[n]*x^n + f[n-1]*x^(n-1) + ... + f[0].
    """
    cdef int s

    s = f[n]
    for i from n > i >= 0:
        s = s * x + f[i]
    return s

cdef double eps_abs, phi, sqrt2
eps_abs = 10.**(-12)
phi = 0.618033988749895
sqrt2 = 1.41421356237310

cdef int easy_is_irreducible(int *a, int n) noexcept:
    r"""
    Very often, polynomials have roots in {+/-1, +/-2, +/-phi, sqrt2}, so we rule
    these out quickly.  Returns 0 if reducible, 1 if inconclusive.
    """
    cdef int s, t, st, i

    # Check if a has a root in {1,-1,2,-2}.
    if eval_seq_as_poly_int(a,n,1) == 0 or eval_seq_as_poly_int(a,n,-1) == 0 or eval_seq_as_poly_int(a,n,2) == 0 or eval_seq_as_poly_int(a,n,-2) == 0:
        return 0

    # Check if f has factors x^2-x-1, x^2+x-1, x^2-2, respectively.
    # Note we only call the ZZx constructor if we're almost certain to reject.
    if fabs(eval_seq_as_poly(a,n,-phi)) < eps_abs:
        s = 2*a[n]
        t = 0
        for i from n > i >= 0:
            st = (s+t)//2
            s = 2*t+st+2*a[i]
            t = st
        if s == 0 and t == 0:
            return 0
    if fabs(eval_seq_as_poly(a,n,phi)) < eps_abs:
        s = 2*a[n]
        t = 0
        for i from n > i >= 0:
            st = (s-t)//2
            s = 2*t-st+2*a[i]
            t = st
        if s == 0 and t == 0:
            return 0
    if fabs(eval_seq_as_poly(a,n,sqrt2)) < eps_abs:
        s = a[n]
        t = 0
        for i from n > i >= 0:
            st = s
            s = 2*t+a[i]
            t = st
        if s == 0 and t == 0:
            return 0

    return 1


def easy_is_irreducible_py(f):
    """
    Used solely for testing easy_is_irreducible.

    EXAMPLES::

      sage: sage.rings.number_field.totallyreal_data.easy_is_irreducible_py(pari('x^2+1'))
      1
      sage: sage.rings.number_field.totallyreal_data.easy_is_irreducible_py(pari('x^2-1'))
      0
    """
    cdef int a[10]

    for i from 0 <= i < len(f):
        a[i] = f[i]
    return easy_is_irreducible(a, len(f)-1)


#****************************************************************************
# Main class and routine
#****************************************************************************

# Global precision to find roots; this should probably depend on the
# architecture in some way.  Algorithm gives provably correct results
# for any eps, but an optimal value of eps will be neither too large
# (which gives trivial bounds on coefficients) nor too small (which
# spends needless time computing higher precision on the roots).
cdef double eps_global
eps_global = 10.**(-4)

from sage.rings.number_field.totallyreal_phc import __lagrange_bounds_phc

cdef class tr_data:
    r"""
    This class encodes the data used in the enumeration of totally real
    fields.

    We do not give a complete description here.  For more information,
    see the attached functions; all of these are used internally by the
    functions in :mod:`.totallyreal`, so see that file for examples and
    further documentation.
    """

    def __init__(self, int n, B, a=[]):
        r"""
        Initialization routine (constructor).

        INPUT:

        - ``n`` -- integer; the degree
        - ``B`` -- integer; the discriminant bound
        - ``a`` -- list (default: ``[]``); the coefficient list to begin with,
          where ``a[len(a)]*x^n + ... + a[0]x^(n-len(a))``

        OUTPUT:

        the data initialized to begin enumeration of totally real fields
        with degree n, discriminant bounded by B, and starting with
        coefficients a.

        EXAMPLES::

            sage: T = sage.rings.number_field.totallyreal_data.tr_data(2,100)
            sage: T.printa()
            k = 0
            a = [0, -1, 1]
            amax = [0, 0, 1]
            beta =  [...]
            gnk =  [...]
        """

        cdef int i

        # Initialize constants.
        self.n = n
        self.B = B
        self.gamma = hermite_constant(n-1)

        # Declare the coefficients of the polynomials (and max such).
        self.a = <int*>sig_malloc(sizeof(int)*(n+1))
        if self.a == NULL: raise MemoryError
        self.amax = <int*>sig_malloc(sizeof(int)*(n+1))
        if self.amax == NULL: raise MemoryError
        # df is memory set aside for the derivative, as
        # used in Newton iteration above.
        self.df = <int*>sig_malloc(sizeof(int)*(n+1))
        if self.df == NULL: raise MemoryError

        for i from 0 <= i < n+1:
            self.a[i] = 0
            self.amax[i] = 0
            self.df[i] = 0

        # beta is an array of arrays (of length n) which list the
        # roots of the derivatives.
        self.beta = <double*>sig_malloc(sizeof(double)*n*(n+1))
        if self.beta == NULL: raise MemoryError
        # gnk is the collection of (normalized) derivatives.
        self.gnk = <int*>sig_malloc(sizeof(int)*(n+1)*n)
        if self.gnk == NULL: raise MemoryError

        for i from 0 <= i < (n+1)*n:
            self.beta[i] = <double>0
            self.gnk[i] = 0

        # Initialize variables.
        if a == []:
            # No starting input, all polynomials will be found; initialize to zero.
            a = [0]*n + [1]
            for i from 0 <= i < n+1:
                self.a[i] = a[i]
                self.amax[i] = a[i]
            self.a[n-1] = -(n//2)
            self.amax[n-1] = 0
            self.k = n-2
        elif len(a) <= n+1:
            # First few coefficients have been specified.
            # The value of k is the largest index of the coefficients of a which is
            # currently unknown; e.g., if k == -1, then we can iterate
            # over polynomials, and if k == n-1, then we have finished iterating.
            if a[len(a)-1] != 1:
                raise ValueError("a[len(a)-1](=%s) must be 1 so polynomial is monic" % a[len(a)-1])

            k = n-len(a)
            self.k = k
            a = [0]*(k+1) + a
            for i from 0 <= i < n+1:
                self.a[i] = a[i]
                self.amax[i] = a[i]

            # Bounds come from an application of Lagrange multipliers in degrees 2,3.
            self.b_lower = -1./n*(self.a[n-1] + (n-1.)*sqrt((1.*self.a[n-1])**2 - 2.*(1+1./(n-1))*self.a[n-2]))
            self.b_upper = -1./n*(self.a[n-1] - (n-1.)*sqrt((1.*self.a[n-1])**2 - 2.*(1+1./(n-1))*self.a[n-2]))
            if k < n-3:
                bminmax = lagrange_degree_3(n,a[n-1],a[n-2],a[n-3])
                if bminmax:
                    self.b_lower = bminmax[0]
                    self.b_upper = bminmax[1]

            # Annoying, but must reverse coefficients for NumPy.
            gnk = [int(binomial(j,k+2))*a[j] for j in range(k+2,n+1)]
            gnk.reverse()
            import numpy
            rts = numpy.roots(gnk).tolist()
            rts.sort()
            self.beta[(k+1)*(n+1)+0] = self.b_lower
            for i from 0 <= i < n-k-2:
                self.beta[(k+1)*(n+1)+(i+1)] = rts[i]
            self.beta[(k+1)*(n+1)+(n-k-1)] = self.b_upper

            # Now to really initialize gnk.
            gnk = [0] + [binomial(j,k+1)*a[j] for j in range(k + 2, n + 1)]
            for i from 0 <= i < n-k:
                self.gnk[(k+1)*n+i] = gnk[i]
        else:
            # Bad input!
            raise ValueError("a has length %s > n+1" % len(a))

    def __dealloc__(self):
        r"""
        Destructor.
        """
        sig_free(self.df)
        sig_free(self.a)
        sig_free(self.amax)
        sig_free(self.beta)
        sig_free(self.gnk)

    def increment(self, verbose=False, haltk=0, phc=False):
        r"""
        'Increment' the totally real data to the next
        value which satisfies the bounds essentially given by Rolle's
        theorem, and return the next polynomial as a sequence of
        integers.

        The default or usual case just increments the constant
        coefficient; then inductively, if this is outside of the
        bounds we increment the next higher coefficient, and so on.

        If there are no more coefficients to be had, returns the zero
        polynomial.

        INPUT:

        - ``verbose`` -- boolean to print verbosely computational details
        - ``haltk`` -- integer; the level at which to halt the inductive
          coefficient bounds
        - ``phc`` -- boolean, if PHCPACK is available, use it when `k = n-5` to
          compute an improved Lagrange multiplier bound

        OUTPUT: the next polynomial, as a sequence of integers

        EXAMPLES::

            sage: T = sage.rings.number_field.totallyreal_data.tr_data(2,100)
            sage: T.increment()
            [-24, -1, 1]
            sage: for i in range(19): _ = T.increment()
            sage: T.increment()
            [-3, -1, 1]
            sage: T.increment()
            [-25, 0, 1]
        """
        cdef int *f_out
        cdef int i

        f_out = <int *>sig_malloc(sizeof(int) * (self.n + 1))
        if f_out == NULL:
            raise MemoryError("unable to allocate coefficient list")
        for i from 0 <= i < self.n:
            f_out[i] = 0
        f_out[self.n] = 1

        self.incr(f_out, verbose, haltk, phc)

        g = [0] * (1 + self.n)
        for i from 0 <= i <= self.n:
            g[i] = f_out[i]
        sig_free(f_out)

        return g

    cdef void incr(self, int *f_out, int verbose, int haltk, int phc) noexcept:
        r"""
        This function 'increments' the totally real data to the next
        value which satisfies the bounds essentially given by Rolle's
        theorem, and returns the next polynomial in the sequence
        f_out.

        The default or usual case just increments the constant
        coefficient; then inductively, if this is outside of the
        bounds we increment the next higher coefficient, and so on.

        If there are no more coefficients to be had, returns the zero
        polynomial.

        INPUT:

        - ``f_out`` -- integer sequence, to be written with the coefficients of
          the next polynomial
        - ``verbose`` -- boolean to print verbosely computational details
        - ``haltk`` -- integer; the level at which to halt the inductive
          coefficient bounds
        - ``phc`` -- boolean; if PHCPACK is available, use it when ``k == n-5``
          to compute an improved Lagrange multiplier bound

        OUTPUT:

            None. The return value is stored in the variable f_out.
        """
        cdef int n, np1, k, i, nk, kz
        cdef int *gnkm
        cdef int *gnkm1
        cdef double *betak
        cdef double bl, br, akmin, akmax, tmp_dbl
        cdef bint maxoutflag

        n = self.n
        np1 = n+1
        k = self.k

        # If k == -1, we have a full polynomial, so we add 1 to the constant coefficient.
        if k == -1:
            self.a[0] += 1
            # Can't have constant coefficient zero!
            if self.a[0] == 0:
                self.a[0] += 1
            if self.a[0] <= self.amax[0] and easy_is_irreducible(self.a, n):
                for i from 0 <= i < n:
                    f_out[i] = self.a[i]
                return
            else:
                if verbose:
                    print(" ", end="")
                    for i from 0 <= i < np1:
                        print(self.a[i], end="")
                    print(">", end="")
                    for i from 0 <= i < np1:
                        print(self.amax[i], end="")
                    print("")

                # Already reached maximum, so "carry the 1" to find the next value of k.
                k += 1
                while k <= n-1 and self.a[k] >= self.amax[k]:
                    k += 1
                self.a[k] += 1
                self.gnk[n*k] = 0
                k -= 1
        # If we are working through an initialization routine, treat that.
        elif haltk and k == haltk-1:
            self.a[k] += 1
            if self.a[k] > self.amax[k]:
                k += 1
                while k <= n-1 and self.a[k] >= self.amax[k]:
                    k += 1
                self.a[k] += 1
                self.gnk[n*k] = 0
                k -= 1

        # If in the previous step we finished all possible values of
        # the lastmost coefficient, so we must compute bounds on the next coefficient.
        # Recall k == n-1 implies iteration is complete.
        while k < n-1:
            # maxoutflag flags a required abort along the way
            maxoutflag = 0

            # Recall k == -1 means all coefficients are good to go.
            while k >= 0 and (not haltk or k >= haltk):
                if verbose:
                    print(k, ":", end="")
                    for i from 0 <= i < np1:
                        print(self.a[i], end="")
                    print("")

                if k == n-2:
                    # We only know the value of a[n-1], the trace.  Need to apply
                    # basic bounds from Hunter's theorem and Siegel's theorem, with
                    # improvements due to Smyth to get bounds on a[n-2].
                    bl = 1./2*(1-1./n)*(1.*self.a[n-1])**2 \
                         - 1./2*self.gamma*(1./n*self.B)**(1./(n-1))
                    self.a[k] = lrint(ceil(bl))
                    br = 1./2*(1.*self.a[n-1])**2 - 0.88595*n
                    self.amax[k] = lrint(floor(br))

                    # If maximum is already greater than the minimum, break!
                    if self.a[k] > self.amax[k]:
                        if verbose:
                            print(" ", end="")
                            for i from 0 <= i < np1:
                                print(self.a[i], end="")
                            print(">", end="")
                            for i from 0 <= i < np1:
                                print(self.amax[i], end="")
                            print("")
                        maxoutflag = 1
                        break

                    # Knowing a[n-1] and a[n-2] means we can apply bounds from
                    # the Lagrange multiplier in degree 2, which can be solved
                    # immediately.
                    self.b_lower = -1./n*(self.a[n-1] + (n-1.)*sqrt((1.*self.a[n-1])**2 - 2.*(1+1./(n-1))*self.a[n-2]))
                    self.b_upper = -1./n*(self.a[n-1] - (n-1.)*sqrt((1.*self.a[n-1])**2 - 2.*(1+1./(n-1))*self.a[n-2]))

                    # Initialize the second derivative.
                    self.beta[k*np1+0] = self.b_lower
                    self.beta[k*np1+1] = -self.a[n-1]*1./n
                    self.beta[k*np1+2] = self.b_upper
                    self.gnk[k*n+0] = 0
                    self.gnk[k*n+1] = (n-1)*self.a[n-1]
                    self.gnk[k*n+2] = n*(n-1)/2

                    if verbose:
                        print(" ", '%.2f' % self.beta[k * np1 + 1])
                else:
                    # Compute the roots of the derivative.
                    self.gnk[(k+1)*n+0] += self.a[k+1]
                    newton_in_intervals(&self.gnk[(k+1)*n], self.df, n-(k+1),
                                        &self.beta[(k+1)*np1],
                                        eps_global, &self.beta[k*np1+1])
                    if verbose:
                        print(" ", end="")
                        for i from 0 <= i < n-k-1:
                            print('%.2f' % self.beta[k * np1 + 1 + i], end="")
                        print("")

                    for i from 0 <= i < n-k-1:
                        if fabs(self.beta[k*np1+i]
                                 - self.beta[k*np1+(i+1)]) < 10*eps_global:
                            # This happens reasonably infrequently, so calling
                            # the Python routine should be sufficiently fast...
                            f = ZZx([self.gnk[(k+1)*n+i] for i in range(n-(k+1)+1)])
                            # Could just take self.gnk(k+2)*n+i, but this may not be initialized...
                            df = ZZx([i*self.gnk[(k+1)*n+i] for i in range(1,n-(k+1)+1)])
                            if gcd(f,df) != 1:
                                if verbose:
                                    print("  gnk has multiple factor!")
                                maxoutflag = 1
                                break
                    if maxoutflag:
                        break

                    # Bounds come from an application of Lagrange multipliers in degrees 2,3.
                    if k == n-3:
                        self.b_lower = -1./n*(self.a[n-1] + (n-1.)*sqrt((1.*self.a[n-1])**2 - 2.*(1+1./(n-1))*self.a[n-2]))
                        self.b_upper = -1./n*(self.a[n-1] - (n-1.)*sqrt((1.*self.a[n-1])**2 - 2.*(1+1./(n-1))*self.a[n-2]))
                    elif k == n-4:
                        # New bounds from Lagrange multiplier in degree 3.
                        bminmax = lagrange_degree_3(n,self.a[n-1],self.a[n-2],self.a[n-3])
                        if bminmax:
                            self.b_lower = bminmax[0]
                            self.b_upper = bminmax[1]
                    elif k == n-5 and phc:
                        # New bounds using phc/Lagrange multiplier in degree 4.
                        bminmax = __lagrange_bounds_phc(n, 4, [self.a[i] for i from 0 <= i <= n])
                        if bminmax:
                            self.b_lower = bminmax[0]
                            self.b_upper = bminmax[1]
                        else:
                            maxoutflag = 1
                            break

                    if verbose:
                        print("  [LM bounds:", '%.2f' % self.b_lower,
                              '%.2f' % self.b_upper, end="")
                        tb = sqrt((1.*self.a[n-1])**2 - 2.*self.a[n-2])
                        print("vs. +/-", '%.2f' % tb, ']')

                    self.beta[k*np1+0] = self.b_lower
                    self.beta[k*np1+n-k] = self.b_upper

                    # Compute next g_(n-(k+1)), k times the formal integral of g_(n-k).
                    gnkm = &self.gnk[k*n]
                    gnkm1 = &self.gnk[(k+1)*n]
                    gnkm[0] = 0
                    for i from 1 <= i < n-k+1:
                        gnkm[i] = gnkm[n+i-1]*(k+1)/i
                    nk = n-(k+1)

                    # Compute upper and lower bounds which guarantee one retains
                    # a polynomial with all real roots.
                    betak = &self.beta[k*np1]
                    akmin = -eval_seq_as_poly(gnkm, n-k, betak[nk+1]) \
                            -fabs(eval_seq_as_poly(gnkm1, nk, betak[nk+1]))*eps_global
                    for i from 1 <= i < (nk+1)/2+1:
                        # Use the fact that f(z) <= f(x)+|f'(x)|eps if |x-z| < eps
                        # for sufficiently small eps, f(z) = 0, and f''(z) < 0.
                        tmp_dbl = -eval_seq_as_poly(gnkm, n-k, betak[nk+1-2*i]) \
                                  -fabs(eval_seq_as_poly(gnkm1, nk, betak[nk+1-2*i]))*eps_global
                        if tmp_dbl > akmin:
                            akmin = tmp_dbl

                    akmax = -eval_seq_as_poly(gnkm, n-k, betak[nk]) \
                            +fabs(eval_seq_as_poly(gnkm1, n-(k+1), betak[nk]))*eps_global
                    for i from 1 <= i < nk/2+1:
                        # Similar calculus statement here.
                        tmp_dbl = -eval_seq_as_poly(gnkm, n-k, betak[nk-2*i]) \
                                  +fabs(eval_seq_as_poly(gnkm1, nk, betak[nk-2*i]))*eps_global
                        if tmp_dbl < akmax:
                            akmax = tmp_dbl

                    self.a[k] = lrint(ceil(akmin))
                    self.amax[k] = lrint(floor(akmax))

                    if self.a[n-1] == 0 and (n-k)%2 == 1:
                        # Can replace alpha by -alpha, so if all
                        # "odd" coefficients are zero, may assume next
                        # "odd" coefficient is positive.
                        kz = n-3
                        while kz > k and self.a[kz] == 0:
                            kz -= 2
                        if kz == k:
                            if self.a[k] < 0:
                                self.a[k] = 0
                    if self.a[k] == 0 and self.a[k+1] == 0:
                        self.a[k] += 1
                    # Can't have constant coefficient zero!
                    if k == 0 and self.a[k] == 0:
                        self.a[k] = 1

                    if self.a[k] > self.amax[k]:
                        if verbose:
                            print(" ", end="")
                            for i from 0 <= i < np1:
                                print(self.a[i], end="")
                            print(">", end="")
                            for i from 0 <= i < np1:
                                print(self.amax[i], end="")
                            print("")
                        maxoutflag = 1
                        break

                self.k -= 1
                k -= 1

            if not maxoutflag and easy_is_irreducible(self.a, n):
                self.k = k
                for i from 0 <= i < n:
                    f_out[i] = self.a[i]
                return
            else:
                k += 1
                while k <= n-1 and self.a[k] >= self.amax[k]:
                    k += 1
                self.a[k] += 1
                self.gnk[n*k] = 0
                k -= 1

        # k == n-1, so iteration is complete; return the zero polynomial (of degree n+1).
        self.k = k
        f_out[n] = 0
        return

    def printa(self):
        """
        Print relevant data for ``self``.

        EXAMPLES::

            sage: T = sage.rings.number_field.totallyreal_data.tr_data(3,2^10)
            sage: T.printa()
            k = 1
            a = [0, 0, -1, 1]
            amax = [0, 0, 0, 1]
            beta =  [...]
            gnk =  [...]
        """
        print("k =", self.k)
        print("a =", [self.a[i] for i in range(self.n + 1)])
        print("amax =", [self.amax[i] for i in range(self.n + 1)])
        print("beta = ", [self.beta[i] for i in range(self.n * (self.n + 1))])
        print("gnk = ", [self.gnk[i] for i in range(self.n * (self.n + 1))])
