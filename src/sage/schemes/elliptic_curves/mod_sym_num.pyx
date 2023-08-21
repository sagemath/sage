# -*- coding: utf-8 -*-
#cdivision=False
#cython: cdivision_warnings=False
#cython: profile=False
r"""
Modular symbols by numerical integration

We describe here the method for computing modular symbols
by numerical approximations of the integral of the modular form on a
path between cusps.

More precisely, let `E` be an elliptic curve and `f` the newform
associated to the isogeny class of `E`. If

.. MATH::

    \lambda(r\to r') = 2\pi i \int_{r}^{r'} f(\tau)  d\tau

then the modular symbol `[r]^{+}` is defined as the quotient of
the real part of `\lambda(\infty\to r)` by the least positive real
period of `E`. Similarly for the negative modular symbol, it is the
quotient of the imaginary part of the above by the smallest positive
imaginary part of a period on the imaginary axis.

The theorem of Manin-Drinfeld shows that the modular symbols are
rational numbers with small denominator. They are used for the
computation of special values of the L-function of `E` twisted by
Dirichlet characters and for the computation of `p`-adic L-functions.

ALGORITHM:

The implementation of modular symbols in eclib and directly in sage
uses the algorithm described in Cremona's book [Cre1997]_ and Stein's
book [St2007]_. First the space of all
modular symbols of the given level is computed, then the space
corresponding to the given newform is determined. Even if these initial
steps may take a while, the evaluation afterwards is instantaneous. All
computations are done with rational numbers and hence are exact.

Instead the method used here (see [Wu2018]_ for details)
is by evaluating the above integrals
`\lambda(r\to r')` by numerical approximation. Since we know precise
bounds on the denominator, we can make rigorous estimates on the
error to guarantee that the result is proven to be the correct rational
number.

The paths over which we integrate are split up and Atkin-Lehner
operators are used to compute the badly converging part of the integrals
by using the Fourier expansion at other cusps than `\infty`.

.. NOTE::

    There is one assumption for the correctness of these computations: The
    Manin constant for the `X_0`-optimal curve should be `1` if the curve
    lies outside the Cremona tables. This is known for all curves in the
    Cremona table, but only conjectured for general curves.

EXAMPLES:

The most likely usage for the code is through the functions
``modular_symbol`` with implementation set to "num" and through
``modular_symbol_numerical``::

    sage: E = EllipticCurve("5077a1")
    sage: M = E.modular_symbol(implementation="num")
    sage: M(0)
    0
    sage: M(1/123)
    4
    sage: Mn = E.modular_symbol_numerical(sign=-1, prec=30)
    sage: Mn(3/123)       # abs tol 1e-11
    3.00000000000018

In more details. A numerical modular symbols ``M`` is created from an
elliptic curve with a chosen ``sign`` (though the other sign will also be
accessible, too)::

    sage: E = EllipticCurve([101,103])
    sage: E.conductor()
    35261176
    sage: M = E.modular_symbol(implementation="num", sign=-1)
    sage: M
    Numerical modular symbol attached to
     Elliptic Curve defined by y^2 = x^3 + 101*x + 103 over Rational Field

We can then compute the value `[13/17]^{-}` and `[1/17]^{+}` by calling
the function ``M``. The value of `[0]^{+}=0` tells us that the rank of
this curve is positive::

    sage: M(13/17)
    -1/2
    sage: M(1/17,sign=+1)
    -3
    sage: M(0, sign=+1)
    0

One can compute the numerical approximation to these rational numbers
to any proven binary precision::

    sage: M.approximative_value(13/17, prec=2) #abs tol 1e-4
    -0.500003172770455
    sage: M.approximative_value(13/17, prec=4) #abs tol 1e-6
    -0.500000296037388
    sage: M.approximative_value(0, sign=+1, prec=6) #abs tol 1e-8
    0.000000000000000


There are a few other things that one can do with ``M``. The Manin
symbol `M(c:d)` for a point `(c:d)` in the projective line can be
computed.::

    sage: M.manin_symbol(1,5)
    -1

In some cases useful, there is a function that returns all `[a/m]^{+}`
for a fixed denominator `m`. This is rather quicker for small `m` than
computing them individually::

    sage: M.all_values_for_one_denominator(7)
    {1/7: 0, 2/7: 3/2, 3/7: 3/2, 4/7: -3/2, 5/7: -3/2, 6/7: 0}

Finally a word of warning. The algorithm is fast when the cusps involved
are unitary. If the curve is semistable, all cusps are unitary. But
rational numbers with a prime `p` dividing the denominator once, but the
conductor more than once, are very difficult. For instance for the above
example, a seemingly harmless command like ``M(1/2)`` would take a very
very long time to return a value. However it is possible to compute them
for smaller conductors::

    sage: E = EllipticCurve("664a1")
    sage: M = E.modular_symbol(implementation="num")
    sage: M(1/2)
    0

The problem with non-unitary cusps is dealt with rather easily when one
can twist to a semistable curve, like in this example::

    sage: C = EllipticCurve("11a1")
    sage: E = C.quadratic_twist(101)
    sage: M = E.modular_symbol(implementation="num")
    sage: M(1/101)
    41

AUTHORS:

- Chris Wuthrich (2013-16)

"""

# ***************************************************************************
#       Copyright (C) 2016 Chris Wuthrich <christian.wuthrich@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from cysignals.memory cimport sig_malloc, sig_free, sig_realloc
from cysignals.signals cimport sig_check

from sage.arith.misc import euler_phi, kronecker as kronecker_symbol
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.misc.verbose import verbose
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.rings.complex_mpfr import ComplexField
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_mpfr cimport RealNumber, RealField


cdef extern from "<math.h>":
    double log(double)
    double exp(double)
    double cos(double)
    double sin(double)
    double ceil(double)
    double sqrt(double)

# doing this before cimport ComplexNumber does not compile
# don't ask me why
cdef extern from "<complex.h>":
    complex cexp(complex)
    complex csqrt(complex)

ctypedef long long llong

cimport sage.rings.fast_arith
import sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_llong arith_llong
fa = sage.rings.fast_arith.arith_llong()

cdef llong llgcd(llong a, llong b) except -1:
    return fa.gcd_longlong(a,b)

cdef llong llinvmod(llong a, llong m):
    return fa.inverse_mod_longlong(a, m)

DEF TWOPI = 6.28318530717958647

# use_partials arguments take three values
# 0 don't use kappa
# 1 use kappa
# 2 decide : use if m^4 < N or m < PARTIALLIMIT
DEF PARTIAL_LIMIT = 82

# ==========================================
# the following are copied from fast_arith.pyx because
# I did not manage to import them properly

cdef llong llabs(llong x) except -1:
    r"""
    Return the absolute value of a long long.
    """
    if x < 0:
        return -x
    return x

cdef llong llsign(llong n) except -2:
    r"""
    Return the sign of a long long
    """
    if n < 0:
        return -1
    return 1

cdef llong llxgcd(llong a, llong b, llong *ss, llong *tt) except -1:
    r"""
    Compute the greatest common divisor `g` of `a` and `b`,
    which is returned as the value. The integers `s` and `t` such
    that `g=s a + t b` are returned via pointers.
    """
    cdef llong psign, qsign, p, q, r, s, c, quot, new_r, new_s
    if a == 0:
        ss[0] = 0
        tt[0] = llsign(b)
        return llabs(b)
    if b == 0:
        ss[0] = llsign(a)
        tt[0] = 0
        return llabs(a)
    psign = 1
    qsign = 1
    if a < 0:
        a = -a
        psign = -1
    if b < 0:
        b = -b
        qsign = -1
    p = 1
    q = 0
    r = 0
    s = 1
    while (b):
        c = a % b
        quot = a/b
        a = b
        b = c
        new_r = p - quot*r
        new_s = q - quot*s
        p = r
        q = s
        r = new_r
        s = new_s
    ss[0] = p * psign
    tt[0] = q * qsign
    return a


def _test_llfunctions(a, b):
    r"""
    Doctest function for the above three functions.
    Given a, b this returns the absolute value of a,
    the sign of b, their gcd g and s and t such that
    g = s*a + b*t

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import _test_llfunctions
        sage: _test_llfunctions(13,19)
        (13, 1, 1, 3, -2)
        sage: _test_llfunctions(-1234,-567)
        (1234, -1, 1, 17, -37)
        sage: _test_llfunctions(144,60)
        (144, 1, 12, -2, 5)
    """
    cdef llong s, t
    a1 = Integer(llabs(a))
    a2 = Integer(llsign(b))
    a3 = Integer(llxgcd(a, b, &s, &t))
    a4 = Integer(s)
    a5 = Integer(t)
    assert a*a4 + b*a5 == a3
    return (a1,a2,a3,a4,a5)


# ================================

# this is a llong version of a function in
# sage.modular.modsym.p1list.pyx

cdef int proj_normalise(llong N, llong u, llong  v,
                        llong* uu, llong* vv) except -1:
    r"""
    Compute the canonical representative of
    `\mathbb{P}^1(\ZZ/N\ZZ)` equivalent to `(u,v)`.

    INPUT:

    - ``N`` -- an integer (the modulus or level)

    - ``u`` -- an integer (the first coordinate of (u:v))

    - ``v`` -- an integer (the second coordinate of (u:v))

    OUTPUT: If gcd(u,v,N) = 1, then returns (in a pointer)

    - ``uu`` -- an integer

    - ``vv`` -- an integer

    if `\gcd(u,v,N) \not= 1`, returns 0, 0, 0.
    """
    cdef llong d, k, g, s, t, min_v, min_t, Ng, vNg
    #verbose("       enter proj_normalise with N=%s, u=%s, v=%s"%(N,u,v),
    #        level=5)
    if N == 1:
        uu[0] = 0
        vv[0] = 0
        return 0
    # the % operator on llongs seems not to work with negatives
    if u < 0:
        u = N - ((-u) % N)
    if v < 0:
        v = N - ((-v) % N)
    u = u % N
    v = v % N
    #verbose("       now N=%s, u=%s, v=%s"%(N,u,v), level=5)
    if u == 0:
        uu[0] = 0
        if llgcd(v, N) == 1:
            vv[0] = 1
        else:
            vv[0] = 0
        return 0
    g = llxgcd(u, N, &s, &t)
    if s < 0:
        s = N - ((-s) % N)
    if t < 0:
        t = N - ((-t) % N)
    s = s % N
    t = t % N
    if llgcd(g, v) != 1:
        uu[0] = 0
        vv[0] = 0
        return 0
    # Now g = s*u + t*N, so s is a "pseudo-inverse" of u mod N
    # Adjust s modulo N/g so it is coprime to N.
    if g != 1:
        d = N / g
        while llgcd(s, N) != 1:
            s = (s+d) % N
    #verbose("       now g=%s, s=%s, t=%s"%(g,s,t), level=5)

    # Multiply [u,v] by s; then [s*u,s*v] = [g,s*v] (mod N)
    u = g
    v = (s*v) % N
    min_v = v
    min_t = 1
    if g != 1:
        Ng = N / g
        vNg = (v * Ng) % N
        t = 1
        k = 2
        while k <= g:
            v = (v + vNg) % N
            t = (t + Ng) % N
            if v < min_v and llgcd(t, N) == 1:
                min_v = v
                min_t = t
            k += 1
    v = min_v
    uu[0] = u
    vv[0] = v
    #verbose("       leaving proj_normalise with s=%s, t=%s"%(u,v), level=5)
    return 0


def _test_proj_normalise(N, u, v):
    r"""
    The doctest function for proj_normalise.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import _test_proj_normalise
        sage: _test_proj_normalise(90,7,77)
        (1, 11)
        sage: _test_proj_normalise(90,7,78)
        (1, 24)
        sage: _test_proj_normalise(24,4,9)
        (4, 3)
        sage: _test_proj_normalise(24,4,8)
        (0, 0)
        sage: _test_proj_normalise(17,-1,-7)
        (1, 7)
     """
    cdef llong uu, vv
    _ = proj_normalise(N,u,v,&uu,&vv)
    return (Integer(uu), Integer(vv))


cdef int best_proj_point(llong u, llong v, llong N,
                         llong* uu, llong* vv) except -1:
    r"""
    Given a point `(u:v)` on the projective line modulo `N`.
    This returns a representation `(x:y)` with
    small `|x|` and `|y|`.

    In many cases this is best possible, in the sense that
    `|x|+|y|` is as small as possible, but not always.
    """
    cdef llong w, p, q, Nnew, r, a, b, si
    cdef llong x0, x1, y0, y1, t0, t1, s0, s1
    #verbose("       enter best_proj_point with N=%s, u=%s, v=%s"%(N,u,v),
    #        level=5)
    if u == 0:
        uu[0] = <llong>0
        vv[0] = <llong>1
        return 0
    if v == 0:
        uu[0] = <llong>1
        vv[0] = <llong>0
        return 0

    if llgcd(u, N) == 1:
        w = (v * llinvmod(u, N) ) % N
        y0 = <llong>0
        y1 = N
        x0 = <llong>1
        x1 = w
    elif llgcd(v, N) == 1:
        w = (u * llinvmod(v,N) ) % N
        y0 = N
        y1 = <llong>0
        x0 = w
        x1 = <llong>1
    else: # cases like (p:q) mod p*q drop here
        p = llgcd(u, N)
        q = llgcd(v, N)
        Nnew = N / p / q
        w = ( (u/p) * llinvmod(v/q, Nnew) ) % Nnew
        y0 = N/q
        y1 = <llong>0
        x0 = w*p
        x1 = q

    # y will always be the longer and x the shorter
    while llabs(x0) + llabs(x1) < llabs(y0)+llabs(y1):
        if llsign(x0) == llsign(x1):
            r = (y0+y1) / (x0+x1)
        else:
            r = (y0-y1) / (x0-x1)
        t0 = y0 - r * x0
        t1 = y1 - r * x1
        s0 =  t0 - x0
        s1 =  t1 - x1
        if llabs(s0)+llabs(s1) < llabs(t0)+llabs(t1):
            t0 = s0
            t1 = s1
        # t is now the shortest vector on the line y + RR x
        #verbose("     reduced vector to (%s,%s)"%(t0,t1), level=4)
        y0 = x0
        y1 = x1
        x0 = t0
        x1 = t1

    if llgcd(y0, y1) == 1:
        uu[0] = y0
        vv[0] = y1
        return 0
    elif llgcd(x0, x1) == 1:
        uu[0] = x0
        vv[0] = x1
        return 0
    else:
        # we fall here if the first two shortest vector are both
        # not permitted, here we do a search until we hit a solution.
        #verbose("   both shortest vectors ((%s,%s) and (%s,%s)) are not "
        #        "permitted. The result is not guaranteed to "
        #        "be best possible."%(x0, x1, y0, y1), level=3)
        r = <llong>2
        a = <llong>1
        b = a
        si = a
        t0 = a * x0 + b * y0
        t1 = a * x1 + b * y1
        while llgcd(t0, t1) != 1:
            b += si
            a -= 1
            if a == 0:
                a = <llong>(-1)
                b = r-1
                si = a
            if b == 0:
                r += 1
                a = r-1
                b = <llong>1
                si = b
            t0 =  a * x0 + b * y0
            t1 =  a * x1 + b * y1
        #verbose("works with t = %s*x+%s*y"%(a, b), level=3)
        uu[0] = t0
        vv[0] = t1
        return 0


def _test_best_proj_point(u, v, N):
    r"""
    Doctest function of best_proj_point.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import _test_best_proj_point
        sage: _test_best_proj_point(11,35,101)
        (-1, 6)
        sage: _test_best_proj_point(123456789,987654321,1000000001)
        (101, 79)
        sage: _test_best_proj_point(3,10,27)
        (3, 1)
        sage: _test_best_proj_point(-1,99,101)
        (1, 2)

    Here an example where the returned value is worse than
    the given::

        sage: _test_best_proj_point(11,1,30)
        (-13, 7)
    """
    cdef llong uu, vv
    a = best_proj_point(u, v, N, &uu, &vv)
    assert a == 0
    return (Integer(uu), Integer(vv))


#======================================================================

cdef class _CuspsForModularSymbolNumerical:
    r"""
    Minimalistic class implementing cusps (not `\infty`).
    Here a cusp is a rational number together with a level.
    This class provides the methods atkin_lehner and
    is_unitary and attaches _width, _a, _m to it.

    It is to only to be used internally.
    """
    cdef public llong _a, _m, _width
    cdef public llong _N_level # trac 29290 renamed
    cdef public Rational _r

    def __init__(self, Rational r, llong N):
        r"""
        The rational (non-infinite) cusp r on X_0(N).

        INPUT:

        - ``r`` -- a rational number

        - ``N`` -- the level as a long long

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import _CuspsForModularSymbolNumerical
            sage: r = _CuspsForModularSymbolNumerical(3/7,99)
        """
        cdef llong a, m, B
        #verbose("       enter __init__ of cusps with r=%s and N=%s"%(r,N),
        #        level=5)
        a = <llong>( r.numerator() )
        m = <llong>( r.denominator() )
        a = a % m
        if 2*a > m:
            a -= m
        self._r = Rational( (a, m) )
        B = llgcd(m, N)
        self._width = N / B
        self._a = a
        self._m = m
        self._N_level = N
        # we could make it inherit from general cusps
        # but there is no need for this here
        # from sage.modular.cusps import Cusp
        # Cusp.__init__(self, a,m)

    cdef public int is_unitary(self):
        r"""
        Return whether the cusp is unitary,
        i.e. whether there exists an Atkin-
        Lehner operator that brings it to i`\infty`.
        """
        cdef llong B
        B = llgcd(self._m, self._N_level)
        return llgcd(self._width, B) == 1

    cdef public int atkin_lehner(self, llong* res) except -1:
        r"""
        If the cusp is unitary, this returns
        an Atkin-Lehner matrix for it. It is
        returned into the pointer as a list of
        four [a,b,c,d] corresponds to
        [[a,b],[c,d]].
        """
        cdef llong Q, B, c, x, y

        #verbose("       enter atkin_lehner for cusp r=%s"%self._r, level=5)
        Q = self._width
        B = llgcd(self._m, self._N_level)
        c = self._m / B
        if llgcd(Q, B) != 1:
            raise ValueError("This cusp is not in the Atkin-Lehner "
                             "orbit of oo.")
        _ = llxgcd( self._a * Q, self._m, &x, &y)
        res[0] = Q * x
        res[1] = y
        res[2] = -c * self._N_level
        res[3] =  Q * self._a
        #verbose("       leaving atkin_lehner with w_Q = "
        #        "[%s, %s, %s, %s]"%(res[0], res[1], res[2], res[3]),
        #        level=5)
        return 0


def _test_cusps(r, N):
    r"""
    Doctest function for the above class.

    Given the rational r and an integer N,
    this gives back its width, whether it is
    unitary and if so an atkin-lehner
    matrix

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num import _test_cusps
        sage: _test_cusps(0/1,11)
        (11, 1, [[0, 1], [-11, 0]])
        sage: _test_cusps(5/7,11)
        (11, 1, [[-11, -3], [-77, -22]])
        sage: _test_cusps(3/11,11)
        (1, 1, [[4, -1], [-11, 3]])
        sage: _test_cusps(1/3,27)
        (9, 0, [[0, 0], [0, 0]])
        sage: _test_cusps(1/9,27)
        (3, 0, [[0, 0], [0, 0]])
        sage: _test_cusps(5/27,27)
        (1, 1, [[11, -2], [-27, 5]])
    """
    cdef llong *wQ = [0L, 0L, 0L, 0L]
    rc = _CuspsForModularSymbolNumerical(r, N)
    a1 = rc._width
    a2 = rc.is_unitary()
    if a2:
        _ = rc.atkin_lehner(wQ)
    a = Integer(wQ[0])
    b = Integer(wQ[1])
    c = Integer(wQ[2])
    d = Integer(wQ[3])
    assert a*d - b*c == a1*a2
    M = [[a, b], [c, d]]
    return (Integer(a1), Integer(a2), M)

# ==============================================


cdef class ModularSymbolNumerical:
    r"""
    This class assigning to an elliptic curve over `\QQ` a modular symbol.
    Unlike the other implementations this class does not precompute a
    basis for this space. Instead at each call, it evaluates the integral
    using numerical approximation. All bounds are very strictly
    implemented and the output is a correct proven rational number.

    INPUT:

    - ``E`` -- an elliptic curve over the rational numbers.

    - ``sign`` -- either -1 or +1 (default). This sets the default
      value of ``sign`` throughout the class. Both are still accessible.

    OUTPUT: a modular symbol

    EXAMPLES::

        sage: E = EllipticCurve("5077a1")
        sage: M = E.modular_symbol(implementation="num")
        sage: M(0)
        0
        sage: M(77/57)
        -1
        sage: M(33/37, -1)
        2
        sage: M = E.modular_symbol(sign=-1, implementation="num")
        sage: M(2/7)
        2

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import ModularSymbolNumerical
        sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
        sage: M(1/3, -1)
        1/2
        sage: M(1/2)
        -4/5
    """
    cdef:
        llong _N_E,  _cut_val, _t_plus, _t_minus
        llong _t_unitary_minus, _t_unitary_plus
        int _lans
        int * _ans
        double * _ans_num
        double _eps_plus, _eps_minus, _eps_unitary_plus, _eps_unitary_minus
        public RealNumber _om1, _om2
        object _E, _epsQs, _Mt, _Epari
        public dict __cached_methods
        Rational _twist_q
        Integer _D
        int _global_sign

        # debug and optimisation
        #public Integer nc_sums # number of calls to summation
        #public Integer nc_direct # number of direct integrations vs
        #public Integer nc_indirect # number of indirect
        #public Integer nc_terms # number of terms summed in total

    def __cinit__(self):
        r"""
        Initialisation function.
        Allocate memory to store the
        Fourier coefficients of the newform.

        EXAMPLES::

            sage: E = EllipticCurve([1,-1])
            sage: M = E.modular_symbol(implementation="num")
            sage: M(12/11) # indirect doctest
            1/2
        """
        self._ans_num = <double *> sig_malloc( 1002 * sizeof(double) )
        self._ans = <int*> sig_malloc(1002 * sizeof(int) )
        if self._ans is NULL or self._ans_num is NULL:
            if self._ans is not NULL: sig_free(self._ans)
            if self._ans_num is not NULL: sig_free(self._ans_num)
            raise MemoryError("Memory.")

    def __init__(self, E, sign=+1):
        r"""
        See the class docstring for full documentation.

        EXAMPLES::

            sage: E = EllipticCurve("27a1")
            sage: M = E. modular_symbol(implementation="num")
            sage: M(1/9)
            1/2
            sage: M(1/3)
            -1/6
            sage: M(1/3, -1)
            1/6
        """
        #verbose("       enter __init_ of modular symbols", level=5)
        self._E = E
        self._Epari= E.pari_mincurve()
        self._global_sign = <int>sign
        self._N_E = <llong>( E.conductor() )
        self._D = -Integer(1)
        self._set_epsQs()
        self._initialise_an_coefficients()
        self._set_den_bounds()
        self.__cached_methods = {}

        # self.nc_sums = Integer(0)
        # self.nc_direct = Integer(0)
        # self.nc_indirect = Integer(0)
        # self.nc_terms = Integer(0)

        # this is a bound to decide when to go directly to ioo
        # rather than using further convergents.
        # see symbol(r) where it is used
        self._cut_val = <llong>( E.conductor().isqrt() // 4 )
        if self._cut_val < 100:
            self._cut_val = 100
        # this is can be used to disable it
        # self._cut_val = <long>(-1)
        # verbose("       leaving __init__", level=5)

    def __dealloc__(self):
        r"""
        Free the memory of the stored Fourier coefficients
        """
        sig_free(self._ans_num)
        sig_free(self._ans)

    # == basics ================

    def __repr__(self):
        """
        String representation of modular symbols.

        EXAMPLES::

            sage: E = EllipticCurve("14a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M
            Numerical modular symbol attached to Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field
        """
        return "Numerical modular symbol attached to %s" % (self._E)

    def elliptic_curve(self):
        r"""
        Return the elliptic curve of this modular symbol.

        EXAMPLES::

            sage: E = EllipticCurve("15a4")
            sage: M = E.modular_symbol(implementation="num")
            sage: M.elliptic_curve()
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + 35*x - 28 over Rational Field
        """
        return self._E

    def __call__(self, r, int sign=0, use_twist=True):
        r"""
        The modular symbol evaluated at rational. It returns a
        rational number.

        INPUT:

        - ``r`` -- a rational (or integer)

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        - ``use_twist`` -- boolean (default: True); decides if we
          allow to use a quadratic twist.

        OUTPUT: a rational number

        EXAMPLES::

            sage: E = EllipticCurve("36a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M(2/5)
            -1/3
            sage: M(2/5, -1)
            1/2

            sage: E = EllipticCurve("54a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M(5/9)
            -1/2

            sage: E = EllipticCurve("5077a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M(234/567)
            0
            sage: M(112379/43568779)
            5
        """
        cdef Rational ra

        if sign == 0:
            sign = self._global_sign

        #verbose("       enter __call__ of modular symbols for r=%s"
        #        "and sign=%s and use_twist=%s"%(r,sign,use_twist), level=5)
        if isinstance(r, Rational):
            ra = r
        elif isinstance(r, Integer):
            ra = Rational( (0,1) )
        elif isinstance(r, sage.rings.infinity.PlusInfinity):
            return Rational(0)
        else:  # who knows
            raise ValueError("The modular symbol can be evaluated at a "
                             "rational number only.")
        if use_twist:
            if self._D == -1:
                self._set_up_twist()
            if self._D != 1:
                return self._twisted_symbol(ra, sign=sign)
        return self._evaluate(ra, sign=sign)

    def approximative_value(self, r, int sign=0, int prec=20, use_twist=True):
        r"""
        The numerical modular symbol evaluated at rational.

        It returns a real number, which should be equal
        to a rational number to the given binary
        precision ``prec``. In practice the precision is
        often much higher. See the examples below.

        INPUT:

        - ``r`` -- a rational (or integer)

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        - ``prec`` -- an integer (default 20)

        - ``use_twist`` -- True (default) allows to use a
          quadratic twist of the curve to lower the conductor.

        OUTPUT: a real number

        EXAMPLES::

            sage: E = EllipticCurve("5077a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M.approximative_value(123/567)  # abs tol 1e-11
            -4.00000000000845
            sage: M.approximative_value(123/567,prec=2) # abs tol 1e-9
            -4.00002815242902

            sage: E = EllipticCurve([11,88])
            sage: E.conductor()
            1715296
            sage: M = E.modular_symbol(implementation="num")
            sage: M.approximative_value(0,prec=2)   # abs tol 1e-11
            -0.0000176374317982166
            sage: M.approximative_value(1/7,prec=2)  # abs tol 1e-11
            0.999981178147778
            sage: M.approximative_value(1/7,prec=10) # abs tol 1e-11
            0.999999972802649
        """
        cdef Rational ra
        cdef ComplexNumber ans
        cdef double eps
        cdef object L
        cdef int cinf

        if sign == 0:
            sign = self._global_sign

        #verbose("       enter approximative_value of modular symbols for r=%s,"
        #        "sign=%s and prec=%s "%(r,sign,prec,), level=5)
        if isinstance(r, Rational):
            ra = r
        elif isinstance(r, Integer):
            ra = Rational( (0,1) )
        else:  # who knows
            raise ValueError("The modular symbol can be evaluated at a"
                             "rational number only.")
        if use_twist:
            if self._D == -1:
                self._set_up_twist()
            if self._D != 1:
                return self._twisted_approx(ra, sign=sign,
                prec=prec)

        eps = <double>2
        eps = eps ** (-prec)

        ans = self._evaluate_approx(ra, eps)

        if prec > self._om1.parent().prec():
            L = self._E.period_lattice().basis(prec=prec)
            self._om1 = L[0]
            self._om2 = L[1].imag()
            cinf = self._E.real_components()
            if cinf == 1:
                self._om2 *= 2

        if sign == 1:
            return ans.real()/ self._om1
        else:
            return ans.imag()/ self._om2

    # == initialisation ========

    def _set_epsQs(self):
        r"""
        This sets the signs of the action by the Atkin-Lehner involutions
        on the modular form

        The eigenvalues are stored in a python dict _epsQs.
        Doctest in _test_init.

        EXAMPLES::

            sage: E = EllipticCurve("20a2")
            sage: M = E.modular_symbol(implementation="num") #indirect doctest
        """
        self._epsQs = {d: prod(self._E.root_number(p)
                               for p in d.prime_divisors())
                       for d in Integer(self._N_E).divisors()}

    def _set_den_bounds(self):
        r"""
        This sets the bounds on the denominators and the allowed errors.
        There are four integers _t_plus, _t_minus, _t_unitary_plus,
        _t_unitary_minus, which are proven upper bounds for the denominator
        of modular symbols (under the assumption that the optimal
        Manin constant is a divisor of 2).
        It also sets _eps_plus etc which are the errors that we are
        allowed to do in computations.

        Doctest in _test_init

        EXAMPLES::

            sage: E = EllipticCurve("240b3")
            sage: M = E.modular_symbol(implementation="num") #indirect doctest
        """
        from sage.databases.cremona import CremonaDatabase

        cdef:
            Integer N, p, delta, t0
            Rational q_plus, q_minus, s
            int co, cinf, E0cinf
            RealNumber E0om1, E0om2, q

        #verbose("       enter _set_bounds", level=5)
        N = Integer( self._N_E )
        E = self._E
        L = E.period_lattice().basis()
        self._om1 = L[0]
        self._om2 = L[1].imag()
        cinf = E.real_components()
        if cinf == 1:
            self._om2 *= 2

        # find the best curve to compare it too.
        # if the curve is in the database,
        # we can compare to the X_0-optimal curve
        isog =  E.isogeny_class()
        if N <= CremonaDatabase().largest_conductor():
            E0 = E.optimal_curve()
        # otherwise, we take a "maximal" curve
        # that the worst that can happen and is sort of the
        # opposite of what we expect, but
        else:
            ff = lambda C: C.period_lattice().complex_area()
            E0 = min(isog.curves, key=ff)
        # E0 has now conjecturally Manin constant = 1

        # now determine the bound for E0 coming from the
        # theorem of Manin and Drinfeld.
        # delta is such that the cusps are all defined over
        # the cyclotomic field zeta_delta
        delta = Integer(1)
        for p in N.prime_divisors():
            delta *= p ** (N.valuation(p)//2)
        if delta % 2 == 0:
            delta *= 2

        # on points of good reduction, the torsion is injective
        # so t0 will be a bound for the denominators of both
        # plus and minus for E0
        p = Integer(1)
        co = 0
        t0 = Integer(0)
        while co < 5 or p < max(100,10*delta) and p < self._lans:
            p += delta
            if p.is_prime() and N % p != 0:
                t0 = t0.gcd( p + 1 - self._ans[p] )
                co += 1
            if (p-2).is_prime() and N % (p-2) != 0:
                t0 = t0.gcd( (p-1)**2 - self._ans[p-2]**2 )
                co += 1
        if E0.real_components() == 1:
            t0 *= Integer(2) # slanted lattice

        # This is a not strictly necessary precaution:
        # Cremona is not always certain to have the optimal
        # curve correctly determined. If not, the index
        # is just 2.
        # For curves outside, we could have in the very worst case
        # that the optimal curve is another maximal curve. Then
        # a factor 2 should be fine, too, but it is not guaranteed.
        t0 *= 2

        # now compare the lattice for E0 with the one for E
        L = E0.period_lattice().basis()
        E0om1 = L[0]
        E0om2 = L[1].imag()
        E0cinf = E0.real_components()
        if E0cinf == 1:
            E0om2 *= Integer(2)

        maxdeg = max(max(x) for x in isog.matrix() )
        q = self._om1 / E0om1 * maxdeg
        q_plus = q.round() / maxdeg
        q = self._om2 / E0om2 * maxdeg
        q_minus = q.round() / maxdeg
        s = q_plus * t0
        self._t_plus = s.numerator()
        s = q_minus * t0
        self._t_minus = s.numerator()

        # now to the bound for the unitary cusps
        # this is a bit better because they
        # are defined over Q
        t0 = E0.torsion_order()
        if cinf == 1:
            t0 *= Integer(2)
        s = q_plus * t0
        self._t_unitary_plus = s.numerator()
        t0 = Integer(2)
        s = q_minus * t0
        self._t_unitary_minus = s.numerator()

        if N.is_squarefree():
            self._t_plus = llgcd(self._t_plus, self._t_unitary_plus)
            self._t_minus = llgcd(self._t_minus, self._t_unitary_minus)

        # set the epsilons
        self._eps_plus = <double>self._om1 / 2 / self._t_plus
        self._eps_minus = <double>self._om2 / 2 / self._t_minus
        self._eps_unitary_plus = <double>self._om1 /2 / self._t_unitary_plus
        self._eps_unitary_minus = <double>self._om2 /2 / self._t_unitary_minus

        # this code checks if the above is ok,
        # we tested quite a few curves with it
        # change the variables to public
        # from sage.schemes.elliptic_curves.mod_sym_num import ModularSymbolNumerical
        # def dens_check(E):
        #     N = E.conductor()
        #     Cu = Gamma0(N).cusps()
        #     m = E.modular_symbol()
        #     d_plus = max( [ denominator(m(r)) for r in Cu if r != oo] )
        #     m = E.modular_symbol(-1)
        #     d_minus = max( [ denominator(m(r)) for r in Cu if r != oo] )
        #     M = ModularSymbolNumerical(E)
        #     print(E.label(), (d_plus, d_minus), (M._t_plus, M._t_minus),
        #           (M._t_unitary_plus, M._t_unitary_plus))
        #     if M._t_plus < d_plus or M._t_minus < d_minus:
        #         print("**** b u g *** ")

    def _set_up_twist(self):
        r"""
        This sets up the minimal twist. This is only called from __call__
        if use_twist is True.

        EXAMPLES::

            sage: E = EllipticCurve("63a2")
            sage: M = E.modular_symbol(implementation="num")
            sage: M(3/4, use_twist=True) # indirect doctest
            -1
        """
        cdef Integer D, ell
        cdef RealNumber qq, Db

        #verbose("       enter _set_up_twist", level=5)
        Et, D = self._E.minimal_quadratic_twist()
        # we now have to make sure that |D| is as small as
        # possible (the above may give a D != 1 and yet
        # the conductor is still the same) and we
        # have to make sure that when twisting by a
        # prime ell, the twisted curve does not have
        # additive reduction. Otherwise, unitary
        # cusps will become non-movable.
        if D != 1:
            Nt = Et.conductor()
            for ell in D.prime_divisors():
                if Nt.valuation(ell) >= 2:
                    D = D.prime_to_m_part(ell)
                    if ell % 4 == 3:
                        D = -D
            if D % 4 == 3:
                D = -D
            Et = self._E.quadratic_twist(D)
        self._D = D
        verbose("  twisting by %s to get conductor %s"%(D,Et.conductor()),
                level=2)
        # now set up period change
        if D != 1:
            self._Mt = ModularSymbolNumerical(Et)
            RR = RealField(53)
            Db = RR(D.abs())
            Db = Db.sqrt()
            # a theorem by vivek pal guarantees that qq below
            # are integers.
            if D > 0:
                qq = self._om1 * Db/ self._Mt._om1 * 2
                self._twist_q = Rational( (qq.round(), 2) )
                qq = self._om2 * Db/ self._Mt._om2 * 2
                assert self._twist_q == Rational( (qq.round(),2) )
            else:
                qq =  self._om2 * Db/self._Mt._om1  * 2
                self._twist_q =  Rational( (qq.round(), 2) )
                qq = self._om1 * Db/self._Mt._om2 * 2
                assert self._twist_q == Rational( ( qq.round(),2))

    def _round(self, RealNumber val, int sign, int unitary):
        r"""
        Round the numerical approximation to the rational.
        A warning is printed if the rounding is off by more
        than 0.1.

        INPUT:

        - ``val`` -- a real number to round

        - ``sign`` -- either +1 or -1

        - ``unitary`` -- a boolean (int)

        OUTPUT: a rational.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
            sage: va = M._from_ioo_to_r_approx(0/1,0.01)
            sage: va = va.real()
            sage: M._round(va,1,True)
            1/5
        """
        cdef:
            Rational res
            llong t, r
            RealNumber q, qt

        #verbose("       enter _round with value=%s"%val, level=5)
        if sign == 1 and unitary:
            q = val/self._om1
            t = self._t_unitary_plus
        elif sign == 1 and not unitary:
            q = val/self._om1
            t = self._t_plus
        elif sign == -1 and unitary:
            q = val/self._om2
            t = self._t_unitary_minus
        elif sign == -1 and not unitary:
            q = val/self._om2
            t = self._t_minus

        qt = q * t
        r = qt.round()
        res = Rational( (r, t) )
        err = (q-res).abs()

        if err > 0.1:
            # the following did not work (compilation failed)
            #from warnings import warn
            #warn(Rounded an error of %s, looks like a bug."%err,
            # RuntimeWarning, stacklevel=5)
            print ( "Warning: Rounded an error of ", err, ", looks like a bug "
                    + "in mod_sym_num.pyx.")
        verbose("    rounding with an error of %s"%err, level=3)
        return res

    def _initialise_an_coefficients(self):
        r"""
        Compute the Fourier coefficients `a_n` for all `n` up to 1000.
        Doctest in _test_init

        EXAMPLES::

            sage: E = EllipticCurve([-11,13])
            sage: M = E.modular_symbol(implementation="num") #indirect doctest
        """
        cdef int n
        #verbose("       enter _initialise_an_coeffients", level=5)
        n = 0
        self_ans = self._E.anlist(1000, python_ints=True)
        while n <= 1000:
            self._ans[n] = self_ans[n]
            n += 1
        self._lans = len(self_ans)
        # the 53 bit numerical version
        n = 1
        self._ans_num[0] = <double>0
        while n <= 1000:
            self._ans_num[n] = <double>(self_ans[n]) / n
            n += 1

    def _add_an_coefficients(self, int T):
        r"""
        Compute the Fourier coefficients `a_n` for all `n`
        up to and including `T`.
        Further doctests in _test_init

        INPUT:

        - ``T`` -- an integer

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
            sage: M._add_an_coefficients(10000)
        """
        cdef llong n
        #verbose("       enter add_an_coeffs with T=%s"%T, level=5)
        # we artificially add 100 extra terms, to avoid calling this
        # function again with only a few new terms
        T += 100

        self._ans_num = <double *> sig_realloc(self._ans_num,
                                                (T+2)*sizeof(double))
        self._ans = <int*> sig_realloc(self._ans, (T+2)*sizeof(int) )
        if self._ans is NULL or self._ans_num is NULL:
            if self._ans is not NULL: sig_free(self._ans)
            if self._ans_num is not NULL: sig_free(self._ans_num)
            raise MemoryError("Memory error with an coefficients.")

        verbose("   not enough precomputed coefficients, "
                "adding %s"%(T - self._lans), level=3)
        # if we add more than 20% new values, redo it from scratch
        if 5* T > 6*self._lans:
            self_ans = self._E.anlist(T+1, python_ints=True)
            n = self._lans # only copy new values
            while n <= T:
                self._ans[n] = self_ans[n]
                n += 1
            n = self._lans
            self._ans_num[0] = <double>0
            while n <= T:
                self._ans_num[n] = <double>(self_ans[n]) / <double>n
                n += 1
            # last n such that ans[n] is allowed
            self._lans = T
        # otherwise, add new values with a_k
        else:
            n = self._lans
            while n <= T:
                self._ans[n] = self._Epari.ellak(n).__int__()
                self._ans_num[n] = <double>(self._ans[n]) / <double>n
                n += 1
            # last n such that ans[n] is allowed
            self._lans = T

    def clear_cache(self):
        r"""
        Clear the cached values in all methods of this class

        EXAMPLES::

            sage: E = EllipticCurve("11a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M(0)
            1/5
            sage: M.clear_cache()
            sage: M(0)
            1/5
        """
        cadi = self.__cached_methods
        for me in cadi:
            cadi[me].clear_cache()

    #================== Low level summation =========

    def _integration_to_tau(self, ComplexNumber tau,
                            int  number_of_terms, int prec):
        r"""
        Given a point `\tau` in the upper half plane
        this returns a complex number that is a close
        approximation to the integral of the modular
        form from `i\infty` to `\tau`.

        INPUT:

        - ``tau`` -- a point in the upper half plane

        - ``number_of_terms`` -- an integer describing
          how many terms to sum

        - ``prec`` -- an integer, setting the precision
          to ``prec`` bits in all of the computation

        OUTPUT:

        - a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: I = ComplexField(53).0
            sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
            sage: M._integration_to_tau(0.01*I, 1000, 53)  # abs tol 1e-11
            0.253841860855911
            sage: M._integration_to_tau(0.01*I, 1000, 200) # abs tol 1e-20
            0.25384186085591068433775876735181198283836641798722...
            sage: M._integration_to_tau(0.001*I, 1000, 200) # abs tol 1e-20
            0.253856919662568106448824346122252246436991427234479...

            sage: E = EllipticCurve("37a1")
            sage: ms = ModularSymbolNumerical(E)
            sage: ms._integration_to_tau(0.0001*I, 1000, 53) # abs tol 1e-11
            -0.0105693920159096
            sage: ms._integration_to_tau(0.3+0.01*I,1000,60) # abs tol 1e-11
            0.41268108621256428 + 0.91370544691462463*I
        """
        #verbose("       enter _integration_to_tau with tau=%s, T=%s,"
        #        "prec=%s"%(tau,number_of_terms,prec), level=5)
        cdef ComplexNumber q, s
        cdef int n

        #self.nc_sums += 1
        #self.nc_terms += Integer(number_of_terms)

        if number_of_terms > 10000000:
            print("Warning: more than 10^7 terms to sum")
        if number_of_terms > self._lans:
            self._add_an_coefficients(number_of_terms)

        CC = ComplexField(prec)
        q = 2 * CC.pi() * CC.gens()[0]
        q *= CC(tau)
        q = q.exp()
        verbose("     start sum over %s terms "%number_of_terms, level=4)
        s = CC(0)
        n = number_of_terms
        # using Horner's rule
        while n > 0:
            sig_check()
            s *= q
            s +=  CC(self._ans[n])/n
            n -= 1
        s  *= q
        #verbose("       leaving integration_to_tau with sum=%s"%s, level=5)
        return s

    # the version using double is 70-80 times faster it seems.
    cdef complex _integration_to_tau_double(self, complex tau,
                                            int number_of_terms):
        r"""
        Given a point `\tau` in the upper half plane
        this returns a complex number that is a close
        approximation to `\lambda(tau)`,
        the integral of the modular
        form from `i\infty` to `\tau`.

        INPUT:

        - ``tau`` -- a point in the upper half plane

        - ``number_of_terms`` -- an integer describing
          how many terms to sum

        OUTPUT:

        - a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("17a1"))
            sage: M._from_ioo_to_r_approx(9/13,0.01,use_partials=0) #abs tol 1e-5
            0.386771552192424 - 2.74574021880459*I
        """
        #verbose("       enter integrations_to_tau_double with tau=%s,"
        #        " T=%s"%(tau,number_of_terms), level=5)
        cdef complex q, s
        cdef int n
        #self.nc_sums += 1
        #self.nc_terms += Integer(number_of_terms)

        if number_of_terms > 10000000:
            print("Warning: more than 10^7 terms to sum")
            #raise Warning("more than 10^7 terms to sum")
        if number_of_terms > self._lans:
            self._add_an_coefficients(number_of_terms)

        q = complex(0,TWOPI) # 2 pi i
        q *= tau
        q = cexp(q)
        verbose("     start sum over %s terms "%number_of_terms, level=4)
        s = 0
        n = number_of_terms
        # using Horner's rule
        while n > 0:
            sig_check()
            s *= q
            s += self._ans_num[n]
            n -= 1
        s *= q
        #verbose("       leaving integration_to_tau_double with sum=%s"%s,
        #        level=5)
        return s

    cdef int _partial_real_sums_double(self, double y, int m,
                                       int number_of_terms,
                                       double* res) except -1:
        r"""
        Given a real positive number `y` (representing
        the imaginary part of a point in the upper half
        plane), this returns a list of approximations
        to the partial sums `\kappa_{j,m}(y)`

        INPUT:

        - ``y`` -- a positive real number

        - ``m`` -- an integer

        - ``number_of_terms`` -- an integer describing
          how many terms to sum

        OUTPUT: a list of real numbers (in a pointer)

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("37b2"))
            sage: M._from_ioo_to_r_approx(0/1,0.01,use_partials=1) #abs tol 1e-5
            0.725681061936153
        """
        #verbose("       enter partial_real_sums_double with y=%s, m=%s,"
        #        " T=%s"%(y,m,number_of_terms), level=5)
        cdef double q, qq
        cdef int n, i
        #self.nc_sums += 1
        #self.nc_terms += Integer(number_of_terms)

        if number_of_terms > 10000000:
            print(" Warning: more than 10^7 terms to sum")

        if number_of_terms > self._lans:
            self._add_an_coefficients(number_of_terms)

        q = <double>(-TWOPI)
        q *= y
        qq = q * <double>m
        q = exp(q)
        qq = exp(qq)
        verbose("     start sum over %s terms "%number_of_terms, level=4)
        i = 0
        while i < m:
            res[i] = 0
            i += 1
        n = number_of_terms
        i = n % m
        while n > 0:
            sig_check()
            res[i] *= qq
            res[i] +=  self._ans_num[n]
            n -= 1
            i -= 1
            if i == -1:
                i = m-1
        i = 1
        while i < m:
            res[i] *= q ** i
            i += 1
        res[0] *= qq
        #verbose("       leaving _partial_real_sums_double with result %s,"
        #        " %s, ... %s"%(res[0], res[1], res[m-1]), level=5)
        return 0

    def _partial_real_sums(self, RealNumber y, int m,
                                       int number_of_terms,
                                       int prec):
        r"""
        Given a real positive number `y` (representing
        the imaginary part of a point in the upper half
        plane), this returns a list of approximations
        to the partial sums `\kappa_{j,m}(y)`

        INPUT:

        - ``y`` -- a positive real number

        - ``m`` -- an integer

        - ``number_of_terms`` -- an integer describing
          how many terms to sum

        - ``prec`` -- an integer, setting the precision
          to ``prec`` bits in all of the computation

        OUTPUT: a list of real numbers

        EXAMPLES::

            sage: from  sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve([1,0,-1064,-1,0])
            sage: m = ModularSymbolNumerical(E)
            sage: m._partial_real_sums(0.2, 7, 10000, 57)
            [-0.00008643214783335204,
             0.2846256879231407,
             -0.04049993474185631,
             -0.01536940717542450,
             -0.001640896372921393,
             4.706713577985852e-8,
             0.0001771330855478248]
        """
        #verbose("       enter partial_real_sums with y=%s, m=%s,"
        #        " T=%s"%(y,m,number_of_terms), level=5)
        cdef RealNumber q, qq
        cdef int n, i
        #self.nc_sums += 1
        #self.nc_terms += Integer(number_of_terms)

        if number_of_terms > 10000000:
            print(" Warning: more than 10^7 terms to sum")

        if number_of_terms > self._lans:
            self._add_an_coefficients(number_of_terms)

        RR = RealField(prec)
        q = - 2 * RR.pi()
        q *= RR(y)
        qq = q * m
        q = q.exp()
        qq = qq.exp()
        verbose("     start sum over %s terms "%number_of_terms, level=4)
        i = 0
        res = []
        while i < m:
            res.append(RR(0))
            i += 1
        n = number_of_terms
        i = n % m
        while n > 0:
            sig_check()
            res[i] *= qq
            res[i] += RR(self._ans[n])/n
            n -= 1
            i -= 1
            if i == -1:
                i = m-1
        i = 1
        while i < m:
            res[i] *= q ** i
            i += 1
        res[0] *= qq
        #verbose("       leaving _partial_real_sums with result"
        #        " %s, %s, ... %s"%(res[0], res[1], res[m-1]), level=5)
        return res

    #================

    def _get_truncation_and_prec(self, double y, double eps):
        r"""
        Compute the numbers of terms needed in the sum approximating
        the integral, and the precision of each term needed.

        INPUT:

        - ``y`` -- a positive real number, representing the imaginary
          part of an element in the upper half plane and

        - ``eps`` -- a positive real number, the maximal allowed error

        OUTPUT:

        two integers `T` and `b`

        If `T` would be larger than `2^31`, the value (-1,-1) is
        returned instead.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("11a2"))
            sage: M._get_truncation_and_prec(0.001,0.0001)
            (2275, 53)
            sage: M._get_truncation_and_prec(0.01,0.01)
            (118, 53)
            sage: M._get_truncation_and_prec(0.00001,0.01)
            (216405, 53)
            sage: M._get_truncation_and_prec(0.00001,0.0001)
            (283246, 58)

            sage: M = ModularSymbolNumerical(EllipticCurve("5077a1"))
            sage: M._get_truncation_and_prec(10^(-11),10^(-2))
            (-1, -1)
        """
        #verbose("       enter _get_truncation_and_prec with y=%s"
        #        " and eps=%s"%(y,eps), level=5)
        # how much of the error comes from truncation and how much
        # from precision.
        DEF split_error_truncations = 0.99
        DEF split_error_prec = 1 - split_error_truncations

        cdef double tt, bb, A, twopiy
        cdef int T, B, T0

        twopiy = <double>TWOPI
        twopiy *=  y
        tt = twopiy
        tt *= split_error_truncations * eps
        tt = log(tt)
        tt = - tt / twopiy
        if tt < 2147483648:
            T = <int>(ceil(tt))
        else:
            T = -1
            return T, T
        #verbose("      now tt =%s, twopiy=%s, T=%s"%(tt,twopiy,T), level=4)

        # the justification for these numbers is explained at the
        # very end of this file
        if T > 4324320:
            A = 6
            T0 = 4324320
        elif T > 2162160:
            A = 5
            T0 = 2162160
        elif T > 831600:
            A = 4
            T0 = 831600
        elif T > 277200:
            A = 3
            T0 = 277200
        elif T > 55440:
            A = 2
            T0 = 55440
        elif T > 10080:
            A = 3/2
            T0 = 10080
        else:  # does not improve
            A = 1
            T0 = T

        tt -= log(A)/twopiy
        T = min(T, max(<int>(ceil(tt)), T0))
        #verbose("    now tt =%s, twopiy=%s, T=%s"%(tt,twopiy,T), level=4)

        bb =  split_error_prec * eps
        bb /=  T + bb
        bb /=  2
        bb /=  T
        bb = - log(bb)/ log(2)
        B = <int>(ceil(bb))
        B = max(B, 53)
        T = max(T, 100)
        #verbose("       leaving _get_truncation_and_prec with T=%s,"
        #        " B=%s"%(T,B), level=5)
        return T, B

# ===============

    @cached_method # its modified below for max eps
    def _kappa(self, llong m, llong z, eps=None):
        r"""
        This returns all `\kappa_{j,m}(1/\sqrt{z})` for a given
        integer `m` and another integer `z`. This function `\kappa` is
        defined by the following sum:

        .. MATH::

           \kappa_{j,m}(y) = \sum_{n\equiv j \bmod m} \frac{a_n}{n}
           e^{- 2\pi n y}

        This is used to compute all the modular symbols `[a/m]` for
        a fixed `m`.

        INPUT:

        - ``m`` -- an integer (long long)

        - ``z`` -- another integer (long long)

        - ``eps`` -- either None (default) or a real number (double),
          describing the precision to which the terms are computed.
          Each term is off by less than ``eps``/``m``. If None, the
          ``eps`` is chosen from the precomputed precision related
          to the periods.

        OUTPUT: a list of `m` real numbers

        EXAMPLES::

            sage: E = EllipticCurve("43a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M._kappa(3,4) # abs tol 1e-11
            [-5.379533671373222e-05, 0.043215661934968536, -0.0018675632930897528]
            sage: M._kappa(3,17) # abs tol 1e-11
            [-0.0068222516409258815, 0.2189879706778559, -0.047856204984566546]
            sage: M._kappa(3,12345,0.01) # abs tol 1e-11
            [-0.04800196513225438, 1.501878908740486, -1.4540035671680258]
            sage: M._kappa(3,12345,0.001) # abs tol 1e-11
            [-0.04790883326924006, 1.5019073235739455, -1.4539982909123526]

        This is to check that the caching works when asked with lower
        precision::

            sage: M._kappa(7,9,0.0001) # abs tol 1e-11
            [-3.848348562241613e-46,
            0.12314471107014528,
            -0.01516461914094593,
            -0.0012449611795634324,
            0.00011498287475216501,
            -2.265525136998248e-05,
            2.3248943281270047e-06]
            sage: M._kappa(7,9,0.1) # abs tol 1e-11
            [-3.848348562241613e-46,
            0.12314471107014528,
            -0.01516461914094593,
            -0.0012449611795634324,
            0.00011498287475216501,
            -2.265525136998248e-05,
            2.3248943281270047e-06]
        """
        #verbose("       enter _kappa with m=%s, z=%s and eps=%s"%(m,z,eps),
        #        level=5)
        cdef:
            int T, prec, j
            double * ra
            object res
            double y, epsi
            RealNumber yr

        y = <double>(z)
        y = sqrt(y)
        y = 1/y
        if eps is None:
            epsi = self._eps_unitary_plus
            if epsi > self._eps_unitary_minus:
                epsi = self._eps_unitary_minus
        else:
            epsi = eps

        # if called with a previous (m,z,eps) but a larger eps,
        # return the cached value
        cac = self.__cached_methods['_kappa'].cache
        for ke in cac:
            mm, zz, eeps = ke[0]
            if mm == m and zz == z:
                if eps is None:
                    if eeps is None or eeps <= epsi:
                        return cac[ke]
                else:
                    if eeps is not None and eeps <= eps:
                        return cac[ke]

        T, prec = self._get_truncation_and_prec(y, epsi)
        if T == -1:
            raise ValueError("Too many terms > 2^31 would have to be"
                             + "summed up. Giving up.")
        T += m
        #verbose("   precision in _kappa set to %s,"
        #        " summing over %s terms"%(prec,T), level=3)

        if prec > 53:
            RR = RealField(prec)
            yr = RR(z)
            yr = yr.sqrt()
            yr = yr**(-1)
            res = self._partial_real_sums(yr, m, T, prec)
            # return a python list of doubles that we cache
            res = [<double>(res[j]) for j in range(m)]
        else:
            ra = <double *> sig_malloc( m * sizeof(double))
            if ra is NULL:
                raise MemoryError
            _ = self._partial_real_sums_double(y, m, T, ra)
            res = [ra[j] for j in range(m)]
            sig_free(ra)
        #verbose("       leaving _kappa with"
        #        " [%s, %s, ... %s]"%(res[0], res[1], res[m-1]), level=5)
        return res

    def _from_ioo_to_r_approx(self, Rational r, double eps,
                              int use_partials=2):
        r"""
        Given a cusp `r` this computes the integral `\lambda(r)`
        from `i\infty` to `r` to the given precision `eps`.

        INPUT:

        - ```r`` -- a rational number

        - ``eps`` -- a positive real number

        - ``use_partials`` -- int: 0 don't use partials,
          1 use them, 2 (default) decide if meaningful

        OUTPUT: a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve("37a1")
            sage: m = ModularSymbolNumerical(E)
            sage: m._from_ioo_to_r_approx(1/7,0.01) # abs tol 1e-11
            2.99345862520910 - 4.24742221394325e-8*I
            sage: m._from_ioo_to_r_approx(2/7,0.01)  # abs tol 1e-11
            2.45138940312063*I
            sage: m._from_ioo_to_r_approx(0/1,0.001)  # abs tol 1e-11
            -2.77555756156289e-17

            sage: E = EllipticCurve("37b1")
            sage: m = ModularSymbolNumerical(E)
            sage: m._from_ioo_to_r_approx(0/1,0.01)  # abs tol 1e-11
            0.725681061936153

            sage: E = EllipticCurve("5077a1")
            sage: m = ModularSymbolNumerical(E)
            sage: m._from_ioo_to_r_approx(-1/7,0.01, use_partials=1)  # abs tol 1e-11
            6.22747859174503 - 1.48055642530800*I
            sage: m._from_ioo_to_r_approx(-1/7,0.01, use_partials=0) #abs tol 1e-11
            6.22747410432385 - 1.48055182979493*I

        This uses 65 bits of precision::

            sage: m._from_ioo_to_r_approx(-1/7,0.0000000001) # abs tol 1e-11
            6.227531974630294568 - 1.480548268241443085*I
        """
        #verbose("       enter _from_ioo_to_r_approx with r=%s"
        #        " and eps=%s"%(r,eps), level=5)
        cdef:
            llong m, Q, epsQ, a, u
            double yy, taui
            int T, prec, j
            double complex  tauc, tauphc, int1c, int2c, twopii, ze1, ze2, su
            ComplexNumber tau, tauph, int1, int2
            llong * wQ = [0L, 0L, 0L, 0L]
            object ka

        rc = _CuspsForModularSymbolNumerical(r, self._N_E)
        Q = rc._width
        _ = rc.atkin_lehner(wQ)
        m = rc._m
        epsQ = self._epsQs[Q]
        r = rc._r

        # now we determine the best place to cut the integration.
        # we will integrate from r to r+i*yy and then to i*oo
        yy = <double>Q
        yy = sqrt(yy)
        yy = 1 / yy / m
        T, prec = self._get_truncation_and_prec(yy, eps/2)
        if T == -1:
            raise ValueError("Too many terms > 2^31 would have to be"
                             + " summed up. Giving up.")
        if m == 1:
            use_partials = 0
        if use_partials == 2:
            use_partials = (prec==53) and ( m**4 < self._N_E or m < PARTIAL_LIMIT)

        if not use_partials and prec > 53:
            CC = ComplexField(prec)
            tau = CC(-Q)
            tau = tau.sqrt()
            tau = r - 1/tau/m
            tauph = (tau * wQ[0]  + wQ[1])/(wQ[2]*tau + wQ[3])
            verbose("  computing integral from i*oo to %s using %s terms "
                    "and precision %s"%(tau, T, prec),level=2)
            int1 = self._integration_to_tau(tau, T, prec)
            verbose("  yields %s "%int1, level=2)
            verbose("  compute integral from %s to %s by computing an "
                    "integral from i*oo to %s"%(r, tau, tauph),level=2)
            int2 = self._integration_to_tau(tauph, T, prec)
            int2 *= -epsQ
            verbose("  yields %s"%int2, level=2)
            return int2 + int1

        elif not use_partials: # prec = 53
            taui = <double>(Q)
            taui = sqrt(taui)
            taui =  1/taui/m
            tauc = complex(r, taui)
            #verbose("act on %s by [[%s,%s],[%s,%s]]"%(tauc, wQ[0],
            # wQ[1], wQ[2], wQ[3]), level =4)
            tauphc = (tauc * wQ[0] + wQ[1])/(wQ[2]*tauc + wQ[3])
            verbose("  computing integral from i*oo to %s using %s terms "
                    "in fast double precision"%(tauc, T),level=2)
            int1c = self._integration_to_tau_double(tauc, T)
            verbose("  yields %s "%int1c, level=2)
            verbose("  compute integral from %s to %s by computing an "
                    "integral from i*oo to %s"%(r, tauc, tauphc),level=2)
            int2c = self._integration_to_tau_double(tauphc, T)
            int2c *= -epsQ
            verbose("  yields %s"%int2c, level=2)
            CC = ComplexField(prec)
            return CC(int2c + int1c)

        else: # use_partial
            verbose("  computing integral from i*oo to %s using "
                    "using partials with "
                    "y =%s"%(r, yy), level=2)
            ka = self._kappa(m, m*m*Q,eps/2)
            a = rc._a
            twopii = TWOPI * complex("j")
            ze1 = twopii / m * a
            ze1 = cexp(ze1)
            u = llinvmod(Q * a, m)
            ze2 = - twopii / m * u
            ze2 =  cexp(ze2)
            j = 0
            su = 0
            verbose("      summing up %s partial sums, having set u = %s,"
                    " z1 =%s, z2=%s"%(m,u,ze1,ze2), level=4)
            while j < m:
                sig_check()
                su += ka[j] * ( (ze1 ** j) - epsQ * (ze2 ** j))
                j += 1
            CC = ComplexField(prec)
            return CC(su)

    cdef _from_r_to_rr_approx_direct(self, Rational r, Rational rr,
                                     Integer epsQ, Integer epsQQ,
                                    llong* wQ, llong* wQQ,
                                    int T, int prec, double eps,
                                    int use_partials=2):
        r"""
        This is just a helper function for _from_r_to_rr_approx. In case
        the integral is evaluated directly this function is called.

        INPUT:

        - ``r``, ``rr`` -- two Rationals

        - ``espQ`` and ``espQ`` -- two Integers

        - `` T `` and ``prec`` -- two ints

        - `` use_partials -- int: 0 don't use partials,
          1 use them, 2 (default) decide if meaningful

        OUTPUT: a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
            sage: M._from_r_to_rr_approx(1/3,0/1,0.000001,use_partials=0) # abs tol 1e-11
            -0.634604652139777 + 1.45881661693850*I
            sage: M._from_r_to_rr_approx(1/3,0/1,0.000001,use_partials=1) # abs tol 1e-11
            -0.634604652139777 + 1.45881661693849*I
        """
        cdef:
            ComplexNumber tau0, tau1, int1, int2
            RealNumber x1, x2, s
            complex tau0c, tau1c, int1c, int2c, ze1, ze2, su, twopii
            llong g, u, v, uu, vv, D, a, aa, m, mm, Q, QQ, z, xi, xixi
            int oi, j
            double x1d, x2d, sd

        #verbose("       enter _from_r_to_rr_approx_direct with r=%s,"
        #        " rr=%s,..."%(r,rr), level=5)
        rc = _CuspsForModularSymbolNumerical(r, self._N_E)
        m = rc._m
        a = rc._a
        Q = rc._width
        rrc = _CuspsForModularSymbolNumerical(rr,self._N_E)
        mm = rrc._m
        aa = rrc._a
        QQ = rrc._width
        oi = llxgcd(Q*a, m, &u, &v)
        oi = llxgcd(QQ*aa, mm, &uu, &vv)
        #verbose("  The inverses are (u,v)=(%s,%s) and
        # (uu,vv)=%s,%s"%(u,v,uu,vv), level=3)

        if use_partials == 2:
            g = llgcd(Q,QQ)
            D = Q * QQ
            D /= g
            D *=  llabs(a*mm-aa*m)
            use_partials = (prec==53) and ( D**4 < self._N_E or D < PARTIAL_LIMIT)

        CC = ComplexField(prec)
        if not use_partials and prec > 53:
            RR = RealField(prec)
            x1 = RR(Q)
            x1 = x1*aa*u + mm*v
            x1 = x1/Q/(a*mm-aa*m)
            x2 = RR(QQ)
            x2 = x2*a*uu + m*vv
            x2 = x2/QQ/(aa*m-a*mm)
            s = RR(Q)
            s = s*QQ
            s = s.sqrt()
            s = s * llabs(a*mm-aa*m)
            s = 1/s
            tau0 = CC(x1, s)
            tau1 = CC(x2, s)
            #verbose("  two points are %s and %s"%(tau0, tau1), level=3)
            verbose("   computing integral from %s to tau by computing "
                    "the integral from i*oo to %s"%(r,tau0), level=3)
            int1 = self._integration_to_tau(tau0, T, prec)
            int1 *= - epsQ
            verbose("   yields %s "%int1, level=3)
            verbose("   computing integral from tau to %s by computing "
                    "the integral from i*oo to %s"%(rr, tau1), level=3)
            int2 = self._integration_to_tau(tau1, T, prec)
            int2 *= epsQQ
            verbose("   yields %s "%int2, level=3)
            ans = int2 + int1
        elif not use_partials:
            x1d = Q
            x1d = x1d*aa*u + mm*v
            x1d = x1d/Q/(a*mm-aa*m)
            x2d = QQ
            x2d = x2d*a*uu + m*vv
            x2d = x2d/QQ/(aa*m-a*mm)
            sd = Q
            sd = sd*QQ
            sd = sqrt(sd)
            sd = sd * llabs(a*mm-aa*m)
            sd = 1/sd
            tau0c = complex(x1d,sd)
            tau1c = complex(x2d,sd)
            verbose("   computing integral from %s to tau by computing "
                    "the integral from i*oo to %s"%(r,
                    tau0c),level=3)
            int1c = self._integration_to_tau_double(tau0c, T)
            int1c *=  - epsQ
            verbose("   yields %s "%int1c, level=3)
            verbose("   computing integral from tau to %s by computing "
                    "the integral from i*oo to "
                    "%s"%(rr, tau1c), level=3)
            int2c = self._integration_to_tau_double(tau1c, T)
            int2c *= epsQQ
            verbose("   yields %s "%int2c, level=3)
            ans = int2c + int1c
        else: # use_partials
            g = llgcd(Q,QQ)
            D = Q * QQ
            D /= g
            D *=  llabs(a*mm-aa*m)
            xi = (Q*aa*u+v*mm) * QQ /g * llsign(a*mm-aa*m)
            xixi = (QQ*a*uu+vv*m) * Q /g * llsign(aa*m-a*mm)
            z = Q * QQ * (a*mm-aa*m)**2
            ka = self._kappa(D, z, eps/2)
            twopii = TWOPI * complex("j")
            ze1 = twopii / D * xixi
            ze1 = cexp(ze1)
            ze2 = twopii / D * xi
            ze2 = cexp(ze2)
            j = 0
            su = 0
            verbose("     summing up %s partial sums"%D, level=4)
            while j < D:
                sig_check()
                su += (ze1**j * epsQQ  - ze2**j * epsQ ) * ka[j]
                j += 1
            ans = su
        return CC(ans)

    def _from_r_to_rr_approx(self, Rational r, Rational rr, double eps,
                             method=None, int use_partials=2):
        r"""
        Given a cusp `r` this computes the integral `\lambda(r\to r')`
        from `r` to `r'` to the given precision ``eps``.

        INPUT:

        - ``r`` -- a rational number `r`

        - ``rr`` -- a second rational number `r'`

        - ``eps`` -- a positive real number

        - ``method`` -- a string or None: either "direct", "indirect",
          "both". When method is not given (default), then the better
          of the two is chosen. "both" raises an error if the two
          methods differ by more than ``eps``.

        - ``use_partials`` -- int: 0 don't use partials,
          1 use them, 2 (default) decide if meaningful

        OUTPUT: a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
            sage: M._from_r_to_rr_approx(0/1,2/5,0.01) # abs tol 1e-11
            1.90381395641933 - 1.45881661693850*I
            sage: M._from_r_to_rr_approx(0/1,2/5,0.001, "both", use_partials=0) #abs tol 1e-11
            1.90381395641931 - 1.45881661693851*I
            sage: M._from_r_to_rr_approx(0/1,2/5,0.001, "both", use_partials=1) #abs tol 1e-11
            1.90381395641933 - 1.45881661693850*I

            sage: M._from_r_to_rr_approx(1/11,1/7,0.001) # abs tol 1e-11
            -0.888446512995687 + 1.45881661693850*I
            sage: M._from_r_to_rr_approx(0/1,44/98761,0.001) # abs tol 1e-11
            0.634604184365293 + 1.45881886531983*I
            sage: M._from_r_to_rr_approx(0/1,123/456,0.0000001) # abs tol 1e-11
            1.26920930437008 - 2.91763323391590*I

            sage: M = ModularSymbolNumerical(EllipticCurve("389a1"))
            sage: M._from_r_to_rr_approx(0/1,1/5,0.0001, "both") #abs tol 1e-5
            -4.98042489791268 - 6.06124058444291e-8*I
            sage: M._from_r_to_rr_approx(1/3,1/7,0.001, "both", use_partials=0) #abs tol 1e-5
            -2.49020330904154 + 5.91520739782657*I
            sage: M._from_r_to_rr_approx(1/3,1/7,0.0001, "both", use_partials=0) #abs tol 1e-5
            -2.49021239793944 + 5.91521298876557*I
            sage: M._from_r_to_rr_approx(1/3,1/7,0.0001, "both") #abs tol 1e-5
            -2.49021247526015 + 5.91521314939661*I


            sage: M = ModularSymbolNumerical(EllipticCurve("5077a1"))
            sage: M._from_r_to_rr_approx(0/1, -7/3179, 0.001) #abs tol 1e-5 #long time
            -6.22753195853020 - 5.92219314384901*I

            sage: E = EllipticCurve([91,127])
            sage: M = ModularSymbolNumerical(E)
            sage: M._from_r_to_rr_approx(7/11,14/75,0.01,"direct")
            Traceback (most recent call last):
            ...
            ValueError: Too many terms > 2^31 would have to be summed up. Giving up.
        """
        cdef:
            llong m, Q, mm, QQ, a, aa
            llong * wQ = [0L, 0L, 0L, 0L]
            llong * wQQ = [0L, 0L, 0L, 0L]
            Integer epsQ, epsQQ
            double s, yy
            ComplexNumber ans, ans2
            int  T=0, prec=0, T1=0, T2=0, oi

        #verbose("       enter _from_r_to_rr_approx with r=%s,"
        #        " rr=%s, "%(r,rr), level=5)
        rc = _CuspsForModularSymbolNumerical(r, self._N_E)
        m = rc._m
        a = rc._a
        Q = rc._width
        oi = rc.atkin_lehner(wQ)
        epsQ = self._epsQs[Q]
        r = rc._r

        rrc = _CuspsForModularSymbolNumerical(rr,self._N_E)
        mm = rrc._m
        aa = rrc._a
        QQ = rrc._width
        oi = rrc.atkin_lehner(wQQ)
        epsQQ = self._epsQs[QQ]
        rr = rrc._r

        if r == rr:
            return ComplexField(53)(0)

        # now we determine the best place to cut the direct integration.
        # we will integrate from r to x+i*y and then to rr
        if method != "indirect":
            # overflow danger
            s = <double>Q
            s *= QQ
            s = sqrt(s)
            s = s * llabs(a*mm-aa*m)
            s = 1/s
            #verbose("    direct method goes to s=%s and eps=%s"%(s,eps),
            #        level=3)
            T, prec = self._get_truncation_and_prec(s, eps/2)
            #verbose("    giving T=%s and prec=%s"%(T,prec), level=3)
            if T == -1:
                if method == "direct" or method == "both":
                    raise ValueError("Too many terms > 2^31 would have to"
                                     + " be summed up. Giving up.")
                else: # method was None
                    method = "indirect"

        # now we compare it to the indirect integration via i*oo
        if method != "direct":
            yy = <double>Q
            yy = sqrt(yy)
            yy = 1/yy/m
            T1, prec1 = self._get_truncation_and_prec(yy, eps/4)
            yy = <double>QQ
            yy = sqrt(yy)
            yy = 1/yy/mm
            T2, prec2 = self._get_truncation_and_prec(yy, eps/4)
            if T1 == -1 or T2 == -1:
                if method == "indirect" or method == "both":
                    raise ValueError("Too many terms > 2^31 would have to"
                                     + " be summed up. Giving up.")
                else: # method was None
                    method = "direct"

        if method is None:
            if prec > 53 and prec1 == 53 and prec2 == 53:
                method = "indirect"
            elif prec == 53 and (prec1>53 or prec2>53):
                method = "direct"
            elif T < T1 + T2:
                method = "direct"
            else:
                method = "indirect"

        if method == "direct" or method == "both":
            verbose(" using the direct integration from %s to %s with "
                    "%s terms to sum"%(r, rr, T), level=2)
            #self.nc_direct += 1
            ans = self._from_r_to_rr_approx_direct(r, rr, epsQ, epsQQ,
                                                   wQ, wQQ, T, prec, eps,
                                                   use_partials)
            if method != "both":
                return ans

        if method == "indirect" or method == "both":
            verbose("  using the indirect integration from %s to %s "
                    "with %s terms to sum"%(r, rr, T1+T2), level=2)
            #self.nc_indirect += 1
            ans2 = ( self._from_ioo_to_r_approx(r, eps/2,
                                                use_partials=use_partials)
                    - self._from_ioo_to_r_approx(rr, eps/2,
                                                use_partials=use_partials) )
            if method != "both":
                return ans2

        if method == "both":
            if not use_partials:
                assert (ans - ans2).abs() < eps, ("Bug in modular symbols. "
                     + "The indirect and direct computation of the modular "
                     + "symbol from %s to %s differ by too much"%(r, rr) )

            return (ans + ans2)/2

    def _transportable_approx(self, Rational r, Rational rr, double eps):
        r"""
        Given two cusps `r` and `r'`, which are `\Gamma_0(N)`-equivalent,
        this computes the integral `\lambda(r\to r')` by transporting
        the path to a path between `\tau` and `\gamma(\tau)` for
        `\gamma\in\Gamma_0(N)` and a well-chosen `\tau`.

        So far we have only implemented that it is transported close
        to `0` or `i\infty`. Note that this method does not require the
        cusps to be unitary.

        INPUT:

        - ``r`` -- a rational number

        - ``rr`` -- another rational number

        - ``eps`` -- a positive real number

        OUTPUT: a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("11a1"))
            sage: M._transportable_approx(0/1,-2/7,0.001) # abs tol 1e-11
            -0.634604652139777 + 1.45881661693850*I
            sage: M._from_r_to_rr_approx(0/1,-2/7,0.001) # abs tol 1e-11
            -0.634604652139776 + 1.45881661693850*I
            sage: M._from_r_to_rr_approx(0/1,-2/7,0.0001) # abs tol 1e-11
            -0.634604652139776 + 1.45881661693850*I
            sage: M._from_r_to_rr_approx(0/1,-2/7,0.00001) # abs tol 1e-11
            -0.634604652139776 + 1.45881661693850*I
            sage: M._from_r_to_rr_approx(0/1,-2/7,0.000001) # abs tol 1e-11
            -0.634604652139776 + 1.45881661693850*I

            sage: M = ModularSymbolNumerical(EllipticCurve("37a1"))
            sage: M._transportable_approx(0/1,-1/19,0.001) # abs tol 1e-11
            -1.14879551982305e-8 + 1.65994273881864e-10*I
            sage: M._transportable_approx(0/1,-4/17,0.001) # abs tol 1e-11
            -2.99345356727791 + 2.45138870627435*I
            sage: M._transportable_approx(0/1,-4/17,0.00001) # abs tol 1e-11
            -2.99345863532461 + 2.45138938269979*I
            sage: M._from_r_to_rr_approx(0/1,-4/17,0.00001) # abs tol 1e-11
            -2.99345862657997 + 2.45138938852658*I

        This goes via i `\infty`::

            sage: M = ModularSymbolNumerical(EllipticCurve("5077a1"))
            sage: M._transportable_approx( 0/1, -35/144, 0.001) #abs tol 1e-11
            -6.22753189644996 + 3.23405342839145e-7*I
            sage: M._from_r_to_rr_approx( 0/1, -35/144, 0.001) # abs tol 1e-10
            -6.22753204310913 - 1.31710951034592e-8*I

        While this one goes via 0::

            sage: M._transportable_approx( 0/1, -7/31798, 0.001) #abs tol 1e-11
            -7.01577418382726e-9 - 7.40274138232394*I
            sage: M._from_r_to_rr_approx( 0/1, -7/31798, 0.001) #abs tol 1e-5 #long time
            -7.02253033502132e-9 - 7.40274138234031*I
        """
        cdef:
            Integer c, a, b, d
            int T, prec
            double yy
            ComplexNumber tau1, tau2, int1, int2
            complex tau1c, tau2c, int1c, int2c

        #verbose("       enter transportable_symbol_approx with r=%s,"
        #        " rr=%s"%(r,rr), level=5)

        #this finds a gamma with smallest |c|
        from sage.modular.cusps import Cusp
        rc = Cusp(r)
        boo, ga = rc.is_gamma0_equiv(rr, self._N_E, "matrix")

        if not boo:
            raise ValueError("The cusps %s and %s are not "
                             "Gamma_0(%s)-equivalent"%(r, rr, self._N_E))

        # now find the same for the move to 0
        c = ga[1][0]
        r0 = - Rational( (c/self._N_E, ga[0][0]) )
        rc0 = Cusp(r0)
        _, ga0 = rc0.is_gamma0_equiv(0, self._N_E, "matrix")

        if c.abs() > ga0[1][0].abs(): # better at 0
            ga = ga0
            c = ga[1][0]
            eN = -self._epsQs[self._N_E]
        else: #better at i oo
            eN = 1

        a = ga[0][0]
        b = ga[0][1]
        d = ga[1][1]

        if (a + d).abs() <= 2:
            # parabolic matrix has a fixed point among the cusps
            # and so we can simply give back 0
            # because ga has then finite order
            CC = ComplexField(53)
            return CC(0)

        yy = <double>(c.abs())
        yy = 1/yy
        T, prec = self._get_truncation_and_prec(yy, eps/2)
        if T == -1:
            raise ValueError("Too many terms > 2^31 would have to be summed"
                             + "up. Giving up.")

        if prec > 53:
            CC = ComplexField(prec)
            tau1 = CC(-d/c, 1/c.abs())
            # computes the integral from tau to i*oo
            verbose("   computing integral from i*oo to %s using %s terms "
                    "and precision %s"%(tau1, T, prec), level=3)
            int1 = self._integration_to_tau(tau1, T, prec)
            verbose("   yields %s "%int1, level=3)
            tau2 = (tau1 * a + b)/(c*tau1 + d)
            verbose("   computing integral from i*oo to %s using %s terms "
                    "and precision %s"%(tau2, T, prec), level=3)
            int2 = self._integration_to_tau(tau2, T, prec)
            ans = eN * (int1 - int2)
        else:
            tau1c = complex(-d/c, 1/c.abs())
            # computes the integral from tau to i*oo
            verbose("   computing integral from i*oo to %s using %s terms "
                    "and fast double precision"%(tau1c, T), level=3)
            int1c = self._integration_to_tau_double(tau1c, T)
            verbose("   yields %s "%int1c, level=3)
            tau2c = (tau1c * a + b)/(c*tau1c + d)
            verbose("   computing integral from i*oo to %s using %s terms "
                    "and fast double precision"%(tau2c, T), level=3)
            int2c = self._integration_to_tau_double(tau2c, T)
            verbose("   yields %s"%int2c, level=3)
            ans = eN * ComplexField(53)(int1c - int2c)
        return ans

#======= precise rationals =====

    # (key=lambda r,sign,use_partials:(r,sign)) lead to a compiler crash
    @cached_method
    def _value_ioo_to_r(self, Rational r, int sign=0,
                        int use_partials=2):
        r"""
        Return `[r]^+` or `[r]^-` for a rational `r`.

        INPUT:

        - ``r`` -- a rational number, which has to be a unitary cusp

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        - ``use_partials`` -- integer: 0 don't use partial summation to do
          the computation, 1 do use them, 2 (default) decide if it is
          meaningful to do so.

        OUTPUT: a rational number.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve("11a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M._value_ioo_to_r(0/1)
            1/5
            sage: M._value_ioo_to_r(3/11,-1)
            1/2

            sage: E = EllipticCurve("5077a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M._value_ioo_to_r(0/1)
            0
            sage: M._value_ioo_to_r(123/456)
            -3
            sage: M._value_ioo_to_r(-9/55,-1)
            -2
        """
        #verbose("       enter _value_ioo_to_r with r=%s, sign=%s"%(r,sign),
        #        level=5)
        cdef:
            double eps
            ComplexNumber la
            RealNumber lap
        if sign == 0:
            sign = self._global_sign
        if sign == 1:
            eps = self._eps_unitary_plus
        else:
            eps = self._eps_unitary_minus
        la = self._from_ioo_to_r_approx(r, eps, use_partials=use_partials)
        if sign == 1:
            lap = la.real()
        else:
            lap = la.imag()
        return self._round(lap, sign, True)

    @cached_method
    def _value_r_to_rr(self, Rational r, Rational rr, int sign=0,
                       int use_partials=2):
        r"""
        Return the rational number `[r']^+ - [r]^+`. However the
        computation may choose to work directly along the path from
        `r` to `r'` rather than going via `i\infty`, depending on what is
        faster.

        INPUT:

        - ``r`` and ``rr`` -- two rational number, both have to be
          uniraty cusps.

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        - ``use_partials`` -- integer: 0 don't use partial summation to do
          the computation, 1 do use them, 2 (default) decide if it is
          meaningful to do so.

        OUTPUT: a rational number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve("57a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M._value_r_to_rr(0/1,8/9)
            0
            sage: M._value_r_to_rr(0/1,8/9,use_partials=0)
            0
            sage: M._value_r_to_rr(17/19,-5/9)
            1/2
            sage: M._value_r_to_rr(14/19,-5/9,-1)
            1/2
            sage: M._value_r_to_rr(174/179,-53/91,-1)
            -1

            sage: E = EllipticCurve([91,127])
            sage: E.conductor()
            55196272
            sage: E.conductor().factor()
            2^4 * 3449767
            sage: M = E.modular_symbol(implementation="num")
            sage: M._value_r_to_rr(0/1,4/5)
            1
            sage: M._value_r_to_rr(7/11,14/75) # long time
            -1
        """
        #verbose("       enter _value_r_to_rr with r=%s, rr=%s,"
        #        " sign=%s"%(r,rr,sign), level=5)
        cdef:
            double eps
            ComplexNumber la
            RealNumber lap

        if sign == 0:
            sign = self._global_sign
        if sign == 1:
            eps = self._eps_unitary_plus
        else:
            eps = self._eps_unitary_minus
        la = self._from_r_to_rr_approx(r, rr, eps, use_partials=use_partials)
        if sign == 1:
            lap = la.real()
        else:
            lap = la.imag()
        return self._round(lap, sign, True)

    @cached_method
    def transportable_symbol(self, Rational r, Rational rr, int sign=0):
        r"""
        Return the symbol `[r']^+ - [r]^+` where `r'=\gamma(r)` for some
        `\gamma\in\Gamma_0(N)`. These symbols can be computed by transporting
        the path into the upper half plane close to one of the unitary cusps.
        Here we have implemented it only to move close to `i\infty` and `0`.

        INPUT:

        - ``r`` and ``rr`` -- two rational numbers

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        OUTPUT: a rational number

        EXAMPLES::

            sage: E = EllipticCurve("11a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M.transportable_symbol(0/1,-2/7)
            -1/2

            sage: E = EllipticCurve("37a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M.transportable_symbol(0/1,-1/19)
            0
            sage: M.transportable_symbol(0/1,-1/19,-1)
            0

            sage: E = EllipticCurve("5077a1")
            sage: M = E.modular_symbol(implementation="num")
            sage: M.transportable_symbol(0/1,-35/144)
            -3
            sage: M.transportable_symbol(0/1,-35/144,-1)
            0
            sage: M.transportable_symbol(0/1, -7/31798)
            0
            sage: M.transportable_symbol(0/1, -7/31798, -1)
            -5
        """
        #verbose("       enter transportable_symbol with r=%s, rr=%s,"
        #        " sign=%s"%(r,rr,sign), level=5)
        cdef:
            double eps
            ComplexNumber la
            RealNumber lap

        if sign == 0:
            sign = self._global_sign
        if sign == 1:
            eps = self._eps_unitary_plus
        else:
            eps = self._eps_unitary_minus
        la = self._transportable_approx(r, rr, eps)
        if sign == 1:
            lap = la.real()
        else:
            lap = la.imag()
        return self._round(lap, sign, True)

    #@cached_method
    def _symbol_non_unitary(self, Rational r, int sign=0):
        r"""
        Given a rational number `r`, this computes the modular symbol
        `[r]^+` or `[r]^-`. There is no assumption here on the cusp `r`,
        so a rather slow method via transportable paths is chosen. This
        should only be used for small denominators that are non unitary
        cusps.

        INPUT:

        - ``r`` -- a rational number representing a unitary cusp

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        OUTPUT: a rational number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("49a1"))
            sage: M._symbol_non_unitary(1/7)
            1/4
            sage: M._symbol_non_unitary(1/7,sign=-1)
            5/28

        Test for :trac:`28476` ::

            sage: M = ModularSymbolNumerical(EllipticCurve("361a1"))
            sage: M._symbol_non_unitary(1/19)
            5/19
        """
        #verbose("       enter _symbol_non_unitary with r=%s,"
        #        " sign=%s"%(r,sign), level=5)
        cdef:
            llong m, B, N_ell, aell, u, N = self._N_E
            Integer ell
            Rational r2, res

        if sign == 0:
            sign = self._global_sign

        rc = _CuspsForModularSymbolNumerical(r, N)
        r = rc._r
        m = rc._m
        B = llgcd(m, N)

        # find a prime congruent to 1 modulo B
        ell = Integer(B) + 1
        while llgcd(ell, N) != 1 or not ell.is_prime():
            ell += B
        if ell > self._lans:
            aell = self._E.ap(ell)
        else:
            aell = Integer(self._ans[ell])
        N_ell = ell + 1 - aell
        # {ell * r , r}
        verbose("     Compute symbol {ell*r -> r} = {%s -> %s}" % (ell*r, r),
                level=4)
        res = self.transportable_symbol(ell * r, r, sign=sign)
        # {(r + u)/ ell, r}
        u = Integer(0)
        while u < ell:
            r2 = (r+u) / ell
            verbose("     Compute symbol {r2-> r} = {%s -> %s}" % (r2, r),
                    level=4)
            res += self.transportable_symbol(r2, r, sign=sign)
            u += 1
        return -res/N_ell

    @cached_method # one call to manin_symbol will set between 4 and 8 values in fact
    def _manin_symbol_with_cache(self, llong u, llong v, int sign):
        r"""
        This helper function is called by manin_symbol below.
        The sign has to be 1 or -1, the coordinates u, v
        should already be in the normalised presentation.

        See the docstring of manin_symbol for detail,
        examples and indirect doctest for this method.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve('37a1')
            sage: ms = ModularSymbolNumerical(E)
            sage: ms.manin_symbol(4,17)
            1
            sage: ms.__cached_methods["_manin_symbol_with_cache"].cache # random
            {((1, 5, 1), ()): 1,
            ((1, 15, 1), ()): -1,
            ((1, 22, 1), ()): -1,
            ((1, 32, 1), ()): 1}
            sage: ms.__cached_methods["_manin_symbol_with_cache"].cache[(1,15,1),()]
            -1
            sage: ms.manin_symbol(4+17,-4)
            0
            sage: ms.__cached_methods["_manin_symbol_with_cache"].cache # random
            {((1, 4, 1), ()): 0,
            ((1, 5, 1), ()): 1,
            ((1, 8, 1), ()): 1,
            ((1, 9, 1), ()): 0,
            ((1, 14, 1), ()): -1,
            ((1, 15, 1), ()): -1,
            ((1, 22, 1), ()): -1,
            ((1, 23, 1), ()): -1,
            ((1, 28, 1), ()): 0,
            ((1, 29, 1), ()): 1,
            ((1, 32, 1), ()): 1,
            ((1, 33, 1), ()): 0}
            sage: ms.__cached_methods["_manin_symbol_with_cache"].cache[ (1,23,1), () ]
            -1
        """
        cdef:
            llong c, d, x, y, N = self._N_E, Mu, Mv, Qu, Qv, du=1, dv=1
            Rational r, rr, res

        #verbose("       enter _manin_symbol_with_cache with u=%s, v=%s,"
        #        " sign =%s"%(u,v,sign), level=5)

        if u == 0:
            verbose("   integrating from 0 to i*oo", level=3)
            r = Rational(0)
            return self._value_ioo_to_r(r, sign, use_partials=0)
        elif v == 0:
            verbose("   integrating from i*oo to 0", level=3)
            r = Rational(0)
            return - self._value_ioo_to_r(r, sign, use_partials=0)
        else:
            # (c:d) = (u:v) but c and d are fairly small
            # in absolute value
            Mu = llgcd(u,N)
            Qu = N/Mu
            Mv = llgcd(v,N)
            Qv = N/Mv
            isunitary = ( llgcd(Qu,Mu) == 1 and llgcd(Qv,Mv) == 1 )
            if isunitary:  # unitary case
                _ = best_proj_point(u, v, self._N_E, &c, &d)
            else:  # at least one of the two cusps is not unitary
                du = llgcd(Qu,Mu)
                dv = llgcd(Qv,Mv)
                NMM = N/Mv/Mu
                if dv == 1:
                    c = Mu
                    d = llinvmod(u/Mu, NMM)
                    d *= v
                    d = d % (N/Mu)
                    while llgcd(c,d) != 1:
                        d += N/Mu
                    d = d % N
                    # now (u:v) = (c:d) with c as small as possible.
                else:
                    d = Mv
                    c = llinvmod(v/Mv, NMM)
                    c *= u
                    c = c % (N/Mv)
                    while llgcd(c,d) != 1:
                        c += N/Mv
                    c = c % N
                    # now (u:v) = (c:d) with d as small as possible.
            #verbose("   better representant on P^1: "
            #        "(%s : %s)"%(c, d), level=3)
            # _, x, y = c.xgcd(d)
            _ = llxgcd(c, d, &x, &y)
            #if above != 1 or (c*v-u*d) % N != 0:
            #    print("BUG: ",u,v,c,d,Mu,Mv)
            x = x % N
            y = y % N
            # [[y -x], [c,d]] has det 1
            if c > 0:
                rr = Rational( (y,c) )
            else:
                rr = - Rational( (y, (-c)) )
            if d > 0:
                r = -Rational( (x,d) )
            else:
                r = Rational( (x, (-d)) )
            if isunitary:
                verbose("   integrate between %s and %s"%(r, rr), level=3)
                return self._value_r_to_rr(r, rr, sign, use_partials=2)
            else:
                if dv > 1:
                    res = self._symbol_non_unitary(r,sign)
                else:
                    res = self._value_ioo_to_r(r,sign, use_partials=2)
                if du > 1:
                    res -= self._symbol_non_unitary(rr,sign)
                else:
                    res -= self._value_ioo_to_r(rr,sign, use_partials=2)
                return res

    def manin_symbol(self, llong u, llong v, int sign=0):
        r"""
        Given a pair `(u,v)` presenting a point in
        `\mathbb{P}^1(\mathbb{Z}/N\mathbb{Z})` and hence a coset of
        `\Gamma_0(N)`, this computes the value of the Manin
        symbol `M(u:v)`.

        INPUT:

        - ``u`` -- an integer

        - ``v`` -- an integer such that `(u:v)` is a projective point
          modulo `N`

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: M = E.modular_symbol(implementation="num")
            sage: M.manin_symbol(1,3)
            -1/2
            sage: M.manin_symbol(1,3, sign=-1)
            -1/2
            sage: M.manin_symbol(1,5)
            1
            sage: M.manin_symbol(1,5)
            1

            sage: E = EllipticCurve('14a1')
            sage: M = E.modular_symbol(implementation="num")
            sage: M.manin_symbol(1,2)
            -1/2
            sage: M.manin_symbol(17,6)
            -1/2
            sage: M.manin_symbol(-1,12)
            -1/2
        """
        #verbose("       enter manin_symbol with u=%s, v=%s,"
        #        " sign =%s"%(u,v,sign), level=5)
        cdef:
            llong un, vn
            int oi
            Rational res, r2

        if sign == 0:
            sign = self._global_sign

        oi = proj_normalise(self._N_E, u, v, &un, &vn)
        #verbose("   normalized representant on P^1: "
        #        "(%s :%s)"%(un, vn), level=3)

        # is it already in the cache ?
        c = self.__cached_methods
        if "_manin_symbol_with_cache" in c:
            c = c["_manin_symbol_with_cache"]
            if c.is_in_cache(un,vn,sign):
                res = self._manin_symbol_with_cache(un,vn,sign)
                return res

        # the actual work in now done in a helper function
        # which does the correct caching
        res = self._manin_symbol_with_cache(un,vn,sign)

        # using the Manin relations coming from
        # complex conjugation, the matrices S and T in SL_2
        # we get for free additional values of Manin
        # symbols that we cache, too.
        # This sets 6 values in average
        c = self.__cached_methods["_manin_symbol_with_cache"]

        # (-v:u) = - (u:v)
        oi = proj_normalise(self._N_E, -v, u, &un, &vn)
        c.set_cache(-res, un, vn, sign)

        # (v:u) = -1 * sign * (u:v)
        oi = proj_normalise(self._N_E, v, u, &un, &vn)
        c.set_cache(-sign*res, un, vn, sign)

        # (-u:v) = sign * (u:v)
        oi = proj_normalise(self._N_E, -u, v, &un, &vn)
        c.set_cache(sign*res, un, vn, sign)

        # (u:v) + ( u+v:-u) +(v,-u-v) = 0
        # is ( u+v:-u) already computed, we set the third
        oi = proj_normalise(self._N_E, u+v, -u, &un, &vn)
        if c.is_in_cache(un,vn,sign):
            r2 = - res - c(un,vn,sign)

            # (v:-u-v) = r2
            oi = proj_normalise(self._N_E, v, -u-v, &un, &vn)
            c.set_cache(r2, un, vn, sign)

            # (u+v:v) = -r2
            oi = proj_normalise(self._N_E, u+v, v, &un, &vn)
            c.set_cache(-r2, un, vn, sign)

            # (-u-v:v) = -1 * sign * r2
            oi = proj_normalise(self._N_E, -u-v, v, &un, &vn)
            c.set_cache(-sign*r2, un, vn, sign)

            # (-v:-u-v) = sign * r2
            oi = proj_normalise(self._N_E, -v, -u-v, &un, &vn)
            c.set_cache(sign*r2, un, vn, sign)

        # is ( v,-u-v) already computed, we set ( u+v:-u)
        oi = proj_normalise(self._N_E, v, -u-v, &un, &vn)
        if c.is_in_cache(un,vn,sign):
            r2 = - res - c(un,vn,sign)

            # (u+v:-u) = r2
            oi = proj_normalise(self._N_E, u+v, -u, &un, &vn)
            c.set_cache(r2, un, vn, sign)

            # (u:u+v) = -r2
            oi = proj_normalise(self._N_E, u, u+v, &un, &vn)
            c.set_cache(-r2, un, vn, sign)

            # (-u:u+v) = -1 * sign * r2
            oi = proj_normalise(self._N_E, -u, u+v, &un, &vn)
            c.set_cache(-sign*r2, un, vn, sign)

            # (-u-v:-u) = sign * r2
            oi = proj_normalise(self._N_E, -u-v, -u, &un, &vn)
            c.set_cache(sign*r2, un, vn, sign)

        return res

    # ===============================

    @cached_method  # not sure this is not a waist
    def _evaluate(self, Rational r, int sign=0):
        r"""
        Given a rational number `r` this computes the modular symbol
        `[r]^+` or `[r]^-`.

        INPUT:

        - ``r`` -- a rational number representing a unitary cusp

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        OUTPUT: a rational number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve('5077a1')
            sage: M = ModularSymbolNumerical(E)
            sage: M._evaluate(1/7)
            3
            sage: M._evaluate(1/7, -1)
            1
            sage: M._evaluate(-1/4112)
            3

            sage: E = EllipticCurve('11a1')
            sage: M = ModularSymbolNumerical(E)
            sage: M._evaluate(1/99999)
            -4/5

            sage: M = ModularSymbolNumerical( EllipticCurve("32a1") )
            sage: M._evaluate(3/5)
            -1/4

        Non-unitary examples::

            sage: M = ModularSymbolNumerical(EllipticCurve("20a1"))
            sage: M._evaluate(1/2)
            -1/6
            sage: M._evaluate(1/2, -1)
            0
            sage: M = ModularSymbolNumerical(EllipticCurve("49a1"))
            sage: M._evaluate(2/7)
            -1/4
            sage: M = ModularSymbolNumerical(EllipticCurve("121a1"))
            sage: M._evaluate(3/11)
            1/2

        This takes a bit longer::

            sage: M = ModularSymbolNumerical(EllipticCurve("78a1"))
            sage: M._evaluate(1/38)
            2
            sage: M._evaluate(5/38)
            1/2
            sage: M._evaluate(1/123456789012345678901234567)
            -1
        """
        #verbose("       enter _evaluate with r=%s, sign=%s"%(r,sign),
        #         level=5)
        cdef:
            llong N = self._N_E, u, v
            Rational r2, res
            Integer a, m, B, Q, y, x, M, uu, vv

        if sign == 0:
            sign = self._global_sign

        a = r.numerator()
        m = r.denominator()
        a = a % m
        if 2*a > m:
            a -= m
        r = a/m

        if N % 4 == 0 and m % 4 != 0 and m % 2 == 0:
            # this follows from T_2
            r2 = (a-m/2)/m
            return -self._evaluate(r2, sign=sign)

        B = m.gcd(N)
        Q = N // B
        #verbose("     cusp is %s/%s of width %s"%(a,m,Q), level=4)

        if r == 0:
            return self._value_ioo_to_r(r, sign=sign)

        if m < self._cut_val:
            # now at some point we go directly to ioo
            M = N//Q
            if Q.gcd(M) == 1:
                res = self._value_ioo_to_r(r, sign=sign)
            else:
                res = self._symbol_non_unitary(r, sign=sign)
        else:
            # so a*y + x*m = 1 and 0 <= |y| < m/2
            _, y, x = a.xgcd(m)
            y = y % m
            if 2*y > m:
                y -= m
            x = (1-y*a) // m
            #verbose("     smallest xgcd is "
            #        + " %s = %s * %s + %s * %s"%(a.gcd(m),a,y,x,m),
            #        level=4)
            # make the cusp -x/y unitary if possible.
            B = y.gcd(N)
            if B.gcd(N//B) != 1:
                if y > 0:
                    y -= m
                    x += a
                else:
                    y += m
                    x -= a
            # Note: it could still be non-unitary.
            # Example: N=36 a=2, m=5
            uu = (-y) % N
            vv = m % N
            u = <llong>uu
            v = <llong>vv
            r2 = - x/y
            verbose("  Next piece: integrate from %s to %s via the Manin"
                    " symbol for (%s : %s)"%(r,r2,u,v), level=2)
            res = self.manin_symbol(u,v,sign=sign)
            res += self._evaluate(r2, sign=sign)

        return res

    @cached_method
    def all_values_for_one_denominator(self, llong m, int sign=0):
        r"""
        Given an integer ``m`` and a ``sign``, this returns the
        modular symbols `[a/m]` for all `a` coprime to `m`
        using partial sums.
        This is much quicker than computing them one by one.

        This will only work if `m` is relatively small and
        if the cusps `a/m` are unitary.

        INPUT:

        - ``m`` -- a natural number

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        OUTPUT: a dictionary of fractions with denominator `m`
        giving rational numbers.

        EXAMPLES::

            sage: E = EllipticCurve('5077a1')
            sage: M = E.modular_symbol(implementation="num")
            sage: M.all_values_for_one_denominator(7)
            {1/7: 3, 2/7: 0, 3/7: -3, 4/7: -3, 5/7: 0, 6/7: 3}
            sage: [M(a/7) for a in [1..6]]
            [3, 0, -3, -3, 0, 3]
            sage: M.all_values_for_one_denominator(3,-1)
            {1/3: 4, 2/3: -4}

            sage: E = EllipticCurve('11a1')
            sage: M = E.modular_symbol(implementation="num")
            sage: M.all_values_for_one_denominator(12)
            {1/12: 1/5, 5/12: -23/10, 7/12: -23/10, 11/12: 1/5}
            sage: M.all_values_for_one_denominator(12, -1)
            {1/12: 0, 5/12: 1/2, 7/12: -1/2, 11/12: 0}

            sage: E = EllipticCurve('20a1')
            sage: M = E.modular_symbol(implementation="num")
            sage: M.all_values_for_one_denominator(4)
            {1/4: 0, 3/4: 0}
            sage: M.all_values_for_one_denominator(8)
            {1/8: 1/2, 3/8: -1/2, 5/8: -1/2, 7/8: 1/2}
        """

        # TODO: This function should really use some multi-radix FFT
        # however for small m there is not much gain
        cdef:
            llong N, Q, z, a, astar, j, epsQ
            double resam, twopim
            RealNumber val
            dict res

        #verbose("       enter all_symbol with m=%s"%m, level=5)

        if sign == 0:
            sign = self._global_sign

        RR = RealField(53)
        N = self._N_E
        Q = N / llgcd(m,N)
        if llgcd(m,Q) > 1:
            raise NotImplementedError("Only implemented for cusps that are "
                                      "in the Atkin-Lehner orbit of oo")
        #verbose("   compute all partial sums with denominator m=%s"%m,
        #        level=3)
        z = Q*m*m
        v = self._kappa(m, z)

        epsQ = self._epsQs[Q]

        res = dict()
        a = 1
        twopim = TWOPI
        twopim = twopim/m
        if sign == 1:
            while a < m:
                if llgcd(a,m) == 1:
                    astar = llinvmod(Q*a, m)
                    j = 0
                    resam = 0
                    while j < m:
                        resam +=  (v[j]*(cos(twopim*j*a)
                                   - epsQ*cos(twopim*j*astar)) )
                        j += 1
                    val = RR(resam)
                    res[Rational( (a,m) )] = self._round(val, 1, True)
                a += 1
        else:
            while a < m:
                if llgcd(a,m) == 1:
                    astar = llinvmod(Q*a, m)
                    j = 0
                    resam = 0
                    while j < m:
                        resam +=  (v[j]*(sin(twopim *j *a)
                                   + epsQ *sin(twopim*j*astar)) )
                        j += 1
                    val = RR(resam)
                    res[Rational( (a,m) )] = self._round(val, -1, True)
                a += 1

        return res

    def _twisted_symbol(self, Rational ra, int sign=0):
        r"""
        Compute the value of the modular symbol by
        using the symbols of the quadratic twist.

        INPUT:

        - ``ra`` -- a rational number

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        OUTPUT: a rational number

        EXAMPLES::

            sage: E = EllipticCurve("735e4")
            sage: M = E.modular_symbol(implementation="num")
            sage: M(1/19, sign=-1, use_twist=False) #indirect doctest
            4
            sage: M(1/19, sign=-1, use_twist=True)
            4
            sage: M(6/19, sign=1, use_twist=False)
            3
            sage: M(6/19, sign=1, use_twist=True)
            3
        """
        cdef Integer D, Da, a
        cdef Rational res, t
        #verbose("       enter _twisted symbol with ra=%s,
        #        sign=%s"%(ra,sign),
        #        level=5)
        if sign == 0:
            sign = self._global_sign
        if self._D == -1:
            self._set_up_twist()
        D = self._D
        Da = D.abs()
        a = Integer(1)
        res = Rational( (0,1) )
        s = sign * D.sign()
        verbose("     start sum of twisted symbols with disc %s"%D, level=4)
        while a < Da:
            if a.gcd(Da) == 1:
                t = self._Mt(ra - a/Da, sign=s, use_twist=False)
                res += kronecker_symbol(D,a) * t
            a += 1
        res = res/ self._twist_q
        if sign == 1 and D < 0:
            res = -res
        return res

#====================== approximative versions

    def _evaluate_approx(self, Rational r, double eps):
        r"""
        Given a rational number `r` this computes the integral
        `\lambda(r)` with maximal error ``eps``.

        INPUT:

        - ``r`` -- a rational number representing a unitary cusp

        - ``eps`` -- a positive real number

        OUTPUT: a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: E = EllipticCurve('5077a1')
            sage: m = ModularSymbolNumerical(E)
            sage: m._evaluate_approx(1/11,0.000001)  # abs tol 1e-11
            9.69540669970570e-10 - 5.80486769763411e-11*I
            sage: m._evaluate_approx(1/17,0.000001) # abs tol 1e-11
            -9.01145713605445e-10 + 7.40274134212215*I

            sage: M = ModularSymbolNumerical( EllipticCurve([-12,79]) )
            sage: M.elliptic_curve().conductor()
            287280
            sage: M._evaluate_approx(0/1,0.01)  # abs tol 1e-11
            0.000000000000000
            sage: M._evaluate_approx(1/17,0.01) # abs tol 1e-11
            1.08712572498569 - 0.548379313090719*I

       Test that is also works for non-unitary cusps (:trac:`29476`) ::

            sage: E = EllipticCurve("20a1")
            sage: m = E.modular_symbol_numerical()
            sage: m(1/2)          #abs tol 1e-4
            -0.166666666666667
        """
        #verbose("       enter _evaluate_approx with r=%s, eps=%s"%(r,eps),
        #        level=5)
        cdef:
            llong N = self._N_E
            ComplexNumber res
            Rational r2
            Integer a, m, Q, x, y, B, M

        a = r.numerator()
        m = r.denominator()
        a = a % m
        if 2*a > m:
            a -= m
        r = a/m
        B = m.gcd(N)
        Q = N // B
        verbose("     cusp is %s/%s of width %s"%(a,m,Q), level=4)

        if r == 0:
            return self._from_ioo_to_r_approx(r, eps, use_partials=0)

        M = N//Q
        if Q.gcd(M) != 1:
            return self._symbol_non_unitary_approx(r, eps)

        if m < self._cut_val:
            # now at some point we go directly to ioo
            return self._from_ioo_to_r_approx(r, eps, use_partials=0)

        _, y, x = a.xgcd(m)
        y = y % m
        if 2*y > m:
            y -= m
        x = (1-y*a) // m
        #verbose("     smallest xgcd is "
        #        + " %s = %s * %s + %s * %s"%(a.gcd(m),a,y,x,m),
        #        level=4)
        # make the cusp -x/y unitary if possible.
        B = y.gcd(N)
        if B.gcd(N // B) != 1:
            if y > 0:
                y -= m
                x += a
            else:
                y += m
                x -= a
        r2 = - x / y
        B = y.gcd(N)
        Q = N // B
        if Q.gcd(N // Q) != 1: # r2 is not unitary
            return self._symbol_non_unitary_approx(r, eps)

        r2 = - x / y
        verbose("Next piece: integrate to the cusp %s " % r2, level=2)
        res = self._from_r_to_rr_approx(r, r2, eps,
                                        use_partials=2)
        res += self._evaluate_approx(r2, eps)
        return res

    def _symbol_non_unitary_approx(self, Rational r, double eps):
        r"""
        Given a rational number `r` this computes the integral
        `\\lambda(r)` to maximal error ``eps``.

        There is no assumption here on the cusp `r`,
        so a rather slow method via transportable paths is chosen. For
        unitary cusps please use ``_symbol_unitary_approx``.

        INPUT:

        - ``r`` -- a rational number representing a unitary cusp

        - ``eps`` -- a positive real

        OUTPUT: a complex number

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.mod_sym_num \
            ....: import ModularSymbolNumerical
            sage: M = ModularSymbolNumerical(EllipticCurve("20a1"))
            sage: M._symbol_non_unitary_approx(1/2,0.0001) # abs tol 1e-11
            -0.470729190326520 + 2.59052039079203e-16*I

            sage: M = ModularSymbolNumerical(EllipticCurve("49a1"))
            sage: M._symbol_non_unitary_approx(2/7,0.000000001) # abs tol 1e-11
            -0.483327926404308 + 0.548042354981878*I

        A bit longer take::

            sage: M = ModularSymbolNumerical(EllipticCurve("78a1"))
            sage: M._symbol_non_unitary_approx(1/38,0.1)  # abs tol 1e-11
            2.90087559068021 + 2.86538550720028e-7*I
            sage: M._symbol_non_unitary_approx(5/38,0.1)  # abs tol 1e-11
            0.725215164486092 - 1.19349741385624*I
         """
        #verbose("       enter _symbol_nonunitary_approx with r=%s,"
        #        " eps=%s"%(r,eps), level=5)
        cdef:
            llong m, B, N_ell, aell, u, N = self._N_E
            Integer ell
            Rational r2
            ComplexNumber res

        rc = _CuspsForModularSymbolNumerical(r, N)
        r = rc._r
        m = rc._m
        B = llgcd(m, N)

        # find a prime congruent to 1 modulo B
        ell = Integer(B) + 1
        while llgcd(ell, N) != 1 or not ell.is_prime():
            ell += B
        if ell > self._lans:
            aell = self._E.ap(ell)
        else:
            aell = Integer(self._ans[ell])
        N_ell = ell + 1 - aell
        # {ell * r , r}
        verbose("     Compute symbol {ell*r -> r} = {%s -> %s}"%(ell*r,r),
                level=4)
        res = self._transportable_approx(ell * r, r, eps)
        # {(r + u)/ ell, r}
        u = Integer(0)
        while u < ell:
            r2 = (r+u) / ell
            verbose("     Compute symbol {r2-> r} = {%s -> %s}"%(r2,r),
                    level=4)
            res += self._transportable_approx(r2, r, eps)
            u += 1
        return -res/N_ell

    def _twisted_approx(self, Rational ra, int sign=0, int prec=20):
        r"""
        Compute the approximate value of the modular
        symbol by using the symbols of the quadratic twist.

        Note that _set_up_twist needs to be called first
        and D must be different from 1.

        INPUT:

        - ``ra`` -- a rational number

        - ``sign`` -- optional either +1 or -1, or 0 (default),
          in which case the sign passed to the class is taken.

        - ``prec`` -- an integer (default 20)

        OUTPUT: a real number

        EXAMPLES::

            sage: E = EllipticCurve("735e4")
            sage: M = E.modular_symbol(implementation="num")
            sage: M.approximative_value(1/19, sign=-1, prec=20, use_twist=False) # indirect doctest abs tol 1e-11
            4.00000000089736
            sage: M.approximative_value(1/19, sign=-1, prec=20, use_twist=True) # abs tol 1e-11
            3.99999999982043

            sage: M.approximative_value(6/19, sign=1, prec=20, use_twist=False) # abs tol 1e-11
            2.99999999944834
            sage: M.approximative_value(6/19, sign=1, prec=20, use_twist=True) # abs tol 1e-11
            3.00000000021937
        """
        cdef Integer D, Da, a, s, precd
        cdef RealNumber res, t
        #verbose("       enter _twisted approx with ra=%s,
        #        eps=%s"%(ra,eps),
        #        level=5)

        if sign == 0:
            sign = self._global_sign
        D = self._D
        s = sign * D.sign()
        Da = D.abs()
        precd = prec + euler_phi(Da).log(2,20).ceil()
        a = Integer(1)
        res = self._Mt.approximative_value(ra - a/Da, s, precd)
        verbose("     start sum of twisted symbols with disc %s"%D, level=4)
        a += 1
        while a < Da:
            if a.gcd(Da) == 1:
                t = self._Mt.approximative_value(ra - a/Da, s, precd, use_twist=False)
                res += kronecker_symbol(D,a) * t
            a += 1
        res = res/self._twist_q
        if sign == 1 and D < 0:
            res = -res
        return res


#==========================
# Doctest functions for the above class

def _test_init(E):
    r"""
    Doctest function for the initialisation of
    ModularSymbolNumerical.

    INPUT:

    - ``E`` -- an elliptic curve

    OUTPUT:

    - a dictionary of eigenvalues for the Atkin-Lehner involutions

    - five integers representing the Fourier coefficients `a_1`,
      `a_2`, `a_3`, `a_389` and `a_2013`. (This will test
      _add_an_coefficients).

    - four integers that are bounds for the denominators as in
      _set_den_bounds

    - four real numbers which are allowed errors in computations

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num import _test_init
        sage: _test_init(EllipticCurve("11a1")) # abs tol 1e-11
        ({1: 1, 11: -1}, [1, -2, -1, -15, -12], [10, 2, 10, 2],
        [0.06346046521397766,
        0.7294083084692475,
        0.06346046521397766,
        0.7294083084692475])
        sage: _test_init(EllipticCurve("11a2")) # abs tol 1e-11
        ({1: 1, 11: -1}, [1, -2, -1, -15, -12], [2, 2, 2, 2],
        [0.06346046521397766,
        0.7294083084692475,
        0.06346046521397766,
        0.7294083084692475])
        sage: _test_init(EllipticCurve("11a3")) # abs tol 1e-11
        ({1: 1, 11: -1}, [1, -2, -1, -15, -12], [50, 2, 50, 2],
        [0.06346046521397768,
        0.7294083084692478,
        0.06346046521397768,
        0.7294083084692478])

        sage: _test_init(EllipticCurve("14a6")) # abs tol 1e-11
        ({1: 1, 2: 1, 7: -1, 14: -1},
         [1, -1, -2, 18, 0],
         [9, 1, 9, 1],
         [0.16511182967224025,
          0.6627456198412432,
          0.16511182967224025,
          0.6627456198412432])

        sage: _test_init(EllipticCurve("20a1")) # abs tol 1e-11
        ({1: 1, 2: -1, 4: -1, 5: 1, 10: -1, 20: -1},
        [1, 0, -2, -6, 0], [48, 48, 12, 2],
        [0.029420574395407434,
        0.023689220823344594,
        0.11768229758162974,
        0.5685412997602702])
        sage: _test_init(EllipticCurve("37a1")) # abs tol 1e-11
        ({1: 1, 37: 1}, [1, -2, -3, 4, -120], [1, 2, 1, 2],
        [1.4967293231159797,
        0.6128473454966975,
        1.4967293231159797,
        0.6128473454966975])

        sage: E = EllipticCurve([91,127])
        sage: E.conductor().factor()
        2^4 * 3449767
        sage: _test_init(E) # abs tol 1e-11
        ({1: 1,
          2: -1,
          4: -1,
          8: -1,
          16: -1,
          3449767: 1,
          6899534: -1,
          13799068: -1,
          27598136: -1,
          55196272: -1},
         [1, 0, 0, 2, 0],
         [4, 4, 2, 2],
         [0.15583810484385163,
          0.14150261234359426,
          0.31167620968770327,
          0.28300522468718853])
     """
    M = ModularSymbolNumerical(E)
    e = M._epsQs
    a1 = Integer(M._ans[1])
    a2 = Integer(M._ans[2])
    a3 = Integer(M._ans[3])
    a4 = Integer(M._ans[389])
    M._add_an_coefficients(2014)
    a5 = Integer(M._ans[2013])
    t1 = Integer(M._t_plus)
    t2 = Integer(M._t_minus)
    t3 = Integer(M._t_unitary_plus)
    t4 = Integer(M._t_unitary_minus)
    e1 = M._eps_plus
    e2 = M._eps_minus
    e3 = M._eps_unitary_plus
    e4 = M._eps_unitary_minus
    return e, [a1, a2, a3, a4, a5], [t1, t2, t3, t4], [e1, e2, e3, e4]


def _test_integration(E, a, b, T):
    r"""
    Doctest for the numerical integration in
    _integration_to_tau_double in the above class.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``a`` and ``b`` -- two real numbers representing real and
      imaginary part of a complex number tau

    - ``T `` -- an integer for the number of terms to use

    OUTPUT: a complex number

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import _test_integration
        sage: E = EllipticCurve("11a1")
        sage: _test_integration(E, 0,0.01,1000) # abs tol 1e-11
        (0.2538418608559108+0j)
        sage: _test_integration(E, 0,0.0001,10000) # abs tol 1e-11
        (0.2538815728257322+0j)

        sage: E = EllipticCurve("37a1")
        sage: _test_integration(E, 0, 0.0001,1000) # abs tol 1e-11
        (-0.0105693920159094+0j)
        sage: _test_integration(E, 0.7, 0.1, 10000) # abs tol 1e-11
        (-0.021614803690068213-0.7770316490609953j)
        sage: _test_integration(E, 0.7, 0.1, 20000) # abs tol 1e-11
        (-0.021614803690068213-0.7770316490609953j)
    """
    M = ModularSymbolNumerical(E)
    c = complex(a,b)
    tt = <int>T
    ans = M._integration_to_tau_double(c,tt)
    return ans


def _test_integration_via_partials(E, y, m, T):
    r"""
    Doctest for the numerical integration in
    _partial_real_sums_double in the above class.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``y`` -- a real number

    - ``m`` -- an integers

    - ``T `` -- an integer for the number of terms to use

    OUTPUT: a list of `m` real numbers

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import _test_integration_via_partials
        sage: E = EllipticCurve("11a1")
        sage: _test_integration_via_partials(E,0.001,3,1000) # abs tol 1e-11
        [-0.16916415619939476, 1.0536872023214188, -0.6306661264594561]

        sage: E = EllipticCurve("121c1")
        sage: _test_integration_via_partials(E,0.03,3,700) # abs tol 1e-11
        [0.49198993741342784, 0.6601504274130793, 0.3177042713926389]
        sage: _test_integration_via_partials(E,0.03,3,7000)  # abs tol 1e-11
        [0.49198993741342784, 0.6601504274130793, 0.3177042713926389]
        """
    cdef int mm = <int>(m)
    cdef double * ra
    ra = <double *> sig_malloc( mm * sizeof(double))
    if ra is NULL:
        raise MemoryError
    M = ModularSymbolNumerical(E)
    _ = M._partial_real_sums_double(y, m, T, ra)
    res = [ra[j] for j in range(m)]
    sig_free(ra)
    return res


def _test_against_table(range_of_conductors, other_implementation="sage", list_of_cusps=[], verb=False):
    r"""
    This test function checks the modular symbols here against the
    ones implemented already. Note that for some curves the current
    implementation does not scale them correctly and so we might be off
    by a small integer.

    INPUT:

    - ``range_of_conductors`` -- a list of integers; all curves with
      conductors in that list will be tested.

    - ``list_of_cusps`` -- a list of rationals to be tested

    - ``verb`` -- if True (default) prints the values

    OUTPUT: Boolean. If False the function also prints information.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.mod_sym_num \
        ....: import _test_against_table
        sage: _test_against_table([11,37]) # long time
        True
    """
    boo = True
    from sage.schemes.elliptic_curves.ell_rational_field import cremona_curves
    for C in cremona_curves(range_of_conductors):
        if verb:
            print("testing curve ", C.label())
        m = C.modular_symbol(implementation=other_implementation)
        m2 = C.modular_symbol(sign=-1, implementation=other_implementation)
        M = ModularSymbolNumerical(C)
        # a few random small rationals
        if len(list_of_cusps)==0:
            list_of_cusps = [Rational((0,1)),Rational((1,1)),Rational((1,2)),
                             Rational((1,3)),Rational((1,4)),Rational((2,5)),
                             Rational((1,6)),Rational((3,7)),Rational((1,8)),
                             Rational((5,9)),Rational((7,10))]
        for r in list_of_cusps:
            mr = m(r)
            m2r = m2(r)
            Mr = M(r)
            M2r = M(r, sign=-1)
            if verb:
                print("r={} : ({},{}),({}, {})".format(r,mr,m2r,Mr,M2r), end="  ", flush=True)
            if mr != Mr or m2r != M2r:
                print (("B u g : curve = {}, cusp = {}, sage's symbols"
                        + "({},{}), our symbols ({}, {})").format(C.label(), r,
                                                                  mr, m2r, Mr, M2r) )
                boo = False
        M.clear_cache()
    return boo


# =========================================
#
# Justification for numbers in _get_truncation_and_prec
# in the estimates of sigma_0(n)
# the following code in sage gives n_0 such that
# sigma0(n) < B * sqrt(n) for all n> n0 for a given B
#
#y = lambda n: number_of_divisors(n)/sqrt(n*1.)
#
#def hidef(B):
#    """finds all n with y(n) > B for some bound B"""
#    li = [1]
#    old = [1]
#    while old:
#        new = []
#        p = 1
#        boo = True
#        while boo:
#            p = next_prime(p)
#            boo = False
#            for n in old:
#                m = n*p
#                if m not in new and y(m) > B:
#                    new.append(m)
#                    boo = True
#        li += new
#        old = new
#    return li
#
#def last_hidef(B):
#    n = max(hidef(B))
#    return n, y(n)
#
#for B in [1,2/3,1/2,1/3,1/4,1/5,1/6]:
#    print(last_hidef(B))
#
#(1260, 1.01418510567422)
#(10080, 0.717137165600636)
#(55440, 0.509647191437626)
#(277200, 0.341881729378914)
#(831600, 0.263180677983908)
#(2162160, 0.217623636951613)
#(4324320, 0.184659779321958)
