"""
Basic arithmetic with c-integers.
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************


###################################################################
# We define the following functions in this file, both
# for int (up to bound = 2**31 - 1) and longlong (up to 2**63 - 1).
# The function definitions are identical except for the types.
# Some of their input can be at most sqrt(bound), since
# it is necessary to multiply numbers and reduce the product
# modulo n, where n is at most bound.
#
#   * abs_int -- absolute value of integer
#   * sign_int -- sign of integer
#   * c_gcd_int -- gcd of two ints
#   * gcd_int -- python export of c_gcd_int
#   * c_xgcd_int -- extended gcd of two ints
#   * c_inverse_mod_int -- inverse of an int modulo another int
#   * c_rational_recon_int -- rational reconstruction of ints
#   * rational_recon_int -- export of rational reconstruction for ints
#
#  The long long functions are the same, except they end in _longlong.
#
###################################################################

# The int definitions

include "sage/ext/gmp.pxi"
include "sage/ext/stdsage.pxi"
include "sage/libs/pari/decl.pxi"

cdef extern from "pari/pari.h":
    cdef long NEXT_PRIME_VIADIFF(long, unsigned char*)

from sage.rings.integer_ring import ZZ
from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.all import pari
from sage.rings.integer cimport Integer

cpdef prime_range(start, stop=None, algorithm="pari_primes", bint py_ints=False):
    r"""
    List of all primes between start and stop-1, inclusive.  If the
    second argument is omitted, returns the primes up to the first
    argument.

    This function is closely related to (and can use) the primes iterator.
    Use algorithm "pari_primes" when both start and stop are not too large,
    since in all cases this function makes a table of primes up to
    stop. If both are large, use algorithm "pari_isprime" instead.

    Algorithm "pari_primes" is faster for most input, but crashes for larger input.
    Algorithm "pari_isprime" is slower but will work for much larger input.

    INPUT:

        - ``start`` -- lower bound

        - ``stop`` -- upper bound

        - ``algorithm`` -- string, one of:

             - "pari_primes": Uses PARI's primes function.  Generates all primes up to stop.
                              Depends on PARI's primepi function.

             - "pari_isprime": Uses a mod 2 wheel and PARI's isprime function by calling
                             the primes iterator.

        - ``py_ints`` -- boolean (default False), return Python ints rather than Sage Integers (faster)


    EXAMPLES:
        sage: prime_range(10)
        [2, 3, 5, 7]
        sage: prime_range(7)
        [2, 3, 5]
        sage: prime_range(2000,2020)
        [2003, 2011, 2017]
        sage: prime_range(2,2)
        []
        sage: prime_range(2,3)
        [2]
        sage: prime_range(5,10)
        [5, 7]
        sage: prime_range(-100,10,"pari_isprime")
        [2, 3, 5, 7]
        sage: prime_range(2,2,algorithm="pari_isprime")
        []
        sage: prime_range(10**16,10**16+100,"pari_isprime")
        [10000000000000061, 10000000000000069, 10000000000000079, 10000000000000099]
        sage: prime_range(10**30,10**30+100,"pari_isprime")
        [1000000000000000000000000000057, 1000000000000000000000000000099]
        sage: type(prime_range(8)[0])
        <type 'sage.rings.integer.Integer'>
        sage: type(prime_range(8,algorithm="pari_isprime")[0])
        <type 'sage.rings.integer.Integer'>

    TESTS:
        sage: len(prime_range(25000,2500000))
        180310
        sage: prime_range(2500000)[-1].is_prime()
        True

    AUTHORS:
      - William Stein (original version)
      - Craig Citro (rewrote for massive speedup)
      - Kevin Stueve (added primes iterator option) 2010-10-16
      - Robert Bradshaw (speedup using Pari prime table, py_ints option)
    """
    cdef Integer z
    cdef long c_start, c_stop, p
    cdef unsigned char* pari_prime_ptr
    if algorithm == "pari_primes":
        if stop is None:
            # In this case, "start" is really stop
            c_start = 1
            c_stop = start
        else:
            c_start = start
            c_stop = stop
            if c_stop <= c_start:
                return []
            if c_start < 1:
                c_start = 1
        if maxprime() < c_stop:
            pari.init_primes(c_stop)
        pari_prime_ptr = <unsigned char*>diffptr
        p = 0
        res = []
        while p < c_start:
            NEXT_PRIME_VIADIFF(p, pari_prime_ptr)
        while p < c_stop:
            if py_ints:
                res.append(p)
            else:
                z = <Integer>PY_NEW(Integer)
                mpz_set_ui(z.value, p)
                res.append(z)
            NEXT_PRIME_VIADIFF(p, pari_prime_ptr)

    elif algorithm == "pari_isprime":
        from sage.rings.arith import primes
        res = list(primes(start, stop))
    else:
        raise ValueError("algorithm argument must be either ``pari_primes`` or ``pari_isprime``")
    return res

cdef class arith_int:
    cdef public int abs_int(self, int x) except -1:
        if x < 0:
            return -x
        return x

    cdef public int sign_int(self, int n) except -2:
        if n < 0:
            return -1
        return 1

    cdef public int c_gcd_int(self, int a, int b) except -1:
        cdef int c
        if a==0:
            return self.abs_int(b)
        if b==0:
            return self.abs_int(a)
        if a<0: a=-a
        if b<0: b=-b
        while(b):
            c = a % b
            a = b
            b = c
        return a


    def gcd_int(self, int a, int b):
        return self.c_gcd_int(a,b)


    cdef public int c_xgcd_int(self, int a, int b, int* ss, int* tt) except -1:
        cdef int psign, qsign, p, q, r, s, c, quot, new_r, new_s

        if a == 0:
            ss[0] = 0
            tt[0] = self.sign_int(b)
            return self.abs_int(b)

        if b == 0:
            ss[0] = self.sign_int(a)
            tt[0] = 0
            return self.abs_int(a)

        psign = 1; qsign = 1

        if a<0: a = -a; psign = -1
        if b<0: b = -b; qsign = -1

        p = 1; q = 0; r = 0; s = 1
        while (b):
            c = a % b; quot = a/b
            a = b; b = c
            new_r = p - quot*r
            new_s = q - quot*s
            p = r; q = s
            r = new_r; s = new_s

        ss[0] = p*psign
        tt[0] = q*qsign

        return a

    def xgcd_int(self, int a, int b):
        cdef int g, s, t
        g = self.c_xgcd_int(a,b, &s, &t)
        return (g,s,t)

    cdef public int c_inverse_mod_int(self, int a, int m) except -1:
        if a == 1 or m<=1: return a%m   # common special case
        cdef int g, s, t
        g = self.c_xgcd_int(a,m, &s, &t)
        if g != 1:
            raise ArithmeticError, "The inverse of %s modulo %s is not defined."%(a,m)
        s = s % m
        if s < 0:
            s = s + m
        return s


    def inverse_mod_int(self, int a, int m):
        return self.c_inverse_mod_int(a, m)

    cdef int c_rational_recon_int(self, int a, int m, int* n, int* d) except -1:
        cdef int u, v, u0, u1, u2, v0, v1, v2, q, t0, t1, t2, x, y
        cdef float bnd

        if m>46340:
            raise OverflowError, "The modulus m(=%s) should be at most 46340"%m
            return -1

        a = a % m

        if a==0 or m == 0:
            n[0] = 0
            d[0] = 1
            return 0

        if m<0: m = -m
        if a<0: a = m - a
        if a==1:
            n[0] = 1
            d[0] = 1
            return 0

        u = m
        v = a
        bnd = sqrt(m/2.0)
        u0=1; u1=0; u2=u
        v0=0; v1=1; v2=v
        while self.abs_int(v2) > bnd:
            q = u2/v2   # floor is implicit
            t0=u0-q*v0; t1=u1-q*v1; t2=u2-q*v2
            u0=v0; u1=v1; u2=v2
            v0=t0; v1=t1; v2=t2;

        x = self.abs_int(v1); y = v2
        if v1<0:  y = -1*y
        if x<=bnd and self.c_gcd_int(x,y)==1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_int(self, int a, int m):
        """
        Rational reconstruction of a modulo m.
        """
        cdef int n, d
        self.c_rational_recon_int(a, m, &n, &d)
        return (n,d)


# The long long versions are next.
cdef class arith_llong:

    cdef public long long abs_longlong(self, long long x) except -1:
        if x < 0:
            return -x
        return x

    cdef public long long sign_longlong(self, long long n) except -2:
        if n < 0:
            return -1
        return 1

    cdef public long long c_gcd_longlong(self, long long a, long long b) except -1:
        cdef long long c
        if a==0:
            return self.abs_longlong(b)
        if b==0:
            return self.abs_longlong(a)
        if a<0: a=-a
        if b<0: b=-b
        while(b):
            c = a % b
            a = b
            b = c
        return a


    def gcd_longlong(self, long long a, long long b):
        return self.c_gcd_longlong(a,b)


    cdef public long long c_xgcd_longlong(self, long long a, long long b,
                                          long long *ss,
                                          long long *tt) except -1:
        cdef long long psign, qsign, p, q, r, s, c, quot, new_r, new_s


        if a == 0:
            ss[0] = 0
            tt[0] = self.sign_longlong(b)
            return self.abs_longlong(b)

        if b == 0:
            ss[0] = self.sign_longlong(a)
            tt[0] = 0
            return self.abs_longlong(a)

        psign = 1; qsign = 1

        if a<0: a = -a; psign = -1
        if b<0: b = -b; qsign = -1

        p = 1; q = 0; r = 0; s = 1
        while (b):
            c = a % b; quot = a/b
            a = b; b = c
            new_r = p - quot*r
            new_s = q - quot*s
            p = r; q = s
            r = new_r; s = new_s

        ss[0] = p*psign
        tt[0] = q*qsign


        return a

    cdef public long long c_inverse_mod_longlong(self, long long a, long long m) except -1:
        cdef long long g, s, t
        g = self.c_xgcd_longlong(a,m, &s, &t)
        if g != 1:
            raise ArithmeticError("The inverse of %s modulo %s is not defined."%(a,m))
        s = s % m
        if s < 0:
            s = s + m
        return s

    def inverse_mod_longlong(self, long long a, long long m):
        return self.c_inverse_mod_longlong(a, m)

    cdef long long c_rational_recon_longlong(self, long long a, long long m,
                                             long long *n, long long *d) except -1:
        cdef long long u, v, u0, u1, u2, v0, v1, v2, q, t0, t1, t2, x, y
        cdef float bnd

        if m > 2147483647:
            raise OverflowError, "The modulus m(=%s) must be at most 2147483647"%m
            return -1

        a = a % m

        if a==0 or m == 0:
            n[0] = 0
            d[0] = 1
            return 0

        if m<0: m = -m
        if a<0: a = m - a
        if a==1:
            n[0] = 1
            d[0] = 1
            return 0

        u = m
        v = a
        bnd = sqrt(m/2.0)
        u0=1; u1=0; u2=u
        v0=0; v1=1; v2=v
        while self.abs_longlong(v2) > bnd:
            q = u2/v2   # floor is implicit
            t0=u0-q*v0; t1=u1-q*v1; t2=u2-q*v2
            u0=v0; u1=v1; u2=v2
            v0=t0; v1=t1; v2=t2;

        x = self.abs_longlong(v1); y = v2
        if v1<0:  y = -1*y
        if x<=bnd and self.gcd_longlong(x,y)==1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_longlong(self, long long a, long long m):
        """
        Rational reconstruction of a modulo m.
        """
        cdef long long n, d
        self.c_rational_recon_longlong(a, m, &n, &d)
        return (n,d)




