# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
ntl_lzz_pX.pyx

Wraps NTL's zz_pX type for SAGE

AUTHORS:
   - Craig Citro
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

include 'misc.pxi'
include 'decl.pxi'
from sage.libs.gmp.mpz cimport *

from cpython.object cimport Py_EQ, Py_NE
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer cimport Integer

from sage.rings.finite_rings.integer_mod cimport IntegerMod_gmp, IntegerMod_int, IntegerMod_int64

from sage.libs.ntl.ntl_lzz_pContext import ntl_zz_pContext
from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

from sage.libs.ntl.ntl_lzz_p import ntl_zz_p
from sage.libs.ntl.ntl_lzz_p cimport ntl_zz_p

ZZ_sage = IntegerRing()

##############################################################################
#
# zz_pX  -- polynomials over the integers modulo p, p small
#
##############################################################################

cdef class ntl_zz_pX():
    r"""
    The class \class{zz_pX} implements polynomial arithmetic modulo `p`,
    for p smaller than a machine word.

    Polynomial arithmetic is implemented using the FFT, combined with
    the Chinese Remainder Theorem.  A more detailed description of the
    techniques used here can be found in [Shoup, J. Symbolic
    Comp. 20:363-397, 1995].

    Small degree polynomials are multiplied either with classical
    or Karatsuba algorithms.
    """
    # See ntl_zz_pX.pxd for definition of data members
    def __init__(self, ls=[], modulus=None):
        """
        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,5,-9],20)
            sage: f
            [1, 2, 5, 11]
            sage: g = ntl.zz_pX([0,0,0],20); g
            []
            sage: g[10]=5
            sage: g
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5]
            sage: g[10]
            5
            sage: f = ntl.zz_pX([10^30+1, 10^50+1], 100); f
            [1, 1]
        """
        if modulus is None:
            raise ValueError("You must specify a modulus.")

        cdef long n
        cdef Py_ssize_t i
        cdef long temp

        if isinstance(modulus, Integer):
            p_sage = modulus
        else:
            p_sage = Integer(self.c.p)

        #self.c.restore_c() ## We did this in __new__

        n = len(ls)
        if (n == 0):
            ## the 0 polynomial is just the empty list;
            ## so in this case, we're done.
            return

        self.x.SetMaxLength(n+1)

        for i from 0 <= i < n:
            a = ls[i]

            if isinstance(a, IntegerMod_int):
                if (self.c.p == (<IntegerMod_int>a)._modulus.int32): ## this is slow
                    zz_pX_SetCoeff_long(self.x, i, (<IntegerMod_int>a).ivalue)
                else:
                    raise ValueError("Mismatched modulus for converting to zz_pX.")
            elif isinstance(a, IntegerMod_int64):
                if (self.c.p == (<IntegerMod_int64>a)._modulus.int64): ## this is slow
                    zz_pX_SetCoeff_long(self.x, i, (<IntegerMod_int64>a).ivalue)
                else:
                    raise ValueError("Mismatched modulus for converting to zz_pX.")
            elif isinstance(a, IntegerMod_gmp):
                if (p_sage == (<IntegerMod_gmp>a)._modulus.sageInteger): ## this is slow
                    zz_pX_SetCoeff_long(self.x, i, mpz_get_si((<IntegerMod_gmp>a).value))
                else:
                    raise ValueError("Mismatched modulus for converting to zz_pX.")
            elif isinstance(a, Integer):
                zz_pX_SetCoeff_long(self.x, i, mpz_fdiv_ui((<Integer>a).value, self.c.p))
            elif isinstance(a, int):
                ## we're lucky that python int is no larger than long
                temp = a
                zz_pX_SetCoeff_long(self.x, i, temp%self.c.p)
            else:
                a = Integer(a)
                zz_pX_SetCoeff_long(self.x, i, mpz_fdiv_ui((<Integer>a).value, self.c.p))

        return

    def __cinit__(self, v=None, modulus=None):
        #################### WARNING ###################
        ## Before creating a zz_pX, you must create a ##
        ## zz_pContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing a zz_pX              ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = ntl_zz_pX.__new__(ntl_zz_pX) without
        ## first restoring a zz_pContext, which could ##
        ## have unfortunate consequences.  See _new  ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None:
            return
        if isinstance(modulus, ntl_zz_pContext_class):
            self.c = <ntl_zz_pContext_class>modulus
        elif isinstance(modulus, Integer):
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)
        elif isinstance(modulus, int):
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)
        else:
            try:
                modulus = int(modulus)
            except Exception:
                raise ValueError("%s is not a valid modulus." % modulus)
            self.c = <ntl_zz_pContext_class>ntl_zz_pContext(modulus)

        ## now that we've determined the modulus, set that modulus.
        self.c.restore_c()

    def __reduce__(self):
        """
        TESTS::

            sage: f = ntl.zz_pX([10,10^30+1], 20)
            sage: f == loads(dumps(f))
            True
        """
        return make_zz_pX, (self.list(), self.c)

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: f = ntl.zz_pX([3,5], 17)
            sage: f.__repr__()
            '[3, 5]'
        """
        return str(self.list())

    def __getitem__(self, long i):
        """
        Return the i-th coefficient of f.

        EXAMPLES::

            sage: f = ntl.zz_pX(range(7), 71)
            sage: f[3] ## indirect doctest
            3

            sage: f[-5]
            0

            sage: f[27]
            0
        """
        cdef ntl_zz_p y
        y = ntl_zz_p.__new__(ntl_zz_p)
        y.c = self.c
        self.c.restore_c()
        y.x = zz_pX_GetCoeff(self.x, i)
        return y

    def __setitem__(self, long i, val):
        """
        Set the i-th coefficient of ``self`` to val. If
        i is out of range, raise an exception.

        EXAMPLES::

            sage: f = ntl.zz_pX([], 7)
            sage: f[3] = 2 ; f
            [0, 0, 0, 2]
            sage: f[-1] = 5
            Traceback (most recent call last):
            ...
            ValueError: index (=-1) is out of range
        """
        if (i < 0):
            raise ValueError("index (=%s) is out of range" % i)
        self.c.restore_c()
        zz_pX_SetCoeff_long(self.x, i, val)
        return

    cdef ntl_zz_pX _new(self):
        """
        Quick and dirty method for creating a new object with the
        same zz_pContext as ``self``.

        EXAMPLES::

            sage: f = ntl.zz_pX([1], 20)
            sage: f.square() ## indirect doctest
            [1]
        """
        cdef ntl_zz_pX y
        y = ntl_zz_pX.__new__(ntl_zz_pX)
        y.c = self.c
        return y

    def __add__(ntl_zz_pX self, other):
        """
        Return ``self + other``.

        EXAMPLES::

            sage: ntl.zz_pX(range(5),20) + ntl.zz_pX(range(6),20) ## indirect doctest
            [0, 2, 4, 6, 8, 5]
            sage: ntl.zz_pX(range(5),20) + ntl.zz_pX(range(6),50)
            Traceback (most recent call last):
            ...
            ValueError: arithmetic operands must have the same modulus.
        """
        cdef ntl_zz_pX y
        if not isinstance(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        y = self._new()
        self.c.restore_c()
        zz_pX_add(y.x, self.x, (<ntl_zz_pX>other).x)
        return y

    def __sub__(ntl_zz_pX self, other):
        """
        Return ``self - other``.

        EXAMPLES::

            sage: ntl.zz_pX(range(5),32) - ntl.zz_pX(range(6),32)
            [0, 0, 0, 0, 0, 27]
            sage: ntl.zz_pX(range(5),20) - ntl.zz_pX(range(6),50) ## indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: arithmetic operands must have the same modulus.
        """
        cdef ntl_zz_pX y
        if not isinstance(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        self.c.restore_c()
        y = self._new()
        zz_pX_sub(y.x, self.x, (<ntl_zz_pX>other).x)
        return y

    def __mul__(ntl_zz_pX self, other):
        """
        EXAMPLES::

            sage: ntl.zz_pX(range(5),20) * ntl.zz_pX(range(6),20) ## indirect doctest
            [0, 0, 1, 4, 10, 0, 10, 14, 11]
            sage: ntl.zz_pX(range(5),20) * ntl.zz_pX(range(6),50)
            Traceback (most recent call last):
            ...
            ValueError: arithmetic operands must have the same modulus.
        """
        cdef ntl_zz_pX y
        if not isinstance(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        self.c.restore_c()
        y = self._new()
        sig_on()
        zz_pX_mul(y.x, self.x, (<ntl_zz_pX>other).x)
        sig_off()
        return y

    def __truediv__(ntl_zz_pX self, other):
        """
        Compute quotient ``self / other``, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3],17) * ntl.zz_pX([4,5],17)**2
            sage: g = ntl.zz_pX([4,5],17)
            sage: f/g ## indirect doctest
            [4, 13, 5, 15]
            sage: ntl.zz_pX([1,2,3],17) * ntl.zz_pX([4,5],17)
            [4, 13, 5, 15]

            sage: f = ntl.zz_pX(range(10),17); g = ntl.zz_pX([-1,0,1],17)
            sage: f/g
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) is not divisible by other (=[16, 0, 1])
            sage: ntl.zz_pX(range(5),20) / ntl.zz_pX(range(6),50)
            Traceback (most recent call last):
            ...
            ValueError: arithmetic operands must have the same modulus.
        """
        cdef long divisible
        cdef ntl_zz_pX q
        if not isinstance(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        self.c.restore_c()
        q = self._new()
        sig_on()
        divisible = zz_pX_divide(q.x, self.x, (<ntl_zz_pX>other).x)
        sig_off()
        if not divisible:
            raise ArithmeticError("self (=%s) is not divisible by other (=%s)" % (self, other))
        return q

    def __mod__(ntl_zz_pX self, other):
        """
        Given polynomials a, b in ZZ[X], there exist polynomials q, r
        in QQ[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns q if q lies in ZZ[X], and otherwise raises an
        Exception.

        EXAMPLES::

            sage: f = ntl.zz_pX([2,4,6],17); g = ntl.zz_pX([2],17)
            sage: f % g   ## indirect doctest
            []

            sage: f = ntl.zz_pX(range(10),17); g = ntl.zz_pX([-1,0,1],17)
            sage: f % g
            [3, 8]
        """
        cdef ntl_zz_pX y
        if not isinstance(other, ntl_zz_pX):
            other = ntl_zz_pX(other, modulus=self.c)
        elif self.c is not (<ntl_zz_pX>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        self.c.restore_c()
        y = self._new()
        sig_on()
        zz_pX_mod(y.x, self.x, (<ntl_zz_pX>other).x)
        sig_off()
        return y

    def __pow__(ntl_zz_pX self, long n, ignored):
        """
        Return the `n`-th nonnegative power of ``self``.

        EXAMPLES::

            sage: g = ntl.zz_pX([-1,0,1],20)
            sage: g**10 ## indirect doctest
            [1, 0, 10, 0, 5, 0, 0, 0, 10, 0, 8, 0, 10, 0, 0, 0, 5, 0, 10, 0, 1]
        """
        if n < 0:
            raise ValueError("Only positive exponents allowed.")
        cdef ntl_zz_pX y = self._new()
        self.c.restore_c()
        sig_on()
        zz_pX_power(y.x, self.x, n)
        sig_off()
        return y

    def quo_rem(ntl_zz_pX self, ntl_zz_pX right):
        """
        Return the quotient and remainder when ``self`` is divided by ``right``.

        Specifically, this returns `r`, `q` such that ``self = q * right + r``.

        EXAMPLES::

            sage: f = ntl.zz_pX(range(7), 19)
            sage: g = ntl.zz_pX([2,4,6], 19)
            sage: f // g
            [1, 1, 15, 16, 1]
            sage: f % g
            [17, 14]
            sage: f.quo_rem(g)
            ([1, 1, 15, 16, 1], [17, 14])
            sage: (f // g) * g + f % g
            [0, 1, 2, 3, 4, 5, 6]
        """
        cdef ntl_zz_pX q = self._new()
        cdef ntl_zz_pX r = self._new()
        self.c.restore_c()
        sig_on()
        zz_pX_divrem(q.x, r.x, self.x, right.x)
        sig_off()
        return q, r

    def __floordiv__(ntl_zz_pX self, ntl_zz_pX right):
        """
        Return the whole part of ``self / right``.

        EXAMPLES::

            sage: f = ntl.zz_pX(range(10), 19); g = ntl.zz_pX([1]*5, 19)
            sage: f // g ## indirect doctest
            [8, 18, 18, 18, 18, 9]
        """
        cdef ntl_zz_pX q = self._new()
        self.c.restore_c()
        sig_on()
        zz_pX_div(q.x, self.x, right.x)
        sig_off()
        return q

    def __lshift__(ntl_zz_pX self, long n):
        """
        Shift this polynomial to the left, which is multiplication by `x^n`.

        EXAMPLES::

            sage: f = ntl.zz_pX([2,4,6], 17)
            sage: f << 2 ## indirect doctest
            [0, 0, 2, 4, 6]
        """
        cdef ntl_zz_pX r = self._new()
        self.c.restore_c()
        zz_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(ntl_zz_pX self, long n):
        """
        Shift this polynomial to the right, which is division by `x^n` (and truncation).

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3], 17)
            sage: f >> 2 ## indirect doctest
            [3]
        """
        cdef ntl_zz_pX r = self._new()
        self.c.restore_c()
        zz_pX_RightShift(r.x, self.x, n)
        return r

    def diff(self):
        """
        The formal derivative of ``self``.

        EXAMPLES::

            sage: f = ntl.zz_pX(range(10), 17)
            sage: f.diff()
            [1, 4, 9, 16, 8, 2, 15, 13, 13]
        """
        cdef ntl_zz_pX r = self._new()
        self.c.restore_c()
        zz_pX_diff(r.x, self.x)
        return r

    def reverse(self):
        """
        Return ``self`` with coefficients reversed, i.e. ``x^n self(x^{-n})``.

        EXAMPLES::

            sage: f = ntl.zz_pX([2,4,6], 17)
            sage: f.reverse()
            [6, 4, 2]
        """
        cdef ntl_zz_pX r = self._new()
        self.c.restore_c()
        zz_pX_reverse(r.x, self.x)
        return r

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: f = ntl.zz_pX([2,0,0,1],20)
            sage: -f
            [18, 0, 0, 19]
        """
        cdef ntl_zz_pX y
        y = self._new()
        self.c.restore_c()
        sig_on()
        zz_pX_negate(y.x, self.x)
        sig_off()
        return y

    def __richcmp__(ntl_zz_pX self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3],20)
            sage: g = ntl.zz_pX([1,2,3,0],20)
            sage: f == g
            True
            sage: g = ntl.zz_pX([0,1,2,3],20)
            sage: f == g
            False
            sage: f != [0]
            True
        """
        self.c.restore_c()

        if op != Py_EQ and op != Py_NE:
            raise TypeError("polynomials are not ordered")

        cdef ntl_zz_pX b
        try:
            b = <ntl_zz_pX?>other
        except TypeError:
            b = ntl_zz_pX(other, self.c)

        return (op == Py_EQ) == (self.x == b.x)

    def list(self):
        """
        Return list of entries as a list of python ints.

        EXAMPLES::

            sage: f = ntl.zz_pX([23, 5,0,1], 10)
            sage: f.list()
            [3, 5, 0, 1]
            sage: type(f.list()[0])
            <... 'int'>
        """
        cdef long i
        self.c.restore_c()
        return [ zz_p_rep(zz_pX_GetCoeff(self.x, i)) for i from 0 <= i <= zz_pX_deg(self.x) ]

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES::

            sage: f = ntl.zz_pX([5,0,1],50)
            sage: f.degree()
            2
            sage: f = ntl.zz_pX(range(100),50)
            sage: f.degree()
            99
            sage: f = ntl.zz_pX([],10)
            sage: f.degree()
            -1
            sage: f = ntl.zz_pX([1],77)
            sage: f.degree()
            0
        """
        self.c.restore_c()
        return zz_pX_deg(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES::

            sage: f = ntl.zz_pX([3,6,9],19)
            sage: f.leading_coefficient()
            9
            sage: f = ntl.zz_pX([],21)
            sage: f.leading_coefficient()
            0
        """
        self.c.restore_c()
        return zz_p_rep(zz_pX_LeadCoeff(self.x))

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES::

            sage: f = ntl.zz_pX([3,6,9],127)
            sage: f.constant_term()
            3
            sage: f = ntl.zz_pX([], 12223)
            sage: f.constant_term()
            0
        """
        self.c.restore_c()
        return zz_p_rep(zz_pX_ConstTerm(self.x))

    def square(self):
        """
        Return f*f.

        EXAMPLES::

            sage: f = ntl.zz_pX([-1,0,1],17)
            sage: f*f
            [1, 0, 15, 0, 1]
        """
        cdef ntl_zz_pX y = self._new()
        self.c.restore_c()
        sig_on()
        zz_pX_sqr(y.x, self.x)
        sig_off()
        return y

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3,4,5],70)
            sage: f.truncate(3)
            [1, 2, 3]
            sage: f.truncate(8)
            [1, 2, 3, 4, 5]
            sage: f.truncate(1)
            [1]
            sage: f.truncate(0)
            []
            sage: f.truncate(-1)
            []
            sage: f.truncate(-5)
            []
        """
        cdef ntl_zz_pX y = self._new()
        self.c.restore_c()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            sig_on()
            zz_pX_trunc(y.x, self.x, m)
            sig_off()
        return y

    def multiply_and_truncate(self, ntl_zz_pX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3,4,5],20)
            sage: g = ntl.zz_pX([10],20)
            sage: f.multiply_and_truncate(g, 2)
            [10]
            sage: g.multiply_and_truncate(f, 2)
            [10]
        """
        cdef ntl_zz_pX y = self._new()
        self.c.restore_c()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            sig_on()
            zz_pX_MulTrunc(y.x, self.x, other.x, m)
            sig_off()
        return y

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3,4,5],20)
            sage: f.square_and_truncate(4)
            [1, 4, 10]
            sage: (f*f).truncate(4)
            [1, 4, 10]
        """
        cdef ntl_zz_pX y = self._new()
        self.c.restore_c()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            sig_on()
            zz_pX_SqrTrunc(y.x, self.x, m)
            sig_off()
        return y

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of ``self`` modulo `x^m`.
        The constant term of ``self`` must be 1 or -1.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3,4,5,6,7],20)
            sage: f.invert_and_truncate(20)
            [1, 18, 1, 0, 0, 0, 0, 8, 17, 2, 13, 0, 0, 0, 4, 0, 17, 10, 9]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 4, 3]
        """
        if m < 0:
            raise ArithmeticError("m (=%s) must be positive" % m)
        n = self.constant_term()
        if n != 1 and n != -1:
            raise ArithmeticError("The constant term of self must be 1 or -1.")

        cdef ntl_zz_pX y = self._new()

        self.c.restore_c()
        if m <= 0:
            y.x = zz_pX_zero()
        else:
            sig_on()
            zz_pX_InvTrunc(y.x, self.x, m)
            sig_off()
        return y

    def is_zero(self):
        """
        Return ``True`` exactly if this polynomial is 0.

        EXAMPLES::

            sage: f = ntl.zz_pX([0,0,0,20],5)
            sage: f.is_zero()
            True
            sage: f = ntl.zz_pX([0,0,1],30)
            sage: f
            [0, 0, 1]
            sage: f.is_zero()
            False
        """
        self.c.restore_c()
        return zz_pX_IsZero(self.x)

    def is_one(self):
        """
        Return ``True`` exactly if this polynomial is 1.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,1],101)
            sage: f.is_one()
            False
            sage: f = ntl.zz_pX([1],2)
            sage: f.is_one()
            True
        """
        self.c.restore_c()
        return zz_pX_IsOne(self.x)

    def is_monic(self):
        """
        Return ``True`` exactly if this polynomial is monic.

        EXAMPLES::

            sage: f = ntl.zz_pX([2,0,0,1],17)
            sage: f.is_monic()
            True
            sage: g = f.reverse()
            sage: g.is_monic()
            False
            sage: g
            [1, 0, 0, 2]
            sage: f = ntl.zz_pX([1,2,0,3,0,2],717)
            sage: f.is_monic()
            False
        """
        self.c.restore_c()
        if zz_pX_IsZero(self.x):
            return False
        return zz_p_rep(zz_pX_LeadCoeff(self.x)) == 1

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES::

            sage: f = ntl.zz_pX([],177)
            sage: f.set_x()
            sage: f
            [0, 1]
            sage: g = ntl.zz_pX([0,1],177)
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory:
            sage: f is g
            False
        """
        self.c.restore_c()
        zz_pX_SetX(self.x)

    def is_x(self):
        """
        ``True`` if this is the polynomial "x".

        EXAMPLES::

            sage: f = ntl.zz_pX([],100)
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = ntl.zz_pX([0,1],383)
            sage: f.is_x()
            True
            sage: f = ntl.zz_pX([1],38)
            sage: f.is_x()
            False
        """
        self.c.restore_c()
        return zz_pX_IsX(self.x)

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3],17)
            sage: f
            [1, 2, 3]
            sage: f.clear()
            sage: f
            []
        """
        self.c.restore_c()
        zz_pX_clear(self.x)

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES::

            sage: f = ntl.zz_pX([1,2,3],17)
            sage: f.preallocate_space(20)
            sage: f
            [1, 2, 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 5]
        """
        self.c.restore_c()
        sig_on()
        self.x.SetMaxLength(n)
        sig_off()
        return


def make_zz_pX(L, context):
    """
    For unpickling.

    TESTS::

        sage: f = ntl.zz_pX(range(16), 12)
        sage: loads(dumps(f)) == f
        True
    """
    return ntl_zz_pX(L, context)
