# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
ntl_lzz_p.pyx

Wraps NTL's zz_p type for Sage

NOTE: This file is essentially useless. While we provide
this wrapper for consistency, this should never be used in
*production* code, i.e. anything intended to be fast. The
reason for this is simple: this is a wrapper for a Python
interface to the zz_p type, which is just a long! Any speed
gains you get from working with longs will be TOTALLY
destroyed by the overhead of having a wrapper.

AUTHORS:

- Craig Citro
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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
from sage.arith.power cimport generic_power_pos

ZZ_sage = IntegerRing()

##############################################################################
#
# zz_pX  -- polynomials over the integers modulo p, p small
#
##############################################################################

cdef class ntl_zz_p():
    r"""
    The class \class{zz_p} implements arithmetic modulo `p`,
    for p smaller than a machine word.

    NOTE: This type is provided mostly for completeness, and
    shouldn't be used in any production code.
    """
    # See ntl_zz_p.pxd for definition of data members
    def __init__(self, a=0, modulus=None):
        """
        EXAMPLES::

            sage: f = ntl.zz_p(5,7)
            sage: f
            5
            sage: g = ntl.zz_p(int(-5),7)
            sage: g
            2
        """
        if modulus is None:
            raise ValueError("You must specify a modulus.")

        if isinstance(modulus, Integer):
            p_sage = modulus
        else:
            p_sage = Integer(self.c.p)

        #self.c.restore_c()   ## This was done in __new__

        if isinstance(a, IntegerMod_int):
            if (self.c.p == (<IntegerMod_int>a)._modulus.int32): ## this is slow
                self.x = (<IntegerMod_int>a).ivalue
            else:
                raise ValueError("Mismatched modulus for converting to zz_p.")

        elif isinstance(a, IntegerMod_int64):
            if (self.c.p == (<IntegerMod_int64>a)._modulus.int64): ## this is slow
                self.x = (<IntegerMod_int64>a).ivalue
            else:
                raise ValueError("Mismatched modulus for converting to zz_p.")

        elif isinstance(a, IntegerMod_gmp):
            if (p_sage == (<IntegerMod_gmp>a)._modulus.sageInteger): ## this is slow
                self.x = mpz_get_si((<IntegerMod_gmp>a).value)
            else:
                raise ValueError("Mismatched modulus for converting to zz_p.")

        elif isinstance(a, Integer):
            self.x = mpz_get_si((<Integer>a).value)%self.c.p

        elif isinstance(a, int):
            ## we're lucky that python int is no larger than long
            self.x = (<long>a)%self.c.p
        else:
            a = Integer(a)
            self.x = mpz_get_si((<Integer>a).value)%self.c.p

        return

    def __cinit__(self, v=None, modulus=None):
        #################### WARNING ###################
        ## Before creating a zz_p, you must create a  ##
        ## zz_pContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing a zz_p               ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = ntl_zz_p.__new__(ntl_zz_p) without
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

    cdef ntl_zz_p _new(self):
        """
        Quick and dirty zz_p object creation.

        EXAMPLES::

            sage: x = ntl.zz_p(23,75)
            sage: y = x*x ## indirect doctest
        """
        cdef ntl_zz_p y
        self.c.restore_c()
        y = ntl_zz_p.__new__(ntl_zz_p)
        y.c = self.c
        return y

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: f = ntl.zz_p(16,244)
            sage: loads(dumps(f)) == f
            True
        """
        return make_zz_p, (zz_p_rep(self.x), self.c)

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ntl.zz_p(3,79).__repr__()
            '3'
        """
        return repr(Integer(zz_p_rep(self.x)))

    def __add__(ntl_zz_p self, other):
        """
        EXAMPLES::

            sage: ntl.zz_p(5,23) + ntl.zz_p(6,23)
            11
        """
        cdef ntl_zz_p y
        if not isinstance(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        self.c.restore_c()
        y = self._new()
        zz_p_add(y.x, self.x, (<ntl_zz_p>other).x)
        return y

    def __sub__(ntl_zz_p self, other):
        """
        EXAMPLES::

            sage: ntl.zz_p(5,23) - ntl.zz_p(6,23)
            22
        """
        cdef ntl_zz_p y
        if not isinstance(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        self.c.restore_c()
        y = self._new()
        zz_p_sub(y.x, self.x, (<ntl_zz_p>other).x)
        return y

    def __mul__(ntl_zz_p self, other):
        """
        EXAMPLES::

            sage: ntl.zz_p(5,23) * ntl.zz_p(6,23)
            7
        """
        cdef ntl_zz_p y
        if not isinstance(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        y = self._new()
        self.c.restore_c()
        zz_p_mul(y.x, self.x, (<ntl_zz_p>other).x)
        return y

    def __truediv__(ntl_zz_p self, other):
        """
        EXAMPLES::

            sage: ntl.zz_p(5,23) / ntl.zz_p(2,23)
            14
        """
        cdef ntl_zz_p q
        if not isinstance(other, ntl_zz_p):
            other = ntl_zz_p(other, modulus=self.c)
        elif self.c is not (<ntl_zz_p>other).c:
            raise ValueError("arithmetic operands must have the same modulus.")
        q = self._new()
        self.c.restore_c()
        sig_on()
        zz_p_div(q.x, self.x, (<ntl_zz_p>other).x)
        sig_off()
        return q

    def __pow__(ntl_zz_p self, long n, ignored):
        """
        Return the `n`-th nonnegative power of ``self``.

        EXAMPLES::

            sage: g = ntl.zz_p(5, 13)
            sage: g ^ 10
            12
            sage: g ^ (-1)
            8
            sage: g ^ (-5)
            8
            sage: g ^ 0
            1
            sage: z = ntl.zz_p(0, 13)
            sage: z ^ 0
            1
            sage: z ^ 1
            0
            sage: z ^ (-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse does not exist
        """
        self.c.restore_c()
        if n == 0:
            return ntl_zz_p(1, self.c)

        if self.is_zero():
            if n > 0:
                return self
            raise ZeroDivisionError("inverse does not exist")

        cdef ntl_zz_p y
        if n > 0:
            return generic_power_pos(self, <unsigned long>n)
        else:
            y = ntl_zz_p.__new__(ntl_zz_p)
            y.c = self.c
            zz_p_inv(y.x, self.x)
            return generic_power_pos(y, -<unsigned long>n)

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: f = ntl.zz_p(5,234)
            sage: -f ## indirect doctest
            229
        """
        cdef ntl_zz_p y
        y = self._new()
        self.c.restore_c()
        zz_p_negate(y.x, self.x)
        return y

    def __richcmp__(ntl_zz_p self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: f = ntl.zz_p(3,20)
            sage: g = ntl.zz_p(2,20)
            sage: h = ntl.zz_p(3,60)
            sage: f == g
            False
            sage: f == f
            True
            sage: f == h
            True
            sage: f == 3
            True
        """
        self.c.restore_c()

        if op != Py_EQ and op != Py_NE:
            raise TypeError("integers mod p are not ordered")

        cdef ntl_zz_p b
        try:
            b = <ntl_zz_p?>other
        except TypeError:
            b = ntl_zz_p(other, self.c)

        return (op == Py_EQ) == (self.x == b.x)

    def __int__(self):
        """
        Return ``self`` as an int.

        EXAMPLES::

            sage: ntl.zz_p(3,next_prime(100)).__int__()
            3
            sage: int(ntl.zz_p(3,next_prime(100)))
            3
            sage: type(int(ntl.zz_p(3,next_prime(100))))
            <... 'int'>
        """
        return zz_p_rep(self.x)

    def square(self):
        """
        Return f*f.

        EXAMPLES::

            sage: f = ntl.zz_p(15,23)
            sage: f*f
            18
        """
        cdef ntl_zz_p y
        y = self._new()
        self.c.restore_c()
        zz_p_sqr(y.x, self.x)
        return y

    def is_zero(self):
        """
        Return ``True`` exactly if this element is 0.

        EXAMPLES::

            sage: f = ntl.zz_p(0,20)
            sage: f.is_zero()
            True
            sage: f = ntl.zz_p(1,20)
            sage: f.is_zero()
            False
        """
        self.c.restore_c()
        return zz_p_rep(self.x) == 0

    def is_one(self):
        """
        Return ``True`` exactly if this element is 1.

        EXAMPLES::

            sage: f = ntl.zz_p(1,11)
            sage: f.is_one()
            True
            sage: f = ntl.zz_p(5,11)
            sage: f.is_one()
            False
        """
        self.c.restore_c()
        return zz_p_rep(self.x) == 1

    def clear(self):
        """
        Reset this element to 0.

        EXAMPLES::

            sage: x = ntl.zz_p(5,102) ; x
            5
            sage: x.clear() ; x
            0
        """
        self.c.restore_c()
        zz_p_clear(self.x)


def make_zz_p(val, context):
    """
    For unpickling.

    TESTS::

        sage: f = ntl.zz_p(1, 12)
        sage: loads(dumps(f)) == f
        True
    """
    return ntl_zz_p(val, context)
