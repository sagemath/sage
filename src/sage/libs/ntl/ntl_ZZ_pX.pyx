# sage_setup: distribution = sagemath-ntl
# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from sage.ext.cplusplus cimport ccrepr, ccreadstr

include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from sage.cpython.string cimport char_to_str
from sage.rings.integer cimport Integer
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZ_p cimport ntl_ZZ_p
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext
from sage.libs.ntl.ntl_ZZ import unpickle_class_args
from sage.misc.randstate cimport current_randstate
from sage.libs.gmp.mpz cimport *


cdef inline make_ZZ_p(ZZ_p_c* x, ntl_ZZ_pContext_class ctx):
    cdef ntl_ZZ_p y
    sig_off()
    y = ntl_ZZ_p(modulus = ctx)
    y.x = x[0]
    del x
    return y


cdef make_ZZ_pX(ZZ_pX_c* x, ntl_ZZ_pContext_class ctx):
    cdef ntl_ZZ_pX y
    y = <ntl_ZZ_pX>ntl_ZZ_pX.__new__(ntl_ZZ_pX)
    y.c = ctx
    y.x = x[0]
    del x
    sig_off()
    return y


##############################################################################
#
# ZZ_pX  -- polynomials over the integers modulo p
#
##############################################################################

cdef class ntl_ZZ_pX():
    r"""
    The class \class{ZZ_pX} implements polynomial arithmetic modulo `p`.

    Polynomial arithmetic is implemented using the FFT, combined with
    the Chinese Remainder Theorem.  A more detailed description of the
    techniques used here can be found in [Shoup, J. Symbolic
    Comp. 20:363-397, 1995].

    Small degree polynomials are multiplied either with classical
    or Karatsuba algorithms.
    """

    # See ntl_ZZ_pX.pxd for definition of data members
    def __init__(self, v=None, modulus=None):
        """
        EXAMPLES::

            sage: c = ntl.ZZ_pContext(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,5,-9], c)
            sage: f
            [1 2 5 11]
            sage: g = ntl.ZZ_pX([0,0,0], c); g
            []
            sage: g[10]=5
            sage: g
            [0 0 0 0 0 0 0 0 0 0 5]
            sage: g[10]
            5
        """
        if modulus is None:
            raise ValueError("You must specify a modulus when creating a ZZ_pX.")

        # self.c._assert_is_current_modulus()  # the context was restored in __new__
        cdef ntl_ZZ_p cc
        cdef Py_ssize_t i

        if isinstance(v, ntl_ZZ_pX) and (<ntl_ZZ_pX>v).c is self.c:
            self.x = (<ntl_ZZ_pX>v).x
        elif isinstance(v, (list, tuple, xrange)):
            for i, x in enumerate(v):
                if not isinstance(x, ntl_ZZ_p):
                    cc = ntl_ZZ_p(x, self.c)
                    self.c.restore_c()
                else:
                    cc = x
                ZZ_pX_SetCoeff(self.x, i, cc.x)
        elif v is not None:
            s = str(v).replace(',', ' ').replace('L', '')  # can change the modulus; Issue #25790
            self.c.restore_c()
            ccreadstr(self.x, s)

    def __cinit__(self, v=None, modulus=None):
        #################### WARNING ###################
        ## Before creating a ZZ_pX, you must create a ##
        ## ZZ_pContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing an ntl_ZZ_p          ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = ntl_ZZ_p.__new__(ntl_ZZ_p) without
        ## first restoring a ZZ_pContext, which could ##
        ## have unfortunate consequences.  See _new   ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None:
            return
        if isinstance(modulus, ntl_ZZ_pContext_class):
            self.c = <ntl_ZZ_pContext_class>modulus
            self.c.restore_c()
        elif modulus is not None:
            self.c = <ntl_ZZ_pContext_class>ntl_ZZ_pContext(ntl_ZZ(modulus))
            self.c.restore_c()

    cdef ntl_ZZ_pX _new(self):
        cdef ntl_ZZ_pX r
        self.c.restore_c()
        r = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        r.c = self.c
        return r

    def __reduce__(self):
        """
        TESTS::

            sage: f = ntl.ZZ_pX([1,2,3,7], 5); f
            [1 2 3 2]
            sage: loads(dumps(f)) == f
            True
        """
        return unpickle_class_args, (ntl_ZZ_pX, (self.list(), self.get_modulus_context()))

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: x = ntl.ZZ_pX([1,0,8],5)
            sage: x
            [1 0 3]
            sage: x.__repr__()
            '[1 0 3]'
        """
        self.c.restore_c()
        return ccrepr(self.x)

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: x = ntl.ZZ_pX([0,5,-3],11)
            sage: y = copy(x)
            sage: x == y
            True
            sage: x is y
            False
            sage: x[0] = 4; y
            [0 5 8]
        """
        cdef ntl_ZZ_pX r = self._new()
        #self.c.restore_c() # restored in _new()
        r.x = self.x
        return r

    def get_modulus_context(self):
        """
        Return the modulus for ``self``.

        EXAMPLES::

            sage: x = ntl.ZZ_pX([0,5,3],17)
            sage: c = x.get_modulus_context()
            sage: y = ntl.ZZ_pX([5],c)
            sage: x+y
            [5 5 3]
        """
        return self.c

    def __setitem__(self, long i, a):
        r"""
        EXAMPLES::

            sage: c = ntl.ZZ_pContext(23)
            sage: x = ntl.ZZ_pX([2, 3, 4], c)
            sage: x[1] = 5
            sage: x
            [2 5 4]
        """
        if i < 0:
            raise IndexError("index (i=%s) must be >= 0" % i)
        cdef ntl_ZZ_p _a
        _a = ntl_ZZ_p(a,self.c)
        self.c.restore_c()
        ZZ_pX_SetCoeff(self.x, i, _a.x)

    cdef void setitem_from_int(ntl_ZZ_pX self, long i, int value) noexcept:
        r"""
        Set `i`-th coefficient to value.

        AUTHOR: David Harvey (2006-08-05)
        """
        self.c.restore_c()
        ZZ_pX_SetCoeff_long(self.x, i, value)

    def _setitem_from_int_doctest(self, i, value):
        r"""
        This method exists solely for automated testing of setitem_from_int().

        TESTS::

            sage: c = ntl.ZZ_pContext(20)
            sage: x = ntl.ZZ_pX([2, 3, 4], c)
            sage: x._setitem_from_int_doctest(5, 42)
            sage: x
            [2 3 4 0 0 2]
        """
        self.c.restore_c()
        self.setitem_from_int(i, value)

    def __getitem__(self, long i):
        r"""
        Return the `i`-th entry of ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: x = ntl.ZZ_pX([2, 3, 4], c)
            sage: x[1]
            3
        """
        cdef ntl_ZZ_p r
        r = ntl_ZZ_p.__new__(ntl_ZZ_p)
        r.c = self.c
        self.c.restore_c()
        if i < 0:
            r.set_from_int(0)
        else:
            r.x = ZZ_pX_coeff( self.x, i)
        return r

    cdef int getitem_as_int(ntl_ZZ_pX self, long i) noexcept:
        r"""
        Return `i`-th coefficient as C int.

        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        self.c.restore_c()
        cdef ZZ_p_c r
        cdef long l = 0
        sig_on()
        r = ZZ_pX_coeff( self.x, i)
        ZZ_conv_to_long(l, ZZ_p_rep(r))
        sig_off()
        return l

    def _getitem_as_int_doctest(self, i):
        r"""
        This method exists solely for automated testing of getitem_as_int().

        TESTS::

            sage: c = ntl.ZZ_pContext(20)
            sage: x = ntl.ZZ_pX([2, 3, 5, -7, 11], c)
            sage: i = x._getitem_as_int_doctest(3)
            sage: i
            13
            sage: type(i)
            <... 'int'>
            sage: x._getitem_as_int_doctest(15)
            0
        """
        # self.c.restore_c() # restored in getitem_as_int()
        return self.getitem_as_int(i)

    def list(self):
        """
        Return list of entries as a list of ntl_ZZ_p.

        EXAMPLES::

            sage: x = ntl.ZZ_pX([1,3,5],11)
            sage: x.list()
            [1, 3, 5]
            sage: type(x.list()[0])
            <class 'sage.libs.ntl.ntl_ZZ_p.ntl_ZZ_p'>
        """
        # could be sped up.
        self.c.restore_c()
        cdef Py_ssize_t i
        return [self[i] for i from 0 <= i <= self.degree()]

    def __add__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: ntl.ZZ_pX(range(5),c) + ntl.ZZ_pX(range(6),c)
            [0 2 4 6 8 5]
        """
        if self.c is not other.c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        #self.c.restore_c() # restored in _new()
        ZZ_pX_add(r.x, self.x, other.x)
        sig_off()
        return r

    def __sub__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: ntl.ZZ_pX(range(5),c) - ntl.ZZ_pX(range(6),c)
            [0 0 0 0 0 15]
        """
        if self.c is not other.c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        #self.c.restore_c() # restored in _new()
        ZZ_pX_sub(r.x, self.x, other.x)
        sig_off()
        return r

    def __mul__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        EXAMPLES::

            sage: c1 = ntl.ZZ_pContext(20)
            sage: alpha = ntl.ZZ_pX(range(5), c1)
            sage: beta = ntl.ZZ_pX(range(6), c1)
            sage: alpha * beta
            [0 0 1 4 10 0 10 14 11]
            sage: c2 = ntl.ZZ_pContext(ntl.ZZ(5))  # we can mix up the moduli
            sage: x = ntl.ZZ_pX([2, 3, 4], c2)
            sage: x^2
            [4 2 0 4 1]
            sage: alpha * beta  # back to the first one and the ntl modulus gets reset correctly
            [0 0 1 4 10 0 10 14 11]
        """
        if self.c is not other.c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        # self.c.restore_c() # restored in _new()
        ZZ_pX_mul(r.x, self.x, other.x)
        sig_off()
        return r

    def __truediv__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        Compute quotient ``self / other``, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,3], c) * ntl.ZZ_pX([4,5], c)**2
            sage: g = ntl.ZZ_pX([4,5], c)
            sage: f/g
            [4 13 5 15]
            sage: ntl.ZZ_pX([1,2,3],c) * ntl.ZZ_pX([4,5],c)
            [4 13 5 15]

            sage: f = ntl.ZZ_pX(range(10), c); g = ntl.ZZ_pX([-1,0,1], c)
            sage: f/g
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0 1 2 3 4 5 6 7 8 9]) is not divisible by other (=[16 0 1])
        """
        if self.c is not other.c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef int divisible
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        #self.c.restore_c() # restored in _new()
        divisible = ZZ_pX_divide(r.x, self.x, other.x)
        sig_off()
        if not divisible:
            raise ArithmeticError("self (=%s) is not divisible by other (=%s)" % (self, other))
        return r

    def __mod__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        Given polynomials a, b in ZZ_p[X], if p is prime, then there exist polynomials q, r
        in ZZ_p[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns r.

        If p is not prime this function may raise a :exc:`RuntimeError`
        due to division by a noninvertible element of ZZ_p.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([2,4,6], c); g = ntl.ZZ_pX([2], c)
            sage: f % g   # 0
            []
            sage: f = ntl.ZZ_pX(range(10), c); g = ntl.ZZ_pX([-1,0,1], c)
            sage: f % g
            [3 8]

            sage: c = ntl.ZZ_pContext(ntl.ZZ(16))
            sage: f = ntl.ZZ_pX([2,4,6], c); g = ntl.ZZ_pX([2], c)
            sage: f % g
            Traceback (most recent call last):
            ...
            NTLError: ZZ_p: division by non-invertible element
        """
        if self.c is not other.c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        #self.c.restore_c() # restored in _new()
        ZZ_pX_rem(r.x, self.x, other.x)
        sig_off()
        return r

    def quo_rem(self, ntl_ZZ_pX other):
        """
        Return the unique integral `q` and `r` such that ``self = q*other +
        r``, if they exist.  Otherwise raises an Exception.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX(range(10), c); g = ntl.ZZ_pX([-1,0,1], c)
            sage: q, r = f.quo_rem(g)
            sage: q, r
            ([3 7 1 4 14 16 8 9], [3 8])
            sage: q*g + r == f
            True
        """
        if self.c is not other.c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef ntl_ZZ_pX r = self._new()
        cdef ntl_ZZ_pX q = self._new()
        sig_on()
        #self.c.restore_c() # restored in _new()
        ZZ_pX_DivRem(q.x, r.x, self.x, other.x)
        sig_off()
        return q,r

    def square(self):
        """
        Return f*f.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([-1,0,1], c)
            sage: f*f
            [1 0 15 0 1]
        """
        #self.c.restore_c() # restored in _new()
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_sqr(r.x, self.x)
        sig_off()
        return r

    def __pow__(self, n, modulus):
        """
        Return the ``n``-th nonnegative power of ``self``.

        If ``modulus`` is not ``None``, the exponentiation is performed
        modulo the polynomial ``modulus``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: g = ntl.ZZ_pX([-1,0,1],c)
            sage: g**10
            [1 0 10 0 5 0 0 0 10 0 8 0 10 0 0 0 5 0 10 0 1]

            sage: x = ntl.ZZ_pX([0,1],c)
            sage: x**10
            [0 0 0 0 0 0 0 0 0 0 1]

        Modular exponentiation::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,0,1],c)
            sage: m = ntl.ZZ_pX([1,0,0,0,0,1],c)
            sage: pow(f, 123**45, m)
            [1 19 3 0 3]

        Modular exponentiation of ``x``::

            sage: f = ntl.ZZ_pX([0,1],c)
            sage: f.is_x()
            True
            sage: m = ntl.ZZ_pX([1,1,0,0,0,1],c)
            sage: pow(f, 123**45, m)
            [15 5 5 11]
        """
        if modulus is None:
            return (<ntl_ZZ_pX>self)._pow(n)
        else:
            return (<ntl_ZZ_pX>self)._powmod(Integer(n), modulus)

    cdef ntl_ZZ_pX _pow(ntl_ZZ_pX self, long n):
        """
        Compute the ``n``-th power of ``self``.
        """
        if n < 0:
            raise NotImplementedError
        cdef long ln = n
        #self.c.restore_c() # restored in _new()
        cdef ntl_ZZ_pX r = self._new()
        if self.is_x() and n >= 1:
            ZZ_pX_LeftShift(r.x, self.x, n-1)
        else:
            sig_on()
            ZZ_pX_power(r.x, self.x, ln)
            sig_off()
        return r

    cdef ntl_ZZ_pX _powmod(ntl_ZZ_pX self, Integer n, ntl_ZZ_pX modulus):
        r"""
        Compute the ``n``-th power of ``self`` modulo a polynomial.
        """
        cdef ntl_ZZ_pX r = self._new()
        cdef ZZ_c n_ZZ
        cdef ZZ_pX_Modulus_c mod
        is_x = self.is_x()
        sig_on()
        mpz_to_ZZ(&n_ZZ, (<Integer>n).value)
        ZZ_pX_Modulus_build(mod, modulus.x)
        if is_x:
            ZZ_pX_PowerXMod_pre(r.x, n_ZZ, mod)
        else:
            ZZ_pX_PowerMod_pre(r.x, self.x, n_ZZ, mod)
        sig_off()
        return r

    def __richcmp__(ntl_ZZ_pX self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,3],c)
            sage: g = ntl.ZZ_pX([1,2,3,0],c)
            sage: f == g
            True
            sage: g = ntl.ZZ_pX([0,1,2,3],c)
            sage: f != g
            True
            sage: f != 0
            True
        """
        self.c.restore_c()

        if op != Py_EQ and op != Py_NE:
            raise TypeError("polynomials are not ordered")

        cdef ntl_ZZ_pX b
        try:
            b = <ntl_ZZ_pX?>other
        except TypeError:
            b = ntl_ZZ_pX(other, self.c)

        return (op == Py_EQ) == (self.x == b.x)

    def is_zero(self):
        """
        Return ``True`` exactly if this polynomial is 0.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([0,0,0,20],c)
            sage: f.is_zero()
            True
            sage: f = ntl.ZZ_pX([0,0,1],c)
            sage: f
            [0 0 1]
            sage: f.is_zero()
            False
        """
        self.c.restore_c()
        return bool(ZZ_pX_IsZero(self.x))

    def is_one(self):
        """
        Return ``True`` exactly if this polynomial is 1.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,1],c)
            sage: f.is_one()
            False
            sage: f = ntl.ZZ_pX([1],c)
            sage: f.is_one()
            True
        """
        self.c.restore_c()
        return bool(ZZ_pX_IsOne(self.x))

    def is_monic(self):
        """
        Return ``True`` exactly if this polynomial is monic.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([2,0,0,1],c)
            sage: f.is_monic()
            True
            sage: g = f.reverse()
            sage: g.is_monic()
            False
            sage: g
            [1 0 0 2]
            sage: f = ntl.ZZ_pX([1,2,0,3,0,2],c)
            sage: f.is_monic()
            False
        """
        self.c.restore_c()
        if ZZ_pX_IsZero(self.x):
            return False
        return bool(ZZ_p_IsOne(ZZ_pX_LeadCoeff(self.x)))

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([2,0,0,1],c)
            sage: -f
            [18 0 0 19]
        """
        cdef ntl_ZZ_pX r = self._new()
        self.c.restore_c()
        ZZ_pX_negate(r.x, self.x)
        return r

    def left_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the left n positions.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([2,0,0,1],c)
            sage: f
            [2 0 0 1]
            sage: f.left_shift(2)
            [0 0 2 0 0 1]
            sage: f.left_shift(5)
            [0 0 0 0 0 2 0 0 1]

        A negative left shift is a right shift.
            sage: f.left_shift(-2)
            [0 1]
        """
        #self.c.restore_c() # restored in _new()
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_LeftShift(r.x, self.x, n)
        sig_off()
        return r

    def right_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the right n positions.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([2,0,0,1],c)
            sage: f
            [2 0 0 1]
            sage: f.right_shift(2)
            [0 1]
            sage: f.right_shift(5)
            []
            sage: f.right_shift(-2)
            [0 0 2 0 0 1]
        """
        #self.c.restore_c() # restored in _new()
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_RightShift(r.x, self.x, n)
        sig_off()
        return r

    def _left_pshift(self, ntl_ZZ n):
        """
        Multiply all coefficients by `n` and the context by `n`.
        """
        cdef ntl_ZZ new_c_p = ntl_ZZ.__new__(ntl_ZZ)
        ZZ_mul(new_c_p.x, (<ntl_ZZ>self.c.p).x, n.x)
        cdef ntl_ZZ_pContext_class new_c = <ntl_ZZ_pContext_class>ntl_ZZ_pContext(new_c_p)
        new_c.restore_c()
        cdef ntl_ZZ_pX ans = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        ans.c = new_c
        ZZ_pX_left_pshift(ans.x, self.x, n.x, new_c.x)
        return ans

    def _right_pshift(self, ntl_ZZ n):
        """
        Divide all coefficients by `n` and the context by `n`.  Only really
        makes sense when `n` divides ``self.c.p``.
        """
        cdef ntl_ZZ new_c_p = ntl_ZZ.__new__(ntl_ZZ)
        ZZ_div(new_c_p.x, (<ntl_ZZ>self.c.p).x, n.x)
        cdef ntl_ZZ_pContext_class new_c = <ntl_ZZ_pContext_class>ntl_ZZ_pContext(new_c_p)
        new_c.restore_c()
        cdef ntl_ZZ_pX ans = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        ans.c = new_c
        ZZ_pX_right_pshift(ans.x, self.x, n.x, new_c.x)
        return ans

    def gcd(self, ntl_ZZ_pX other):
        """
        Return the gcd d = gcd(a, b), where by convention the leading coefficient
        of d is >= 0.  We uses a multi-modular algorithm.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,3],c) * ntl.ZZ_pX([4,5],c)**2
            sage: g = ntl.ZZ_pX([1,1,1],c)**3 * ntl.ZZ_pX([1,2,3],c)
            sage: f.gcd(g)
            [6 12 1]
            sage: g.gcd(f)
            [6 12 1]
        """
        #self.c.restore_c() # restored in _new()
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_GCD(r.x, self.x, other.x)
        sig_off()
        return r

    def xgcd(self, ntl_ZZ_pX other, plain=True):
        """
        Return `r`, `s`, `t` such that ``r = s*self + t*other``.

        Here r is the resultant of a and b; if r != 0, then this
        function computes s and t such that: a*s + b*t = r; otherwise
        s and t are both 0.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,3],c) * ntl.ZZ_pX([4,5],c)**2
            sage: g = ntl.ZZ_pX([1,1,1],c)**3 * ntl.ZZ_pX([1,2,3],c)
            sage: f.xgcd(g)   # nothing since they are not coprime
            ([6 12 1], [15 13 6 8 7 9], [4 13])

        In this example the input quadratic polynomials have a common
        root modulo 13::

            sage: f = ntl.ZZ_pX([5,0,1],c)
            sage: g = ntl.ZZ_pX([18,0,1],c)
            sage: f.xgcd(g)
            ([1], [13], [4])
            """
        self.c.restore_c()
        cdef ntl_ZZ_pX s = self._new()
        cdef ntl_ZZ_pX t = self._new()
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        if plain:
            ZZ_pX_PlainXGCD(r.x, s.x, t.x, self.x, other.x)
        else:
            ZZ_pX_XGCD(r.x, s.x, t.x, self.x, other.x)
        sig_off()
        return (r,s,t)

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([5,0,1],c)
            sage: f.degree()
            2
            sage: f = ntl.ZZ_pX(range(100),c)
            sage: f.degree()
            99
            sage: f = ntl.ZZ_pX([], c)
            sage: f.degree()
            -1
            sage: f = ntl.ZZ_pX([1],c)
            sage: f.degree()
            0
        """
        self.c.restore_c()
        return ZZ_pX_deg(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([3,6,9],c)
            sage: f.leading_coefficient()
            9
            sage: f = ntl.ZZ_pX([],c)
            sage: f.leading_coefficient()
            0
        """
        self.c.restore_c()
        #return ZZ_pX_LeadCoeff(self.x)
        cdef long i
        i = ZZ_pX_deg(self.x)
        return self[i]

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([3,6,9],c)
            sage: f.constant_term()
            3
            sage: f = ntl.ZZ_pX([],c)
            sage: f.constant_term()
            0
        """
        self.c.restore_c()
        #return ZZ_pX_ConstTerm(self.x)
        return self[0]

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([],c)
            sage: f.set_x()
            sage: f
            [0 1]
            sage: g = ntl.ZZ_pX([0,1],c)
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory::

            sage: f is g
            False
        """
        self.c.restore_c()
        ZZ_pX_SetX(self.x)
        #ZZ_pX_clear(self.x)
        #ZZ_pX_SetCoeff_long(self.x, 1, 1)

    def is_x(self):
        """
        ``True`` if this is the polynomial "x".

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([],c)
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = ntl.ZZ_pX([0,1],c)
            sage: f.is_x()
            True
            sage: f = ntl.ZZ_pX([1],c)
            sage: f.is_x()
            False
        """
        return bool(ZZ_pX_IsX(self.x))

    def convert_to_modulus(self, ntl_ZZ_pContext_class c):
        """
        Return a new ntl_ZZ_pX which is the same as self, but considered modulo a different p.

        In order for this to make mathematical sense, c.p should divide self.c.p
        (in which case ``self`` is reduced modulo c.p) or self.c.p should divide c.p
        (in which case ``self`` is lifted to something modulo c.p congruent to ``self`` modulo self.c.p)

        EXAMPLES::

            sage: a = ntl.ZZ_pX([412,181,991],5^4)
            sage: a
            [412 181 366]
            sage: b = ntl.ZZ_pX([198,333,91],5^4)
            sage: ap = a.convert_to_modulus(ntl.ZZ_pContext(5^2))
            sage: bp = b.convert_to_modulus(ntl.ZZ_pContext(5^2))
            sage: ap
            [12 6 16]
            sage: bp
            [23 8 16]
            sage: ap*bp
            [1 9 8 24 6]
            sage: (a*b).convert_to_modulus(ntl.ZZ_pContext(5^2))
            [1 9 8 24 6]
        """
        c.restore_c()
        cdef ntl_ZZ_pX ans = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
        ZZ_pX_conv_modulus(ans.x, self.x, c.x)
        ans.c = c
        return ans

    def derivative(self):
        """
        Return the derivative of this polynomial.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,7,0,13],c)
            sage: f.derivative()
            [7 0 19]
        """
        cdef ntl_ZZ_pX r = self._new()  # restores context
        sig_on()
        ZZ_pX_diff(r.x, self.x)
        sig_off()
        return r

    def factor(self, verbose=False):
        """
        Return the factorization of ``self``. Assumes ``self`` is
        monic.

        NOTE: The roots are returned in a random order.

        EXAMPLES:

        These computations use pseudo-random numbers, so we set the
        seed for reproducible testing. ::

            sage: set_random_seed(12)
            sage: ntl.ZZ_pX([-1,0,0,0,0,1],5).factor()
            [([4 1], 5)]
            sage: ls = ntl.ZZ_pX([-1,0,0,0,1],5).factor()
            sage: ls
            [([4 1], 1), ([2 1], 1), ([1 1], 1), ([3 1], 1)]
            sage: prod( [ x[0] for x in ls ] )
            [4 0 0 0 1]
            sage: ntl.ZZ_pX([3,7,0,1], 31).factor()
            [([3 7 0 1], 1)]
            sage: ntl.ZZ_pX([3,7,1,8], 28).factor()
            Traceback (most recent call last):
            ...
            ValueError: self must be monic.
        """
        current_randstate().set_seed_ntl(False)

        self.c.restore_c()
        cdef ZZ_pX_c** v
        cdef ntl_ZZ_pX r
        cdef long* e
        cdef long i, n
        if not self.is_monic():
            raise ValueError("self must be monic.")
        sig_on()
        ZZ_pX_factor(&v, &e, &n, &self.x, verbose)
        sig_off()
        F = []
        for i from 0 <= i < n:
            r = self._new()
            r.x = v[i][0]
            F.append((r, e[i]))
            #F.append((make_ZZ_pX(v[i], self.c), e[i]))
        for i from 0 <= i < n:
            del v[i]
        sig_free(v)
        sig_free(e)
        return F

    def linear_roots(self):
        """
        Assumes that input is monic, and has deg(f) roots.
        Returns the list of roots.

        NOTE: This function will go into an infinite loop if you
        give it a polynomial without deg(f) linear factors. Note
        also that the result is not deterministic, i.e. the
        order of the roots returned is random.

        EXAMPLES::

            sage: ntl.ZZ_pX([-1,0,0,0,0,1],5).linear_roots()
            [1, 1, 1, 1, 1]
            sage: roots = ntl.ZZ_pX([-1,0,0,0,1],5).linear_roots()
            sage: [ ntl.ZZ_p(i,5) in roots for i in [1..4] ]
            [True, True, True, True]
            sage: ntl.ZZ_pX([3,7,1,8], 28).linear_roots()
            Traceback (most recent call last):
            ...
            ValueError: self must be monic.
        """
        current_randstate().set_seed_ntl(False)

        self.c.restore_c()
        cdef ZZ_p_c** v
        cdef ntl_ZZ_p r
        cdef long i, n
        if not self.is_monic():
            raise ValueError("self must be monic.")
        sig_on()
        ZZ_pX_linear_roots(&v, &n, &self.x)
        sig_off()
        F = []
        for i from 0 <= i < n:
            r = ntl_ZZ_p.__new__(ntl_ZZ_p)
            r.c = self.c
            r.x = v[i][0]
            F.append(r)
            #F.append(make_ZZ_p(v[i], self.c))
        for i from 0 <= i < n:
            del v[i]
        sig_free(v)
        return F

    def reverse(self, hi=None):
        """
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.  If hi is set then this function behaves
        as if this polynomial has degree hi.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,3,4,5],c)
            sage: f.reverse()
            [5 4 3 2 1]
            sage: f.reverse(hi=10)
            [0 0 0 0 0 0 5 4 3 2 1]
            sage: f.reverse(hi=2)
            [3 2 1]
            sage: f.reverse(hi=-2)
            []
        """
        cdef ntl_ZZ_pX r = self._new()
        if hi is None:
            ZZ_pX_reverse(r.x, self.x)
        else:
            ZZ_pX_reverse_hi(r.x, self.x, int(hi))
        return r

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,3,4,5],c)
            sage: f.truncate(3)
            [1 2 3]
            sage: f.truncate(8)
            [1 2 3 4 5]
            sage: f.truncate(1)
            [1]
            sage: f.truncate(0)
            []
            sage: f.truncate(-1)
            []
            sage: f.truncate(-5)
            []
        """
        cdef ntl_ZZ_pX r = self._new()
        if m <= 0:
            return r
        sig_on()
        ZZ_pX_trunc(r.x, self.x, m)
        sig_off()
        return r

    def multiply_and_truncate(self, ntl_ZZ_pX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,3,4,5],c)
            sage: g = ntl.ZZ_pX([10],c)
            sage: f.multiply_and_truncate(g, 2)
            [10]
            sage: g.multiply_and_truncate(f, 2)
            [10]
        """
        cdef ntl_ZZ_pX r = self._new()
        if m <= 0:
            return r
        sig_on()
        ZZ_pX_MulTrunc(r.x, self.x, other.x, m)
        sig_off()
        return r

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,3,4,5],c)
            sage: f.square_and_truncate(4)
            [1 4 10]
            sage: (f*f).truncate(4)
            [1 4 10]
        """
        cdef ntl_ZZ_pX r = self._new()
        if m <= 0:
            return r
        sig_on()
        ZZ_pX_SqrTrunc(r.x, self.x, m)
        sig_off()
        return r

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of ``self`` modulo `x^m`.

        The constant term of ``self`` must be a unit.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,3,4,5,6,7],c)
            sage: f.invert_and_truncate(20)
            [1 18 1 0 0 0 0 8 17 2 13 0 0 0 4 0 17 10 9]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 4 4 3]
        """
        if m < 0:
            raise ArithmeticError("m (=%s) must be positive" % m)
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_InvTrunc(r.x, self.x, m)
        sig_off()
        return r

    def invmod(self, ntl_ZZ_pX modulus):
        """
        Return the inverse of ``self`` modulo the modulus using NTL's InvMod.
        """
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_InvMod(r.x, self.x, modulus.x)
        sig_off()
        return r

    def invmod_newton(self, ntl_ZZ_pX modulus):
        """
        Return the inverse of ``self`` modulo the modulus using Newton lifting.
        Only works if modulo a power of a prime, and if modulus is either
        unramified or Eisenstein.
        """
        cdef Integer pn = Integer(self.c.p)
        cdef ntl_ZZ_pX ans = self._new()
        F = pn.factor()
        if len(F) > 1:
            raise ValueError("must be modulo a prime power")
        p = F[0][0]
        cdef ntl_ZZ pZZ = <ntl_ZZ>ntl_ZZ(p)
        cdef ZZ_pX_Modulus_c mod
        ZZ_pX_Modulus_build(mod, modulus.x)
        cdef ntl_ZZ_pX mod_prime
        cdef ntl_ZZ_pContext_class ctx
        cdef long mini = 0, minval = 0
        if Integer(modulus[0].lift()).valuation(p) == 1:
            eisenstein = True
            for c in modulus.list()[1:-1]:
                if not p.divides(Integer(c.lift())):
                    eisenstein = False
                    break
            if eisenstein:
                if p.divides(Integer(self[0])):
                    raise ZeroDivisionError("cannot invert element")
                ZZ_pX_InvMod_newton_ram(ans.x, self.x, mod, self.c.x)
            else:
                raise ValueError("not eisenstein or unramified")
        else:
            ctx = <ntl_ZZ_pContext_class>ntl_ZZ_pContext(p)
            mod_prime = ntl_ZZ_pX.__new__(ntl_ZZ_pX)
            ZZ_pX_conv_modulus(mod_prime.x, modulus.x, ctx.x)
            mod_prime.c = ctx
            F = mod_prime.factor()
            if len(F) == 1 and F[0][1] == 1:
                ZZ_pX_min_val_coeff(minval, mini, self.x, pZZ.x)
                if minval > 0:
                    raise ZeroDivisionError("cannot invert element")
                ZZ_pX_InvMod_newton_unram(ans.x, self.x, mod, self.c.x, ctx.x)
            else:
                raise ValueError("not eisenstein or unramified")
        return ans

    def multiply_mod(self, ntl_ZZ_pX other, ntl_ZZ_pX modulus):
        """
        Return ``self*other % modulus``.  The modulus must be monic with
        ``deg(modulus) > 0``, and ``self`` and ``other`` must have smaller degree.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(ntl.ZZ(20))
            sage: modulus = ntl.ZZ_pX([1,2,0,1],c)    # must be monic
            sage: g = ntl.ZZ_pX([-1,0,1],c)
            sage: h = ntl.ZZ_pX([3,7,13],c)
            sage: h.multiply_mod(g, modulus)
            [10 6 4]
        """
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_MulMod(r.x, self.x, other.x, modulus.x)
        sig_off()
        return r

    def trace_mod(self, ntl_ZZ_pX modulus):
        """
        Return the trace of this polynomial modulus the modulus.
        The modulus must be monic, and of positive degree bigger
        than the degree of ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,0,3],c)
            sage: mod = ntl.ZZ_pX([5,3,-1,1,1],c)
            sage: f.trace_mod(mod)
            3
        """
        self.c.restore_c()
        cdef ntl_ZZ_p r = ntl_ZZ_p.__new__(ntl_ZZ_p)
        r.c = self.c
        sig_on()
        ZZ_pX_TraceMod(r.x, self.x, modulus.x)
        sig_off()
        return r

    def trace_list(self):
        """
        Return the list of traces of the powers `x^i` of the
        monomial x modulo this polynomial for i = 0, ..., deg(f)-1.

        This polynomial must be monic.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,0,3,0,1],c)
            sage: f.trace_list()
            [5, 0, 14, 0, 10]

        The input polynomial must be monic or a :exc:`ValueError` is raised::

            sage: c = ntl.ZZ_pContext(20)
            sage: f = ntl.ZZ_pX([1,2,0,3,0,2],c)
            sage: f.trace_list()
            Traceback (most recent call last):
            ...
            ValueError: polynomial must be monic.
        """
        # This function should be redone to use TraceVec, which requires improving the wrapper for vec_ZZ_p
        self.c.restore_c()
        if not self.is_monic():
            raise ValueError("polynomial must be monic.")
        sig_on()
        cdef char* t
        t = ZZ_pX_trace_list(&self.x)
        r = eval(char_to_str(t).replace(' ', ','))
        string_delete(t)
        return r

    def resultant(self, ntl_ZZ_pX other):
        """
        Return the resultant of ``self`` and ``other``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([17,0,1,1],c)
            sage: g = ntl.ZZ_pX([34,-17,18,2],c)
            sage: f.resultant(g)
            0
        """
        self.c.restore_c()
        cdef ntl_ZZ_p r = ntl_ZZ_p.__new__(ntl_ZZ_p)
        r.c = self.c
        sig_on()
        ZZ_pX_resultant(r.x, self.x, other.x)
        sig_off()
        return r

    def norm_mod(self, ntl_ZZ_pX modulus):
        """
        Return the norm of this polynomial modulo the modulus.  The
        modulus must be monic, and of positive degree strictly greater
        than the degree of ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,0,3],c)
            sage: mod = ntl.ZZ_pX([-5,2,0,0,1],c)
            sage: f.norm_mod(mod)
            11

        The norm is the constant term of the characteristic polynomial::

            sage: f.charpoly_mod(mod)
            [11 1 8 14 1]
        """
        self.c.restore_c()
        cdef ntl_ZZ_p r = ntl_ZZ_p.__new__(ntl_ZZ_p)
        r.c = self.c
        sig_on()
        ZZ_pX_NormMod(r.x, self.x, modulus.x)
        sig_off()
        return r

    def discriminant(self):
        r"""
        Return the discriminant of a=self, which is by definition
        $$
                (-1)^{m(m-1)/2} {\mbox{\tt resultant}}(a, a')/lc(a),
        $$
        where m = deg(a), and lc(a) is the leading coefficient of a.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,0,3],c)
            sage: f.discriminant()
            1
        """
        self.c.restore_c()
        cdef long m

        c = ~self.leading_coefficient()
        m = self.degree()
        if (m*(m-1) // 2) % 2:
            c = -c
        return c*self.resultant(self.derivative())

    def charpoly_mod(self, ntl_ZZ_pX modulus):
        """
        Return the characteristic polynomial of this polynomial modulo
        the modulus.  The modulus must be monic of degree bigger than
        ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,0,3],c)
            sage: mod = ntl.ZZ_pX([-5,2,0,0,1],c)
            sage: f.charpoly_mod(mod)
            [11 1 8 14 1]
        """
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_CharPolyMod(r.x, self.x, modulus.x)
        sig_off()
        return r

    def minpoly_mod(self, ntl_ZZ_pX modulus):
        """
        Return the minimal polynomial of this polynomial modulo the
        modulus.  The modulus must be monic of degree bigger than
        ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([0,0,1],c)
            sage: g = f*f
            sage: f.charpoly_mod(g)
            [0 0 0 0 1]

        However, since `f^2 = 0` modulo `g`, its minimal polynomial
        is of degree `2`::

            sage: f.minpoly_mod(g)
            [0 0 1]
        """
        cdef ntl_ZZ_pX r = self._new()
        sig_on()
        ZZ_pX_MinPolyMod(r.x, self.x, modulus.x)
        sig_off()
        return r

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,3],c)
            sage: f
            [1 2 3]
            sage: f.clear()
            sage: f
            []
        """
        self.c.restore_c()
        ZZ_pX_clear(self.x)

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,3],c)
            sage: f.preallocate_space(20)
            sage: f
            [1 2 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1 2 3 0 0 0 0 0 0 0 5]
        """
        self.c.restore_c()
        sig_on()
        self.x.SetMaxLength(n)
        #ZZ_pX_preallocate_space(&self.x, n)
        sig_off()

    def compose_mod(self, ntl_ZZ_pX other, ntl_ZZ_pX modulus):
        r"""
        Compute `f(g) \bmod h`.

        To be precise about the order fo compostion, given ``self``, ``other``
        and ``modulus`` as `f(x)`, `g(x)` and `h(x)` compute `f(g(x)) \bmod h(x)`.

        INPUT:

        - ``other`` -- a polynomial `g(x)`
        - ``modulus`` -- a polynomial `h(x)`

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(17)
            sage: f = ntl.ZZ_pX([1,2,3],c)
            sage: g = ntl.ZZ_pX([3,2,1],c)
            sage: h = ntl.ZZ_pX([1,0,1],c)
            sage: f.compose_mod(g, h)
            [5 11]

        AUTHORS:

        - Giacomo Pope (2024-08) initial implementation
        """
        cdef ntl_ZZ_pX r = self._new()
        cdef ZZ_pX_Modulus_c mod

        sig_on()
        # NTL requires that g is reduced
        other = other % modulus

        # Build the modulus type
        ZZ_pX_Modulus_build(mod, modulus.x)

        # Compute f(g) % h
        ZZ_pX_CompMod(r.x, self.x, other.x, mod)
        sig_off()
        return r

cdef class ntl_ZZ_pX_Modulus():
    """
    Thin holder for ZZ_pX_Moduli.
    """
    def __cinit__(self, ntl_ZZ_pX poly):
        ZZ_pX_Modulus_build(self.x, poly.x)
        self.poly = poly

    def __repr__(self):
        return "NTL ZZ_pXModulus %s (mod %s)"%(self.poly, self.poly.c.p)

    def degree(self):
        cdef Integer ans = Integer.__new__(Integer)
        mpz_set_ui(ans.value, ZZ_pX_Modulus_deg(self.x))
        return ans

    ## TODO: NTL's ZZ_pX has minpolys of linear recurrence sequences!!!
