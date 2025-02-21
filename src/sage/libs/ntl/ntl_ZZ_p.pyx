# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

#*****************************************************************************
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
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from sage.ext.cplusplus cimport ccrepr, ccreadstr

include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer cimport Integer
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.rings.rational cimport Rational

from sage.libs.ntl.ntl_ZZ import unpickle_class_args
from sage.libs.ntl.convert cimport PyLong_to_ZZ, mpz_to_ZZ

from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pContext import ntl_ZZ_pContext

from sage.misc.randstate cimport current_randstate


ZZ_sage = IntegerRing()


def ntl_ZZ_p_random_element(v):
    """
    Return a random number modulo p.

    EXAMPLES::

        sage: a = sage.libs.ntl.ntl_ZZ_p.ntl_ZZ_p_random_element(17)
        sage: type(a)
        <class 'sage.libs.ntl.ntl_ZZ_p.ntl_ZZ_p'>
        sage: a.modulus()
        17
    """
    current_randstate().set_seed_ntl(False)

    cdef ntl_ZZ_p y
    v = ntl_ZZ_pContext(v)
    y = ntl_ZZ_p(0,v)
    sig_on()
    ZZ_p_random(y.x)
    sig_off()
    return y


##############################################################################
#
# ZZ_p_c: integers modulo p
#
##############################################################################
cdef class ntl_ZZ_p():
    r"""
    The \class{ZZ_p} class is used to represent integers modulo `p`.
    The modulus `p` may be any positive integer, not necessarily prime.

    Objects of the class \class{ZZ_p} are represented as a \code{ZZ} in the
    range `0, \ldots, p-1`.

    Each \class{ZZ_p} contains a pointer of a \class{ZZ_pContext} which
    contains pre-computed data for NTL.  These can be explicitly constructed
    and passed to the constructor of a \class{ZZ_p} or the \class{ZZ_pContext}
    method \code{ZZ_p} can be used to construct a \class{ZZ_p} element.

    This class takes care of making sure that the C++ library NTL global
    variable is set correctly before performing any arithmetic.
    """

    def __init__(self, v=None, modulus=None):
        r"""
        Initialize an NTL integer mod p.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(11)
            sage: ntl.ZZ_p(12r, c)
            1
            sage: ntl.ZZ_p(Integer(95413094), c)
            7
            sage: ntl.ZZ_p('-1', c)
            10

        AUTHOR: Joel B. Mohler (2007-06-14)
        """
        if modulus is None:
            raise ValueError("You must specify a modulus when creating a ZZ_p.")

        # self.c._assert_is_current_modulus()  # The context was restored in __new__

        cdef ZZ_c temp, num, den
        cdef long failed
        if v is not None:
            if isinstance(v, ntl_ZZ_p):
                self.x = (<ntl_ZZ_p>v).x
            elif isinstance(v, int):
                PyLong_to_ZZ(&temp, v)
                self.x = ZZ_to_ZZ_p(temp)
            elif isinstance(v, Integer):
                mpz_to_ZZ(&temp, (<Integer>v).value)
                self.x = ZZ_to_ZZ_p(temp)
            elif isinstance(v, Rational):
                mpz_to_ZZ(&num, (<Integer>v.numerator()).value)
                mpz_to_ZZ(&den, (<Integer>v.denominator()).value)
                ZZ_p_div(self.x, ZZ_to_ZZ_p(num), ZZ_to_ZZ_p(den))
            else:
                str_v = str(v)  # can cause modulus to change; Issue #25790
                self.c.restore_c()
                ccreadstr(self.x, str_v)

    def __cinit__(self, v=None, modulus=None):
        #################### WARNING ###################
        ## Before creating a ZZ_p, you must create a  ##
        ## ZZ_pContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing an ntl_ZZ_p          ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = ntl_ZZ_p.__new__(ntl_ZZ_p) without
        ## first restoring a ZZ_pContext, which could ##
        ## have unfortunate consequences.  See _new  ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None:
            return
        if isinstance(modulus, ntl_ZZ_pContext_class):
            self.c = <ntl_ZZ_pContext_class>modulus
            self.c.restore_c()
        else:
            self.c = <ntl_ZZ_pContext_class>ntl_ZZ_pContext(modulus)
            self.c.restore_c()

    cdef ntl_ZZ_p _new(self):
        cdef ntl_ZZ_p r
        self.c.restore_c()
        r = ntl_ZZ_p.__new__(ntl_ZZ_p)
        r.c = self.c
        return r

    def __reduce__(self):
        """
        EXAMPLES::

            sage: a = ntl.ZZ_p(4,7)
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_class_args, (ntl_ZZ_p, (self.lift(), self.modulus_context()))

    def modulus_context(self):
        """
        Return the modulus for ``self``.

        EXAMPLES::

            sage: x = ntl.ZZ_p(5,17)
            sage: c = x.modulus_context()
            sage: y = ntl.ZZ_p(3,c)
            sage: x+y
            8
            sage: c == y.modulus_context()
            True
            sage: c == ntl.ZZ_p(7,17).modulus_context()
            True
        """
        return self.c

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ntl.ZZ_p(7,192).__repr__()
            '7'
        """
        self.c.restore_c()
        return ccrepr(self.x)

    def __richcmp__(ntl_ZZ_p self, other, int op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(11)
            sage: ntl.ZZ_p(12r, c) == ntl.ZZ_p(1, c)
            True
            sage: ntl.ZZ_p(35r, c) == ntl.ZZ_p(1, c)
            False
            sage: "2" == ntl.ZZ_p(35r, c)
            True
            sage: ntl.ZZ_p(35r, c) == 2
            True
        """
        self.c.restore_c()

        if op != Py_EQ and op != Py_NE:
            raise TypeError("integers mod p are not ordered")

        cdef ntl_ZZ_p b
        try:
            b = <ntl_ZZ_p?>other
        except TypeError:
            b = ntl_ZZ_p(other, self.c)

        return (op == Py_EQ) == (self.x == b.x)

    def __invert__(ntl_ZZ_p self):
        r"""
        EXAMPLES::

            sage: c=ntl.ZZ_pContext(11)
            sage: ~ntl.ZZ_p(2r,modulus=c)
            6
        """
        cdef ntl_ZZ_p r = self._new()
        sig_on()
        self.c.restore_c()
        ZZ_p_inv(r.x, self.x)
        sig_off()
        return r

    def __mul__(ntl_ZZ_p self, other):
        """
        EXAMPLES::

            sage: x = ntl.ZZ_p(5,31) ; y = ntl.ZZ_p(8,31)
            sage: x*y ## indirect doctest
            9
        """
        cdef ntl_ZZ_p y
        cdef ntl_ZZ_p r = self._new()
        if not isinstance(other, ntl_ZZ_p):
            other = ntl_ZZ_p(other,self.c)
        elif self.c is not (<ntl_ZZ_p>other).c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        y = other
        self.c.restore_c()
        ZZ_p_mul(r.x, self.x, y.x)
        return r

    def __sub__(ntl_ZZ_p self, other):
        """
        EXAMPLES::

            sage: x = ntl.ZZ_p(5,31) ; y = ntl.ZZ_p(8,31)
            sage: x-y ## indirect doctest
            28
            sage: y-x
            3
        """
        if not isinstance(other, ntl_ZZ_p):
            other = ntl_ZZ_p(other,self.c)
        elif self.c is not (<ntl_ZZ_p>other).c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        cdef ntl_ZZ_p r = self._new()
        self.c.restore_c()
        ZZ_p_sub(r.x, self.x, (<ntl_ZZ_p>other).x)
        return r

    def __add__(ntl_ZZ_p self, other):
        """
        EXAMPLES::

            sage: x = ntl.ZZ_p(5,31) ; y = ntl.ZZ_p(8,31)
            sage: x+y ## indirect doctest
            13
        """
        cdef ntl_ZZ_p y
        cdef ntl_ZZ_p r = ntl_ZZ_p(modulus=self.c)
        if not isinstance(other, ntl_ZZ_p):
            other = ntl_ZZ_p(other,modulus=self.c)
        elif self.c is not (<ntl_ZZ_p>other).c:
            raise ValueError("You cannot perform arithmetic with elements of different moduli.")
        y = other
        sig_on()
        self.c.restore_c()
        ZZ_p_add(r.x, self.x, y.x)
        sig_off()
        return r

    def __neg__(ntl_ZZ_p self):
        """
        EXAMPLES::

            sage: x = ntl.ZZ_p(5,31)
            sage: -x ## indirect doctest
            26
        """
        cdef ntl_ZZ_p r = ntl_ZZ_p(modulus=self.c)
        sig_on()
        self.c.restore_c()
        ZZ_p_negate(r.x, self.x)
        sig_off()
        return r

    def __pow__(ntl_ZZ_p self, long e, ignored):
        """
        EXAMPLES::

            sage: x = ntl.ZZ_p(5,31)
            sage: x**3 ## indirect doctest
            1
        """
        cdef ntl_ZZ_p r = ntl_ZZ_p(modulus=self.c)
        sig_on()
        self.c.restore_c()
        ZZ_p_power(r.x, self.x, e)
        sig_off()
        return r

    def __int__(self):
        """
        Return ``self`` as an int.

        EXAMPLES::

            sage: x = ntl.ZZ_p(3,8)
            sage: x.__int__()
            3
            sage: type(x.__int__())
            <... 'int'>
        """
        return self.get_as_int()

    cdef int get_as_int(ntl_ZZ_p self) noexcept:
        r"""
        Return value as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        self.c.restore_c()
        return ZZ_p_to_int(self.x)

    def _get_as_int_doctest(self):
        r"""
        This method exists solely for automated testing of get_as_int().

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: x = ntl.ZZ_p(42,modulus=c)
            sage: i = x._get_as_int_doctest()
            sage: i
            2
            sage: type(i)
            <... 'int'>
        """
        self.c.restore_c()
        return self.get_as_int()

    cdef void set_from_int(ntl_ZZ_p self, int value) noexcept:
        r"""
        Set the value from a C int.

        AUTHOR: David Harvey (2006-08-05)
        """
        self.c.restore_c()
        self.x = int_to_ZZ_p(value)

    def _set_from_int_doctest(self, value):
        r"""
        This method exists solely for automated testing of set_from_int().

        EXAMPLES::

            sage: c=ntl.ZZ_pContext(ntl.ZZ(20))
            sage: x = ntl.ZZ_p(modulus=c)
            sage: x._set_from_int_doctest(42)
            sage: x
            2
            sage: x = ntl.ZZ_p(7,81)
            sage: x._set_from_int_doctest(int(3))
            sage: x
            3
        """
        self.c.restore_c()
        self.set_from_int(int(value))

    def lift(self):
        """
        Return a lift of ``self`` as an ntl.ZZ object.

        EXAMPLES::

            sage: x = ntl.ZZ_p(8,18)
            sage: x.lift()
            8
            sage: type(x.lift())
            <class 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
        """
        cdef ntl_ZZ r = ntl_ZZ()
        self.c.restore_c()
        r.x = ZZ_p_rep(self.x)
        return r

    def modulus(self):
        r"""
        Return the modulus as an NTL ZZ.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(ntl.ZZ(20))
            sage: n = ntl.ZZ_p(2983,c)
            sage: n.modulus()
            20
        """
        cdef ntl_ZZ r = ntl_ZZ()
        self.c.restore_c()
        ZZ_p_modulus( &r.x, &self.x )
        return r

    def lift_centered(self):
        """
        Compute a representative of ``self`` in `(-n/2 , n/2]` as an
        ``ntl.ZZ`` object.

        OUTPUT:

        - A ``ntl.ZZ`` object `r` such that  `-n/2 < r \\leq n/2` and `Mod(r, n) == self`.

        EXAMPLES::

            sage: x = ntl.ZZ_p(8, 18)
            sage: x.lift_centered()
            8
            sage: type(x.lift_centered())
            <class 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
            sage: x = ntl.ZZ_p(12, 18)
            sage: x.lift_centered()
            -6
            sage: type(x.lift_centered())
            <class 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
        """
        cdef ntl_ZZ r = self.lift()
        cdef ntl_ZZ m = self.modulus()
        if r*2 > m:
            r -= m
        return r

    def _integer_(self, ZZ=None):
        """
        Return a lift of ``self`` as a Sage integer.

        EXAMPLES::

            sage: x = ntl.ZZ_p(8,188)
            sage: x._integer_()
            8

            sage: type(x._integer_())
            <class 'sage.rings.integer.Integer'>
        """
        self.c.restore_c()
        cdef ZZ_c rep = ZZ_p_rep(self.x)
        cdef Integer ans = Integer.__new__(Integer)
        ZZ_to_mpz(ans.value, &rep)
        return ans

    def _sage_(self):
        r"""
        Return the value as a sage IntegerModRing.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(20)
            sage: n = ntl.ZZ_p(2983, c)
            sage: type(n._sage_())
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: n
            3

        AUTHOR: Joel B. Mohler
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

        cdef ZZ_c rep
        self.c.restore_c()
        rep = ZZ_p_rep(self.x)
        cdef Integer ans = Integer.__new__(Integer)
        ZZ_to_mpz(ans.value, &rep)
        return IntegerModRing(self.modulus()._integer_())(ans)
