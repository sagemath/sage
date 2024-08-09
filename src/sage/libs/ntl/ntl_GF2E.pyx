# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
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

from sage.ext.cplusplus cimport ccrepr

include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_GF2 cimport ntl_GF2
from sage.libs.ntl.ntl_GF2X cimport ntl_GF2X
from sage.libs.ntl.ntl_GF2EContext cimport ntl_GF2EContext_class
from sage.libs.ntl.ntl_GF2EContext import ntl_GF2EContext
from sage.libs.ntl.ntl_ZZ import unpickle_class_args
from sage.misc.randstate cimport current_randstate


##############################################################################
#
# ntl_GF2E: GF(2**x) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on code by William Stein)
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2007-10: adapted to new conventions
#
##############################################################################

def ntl_GF2E_random(ntl_GF2EContext_class ctx):
    """
    Return a random element from GF2E modulo the current modulus.

    INPUT:

    - ``ctx`` -- the GF2E context for which a random element should be created

    EXAMPLES::

        sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
        sage: ntl.GF2E_random(ctx).modulus_context() is ctx
        True
    """
    current_randstate().set_seed_ntl(False)

    cdef ntl_GF2E r
    ctx.restore_c()
    r = ntl_GF2E.__new__(ntl_GF2E)
    r.c = ctx
    r.x = GF2E_random()
    return r


cdef class ntl_GF2E():
    r"""
    The :class:`GF2E` represents a finite extension field over GF(2)
    using NTL. Elements are represented as polynomials over GF(2)
    modulo a modulus.

    This modulus must be set by creating a GF2EContext first and pass
    that context to the constructor of all elements.
    """

    def __init__(self, x=None, modulus=None):
        """
        Construct a new finite field element in GF(2**x).

        If you pass a string to the constructor please note that byte
        sequences and the hexadecimal notation are Little Endian in
        NTL.  So e.g. '[0 1]' == '0x2' == x.

        INPUT:

        - ``x`` -- value to be assigned to this element. Same types as
          ``ntl.GF2X()`` are accepted.
        - ``modulus`` -- the context/modulus of the field

        OUTPUT: a new ``ntl.GF2E`` element

        EXAMPLES::

            sage: k.<a> = GF(2^8)
            sage: e = ntl.GF2E(a,k); e
            [0 1]
            sage: ctx = e.modulus_context()
            sage: ntl.GF2E('0x1c', ctx)
            [1 0 0 0 0 0 1 1]
            sage: ntl.GF2E('[1 0 1 0]', ctx)
            [1 0 1]
            sage: ntl.GF2E([1,0,1,0], ctx)
            [1 0 1]
            sage: ntl.GF2E(ntl.GF2(1),ctx)
            [1]
        """
        if modulus is None:
            raise ValueError("You must specify a modulus when creating a GF2E.")

        cdef ntl_GF2X _x

        if isinstance(x, ntl_GF2X):
            GF2E_conv_GF2X(self.x, (<ntl_GF2X>x).x)

        elif isinstance(x, int):
            GF2E_conv_long(self.x, x)

        elif isinstance(x, ntl_ZZ):
            GF2E_conv_ZZ(self.x, (<ntl_ZZ>x).x)

        elif isinstance(x, ntl_GF2):
            GF2E_conv_GF2(self.x, (<ntl_GF2>x).x)
        else:
            _x = ntl_GF2X(x)
            GF2E_conv_GF2X(self.x, _x.x)

    def __cinit__(self, x=None, modulus=None):
        #################### WARNING ###################
        ## Before creating a GF2E, you must create a  ##
        ## GF2EContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing an ntl_GF2E          ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = ntl_GF2E.__new__(ntl_GF2E) without
        ## first restoring a GF2EContext, which could ##
        ## have unfortunate consequences.  See _new   ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None:
            return
        if isinstance(modulus, ntl_GF2EContext_class):
            self.c = <ntl_GF2EContext_class>modulus
            self.c.restore_c()
        else:
            self.c = <ntl_GF2EContext_class>ntl_GF2EContext(modulus)
            self.c.restore_c()

    cdef ntl_GF2E _new(self):
        cdef ntl_GF2E r
        self.c.restore_c()
        r = ntl_GF2E.__new__(ntl_GF2E)
        r.c = self.c
        return r

    def __reduce__(self):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext( ntl.GF2X([1,1,0,1,1,0,0,0,1]) )
            sage: a = ntl.GF2E(ntl.ZZ_pX([1,1,3],2), ctx)
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_class_args, (ntl_GF2E, (str(self), self.modulus_context()))

    def modulus_context(self):
        """
        Return the structure that holds the underlying NTL GF2E modulus.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext( ntl.GF2X([1,1,0,1,1,0,0,0,1]) )
            sage: a = ntl.GF2E(ntl.ZZ_pX([1,1,3],2), ctx)
            sage: cty = a.modulus_context(); cty
            NTL modulus [1 1 0 1 1 0 0 0 1]
            sage: ctx == cty
            True
        """
        return self.c

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: ntl.GF2E([1,1,0,1], ctx) # indirect doctest
            [1 1 0 1]
        """
        self.c.restore_c()
        return ccrepr(self.x)

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: x = ntl.GF2E([0,1,1],GF(2^4,'a'))
            sage: y = copy(x)
            sage: x == y
            True
            sage: x is y
            False
        """
        cdef ntl_GF2E r = self._new()
        r.x = self.x
        return r

    def __mul__(ntl_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([1,1,0,1,1], ctx)
            sage: x*y ## indirect doctest
            [0 0 1 1 1 0 1 1]
        """
        cdef ntl_GF2E r
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other,self.c)
        elif self.c is not (<ntl_GF2E>other).c:
            raise ValueError("You cannot perform arithmetic with elements in different fields.")
        r = self._new()
        GF2E_mul(r.x, self.x, (<ntl_GF2E>other).x)
        return r

    def __sub__(ntl_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([1,1,0,1,1], ctx)
            sage: x - y ## indirect doctest
            [0 1 1 1]
        """
        cdef ntl_GF2E r
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other,self.c)
        elif self.c is not (<ntl_GF2E>other).c:
            raise ValueError("You cannot perform arithmetic with elements in different fields.")
        r = self._new()
        GF2E_sub(r.x, self.x, (<ntl_GF2E>other).x)
        return r

    def __add__(ntl_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([1,1,0,1,1], ctx)
            sage: x+y ## indirect doctest
            [0 1 1 1]
        """
        cdef ntl_GF2E r
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other,self.c)
        elif self.c is not (<ntl_GF2E>other).c:
            raise ValueError("You cannot perform arithmetic with elements in different fields.")
        r = self._new()
        GF2E_add(r.x, self.x, (<ntl_GF2E>other).x)
        return r

    def __truediv__(ntl_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([1,1,0,1,1], ctx)
            sage: x/y ## indirect doctest
            [1 0 1 0 0 1 0 1]
        """
        cdef ntl_GF2E r
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other, self.c)
        elif self.c is not (<ntl_GF2E>other).c:
            raise ValueError("You cannot perform arithmetic with elements in different fields.")
        r = self._new()
        GF2E_div(r.x, self.x, (<ntl_GF2E>other).x)
        return r

    def __neg__(ntl_GF2E self):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx)
            sage: -x ## indirect doctest
            [1 0 1 0 1]
        """
        cdef ntl_GF2E r = self._new()
        r.x = self.x
        return r

    def __pow__(ntl_GF2E self, long e, ignored):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx)
            sage: x**2 ## indirect doctest
            [0 1 0 1]
        """
        cdef ntl_GF2E r = self._new()
        GF2E_power(r.x, self.x, e)
        return r

    def __richcmp__(ntl_GF2E self, other, int op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([1,1,0,1,1], ctx)
            sage: x == x
            True
            sage: x == y
            False
            sage: ntl.GF2E(0,ctx) == 0
            True
            sage: a = ntl.GF2E([0,1],GF(2^2,'a'))
            sage: a == x
            False
        """
        self.c.restore_c()

        if op != Py_EQ and op != Py_NE:
            raise TypeError("elements of GF(2^e) are not ordered")

        cdef ntl_GF2E b
        try:
            b = <ntl_GF2E?>other
        except TypeError:
            b = ntl_GF2E(other, self.c)

        return (op == Py_EQ) == (self.x == b.x)

    def IsZero(ntl_GF2E self):
        """
        Return ``True`` if this element equals zero, ``False`` otherwise.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([1,1,0,1,1,0,0,0,1], ctx)
            sage: x.IsZero()
            False
            sage: y.IsZero()
            True
        """
        return bool(GF2E_IsZero(self.x))

    def IsOne(ntl_GF2E self):
        """
        Return ``True`` if this element equals one, ``False`` otherwise.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([0,1,0,1,1,0,0,0,1], ctx)
            sage: x.IsOne()
            False
            sage: y.IsOne()
            True
        """
        return bool(GF2E_IsOne(self.x))

    def trace(ntl_GF2E self):
        """
        Return the trace of this element.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: x = ntl.GF2E([1,0,1,0,1], ctx) ; y = ntl.GF2E([0,1,1,0,1,1], ctx)
            sage: x.trace()
            0
            sage: y.trace()
            1
        """
        cdef ntl_GF2 x = ntl_GF2.__new__(ntl_GF2)
        x.x = GF2E_trace(self.x)
        return x

    def rep(ntl_GF2E self):
        """
        Return a ntl.GF2X copy of this element.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: a = ntl.GF2E('0x1c', ctx)
            sage: a.rep()
            [1 0 0 0 0 0 1 1]
            sage: type(a.rep())
            <class 'sage.libs.ntl.ntl_GF2X.ntl_GF2X'>
        """
        cdef ntl_GF2X x = ntl_GF2X.__new__(ntl_GF2X)
        x.x = GF2E_rep(self.x)
        return x

    def list(ntl_GF2E self):
        """
        Represent this element as a list of binary digits.

        EXAMPLES::

             sage: e=ntl.GF2E([0,1,1],GF(2^4,'a'))
             sage: e.list()
             [0, 1, 1]
             sage: e=ntl.GF2E('0xff',GF(2^8,'a'))
             sage: e.list()
             [1, 1, 1, 1, 1, 1, 1, 1]

        OUTPUT: list of digits representing the coefficients in this element's
        polynomial representation
        """
        cdef int i
        cdef GF2X_c x = GF2E_rep(self.x)
        cdef ntl_GF2 b

        l = []

        for i from 0 <= i <= GF2X_deg(x):
            b = ntl_GF2.__new__(ntl_GF2)
            b.x = GF2X_coeff(x,i)
            l.append(b)
        return l

    def _sage_(ntl_GF2E self, k=None):
        r"""
        Return a :class:`FiniteFieldElement` representation of this element.

        If a :class:`FiniteField` `k` is provided, it is constructed
        in this field if possible. A :class:`FiniteField` will be
        constructed if none is provided.

        INPUT:

        - ``k`` -- (optional) a field `\GF{2^d}`

        OUTPUT: :class:`FiniteFieldElement` over `k`

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: e = ntl.GF2E([0,1], ctx)
            sage: a = e._sage_(); a
            a
        """
        cdef int i
        cdef int length

        self.c.restore_c()

        e = GF2E_degree()

        if k is None:
            from sage.rings.finite_rings.finite_field_constructor import FiniteField
            f = self.c.m._sage_()
            k = FiniteField(2**e, name='a', modulus=f)

        a=k.gen()
        l = self.list()

        length = len(l)
        ret = 0

        for i from 0 <= i < length:
            if l[i]==1:
                ret = ret + a**i

        return ret
