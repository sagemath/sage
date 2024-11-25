# distutils: libraries = NTL_LIBRARIES gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

# ****************************************************************************
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

from sage.ext.cplusplus cimport ccrepr, ccreadstr

include 'misc.pxi'
include 'decl.pxi'

from cpython.object cimport Py_EQ, Py_NE
from sage.rings.integer cimport Integer

##############################################################################
# GF2: Bits
##############################################################################

cdef class ntl_GF2():
    r"""
    The \class{GF2} represents the field GF(2). Computationally
    speaking, it is not a particularly useful class.  Its main use is
    to make the interfaces to the various finite field classes as
    uniform as possible.
    """
    def __init__(self, v=None):
        r"""
        Initialize a NTL bit.

        EXAMPLES::

            sage: ntl.GF2(1)
            1
            sage: ntl.GF2(int(2))
            0
            sage: ntl.GF2('1')
            1
        """
        if isinstance(v, ntl_GF2):
            self.x = (<ntl_GF2>v).x
        elif isinstance(v, (int, Integer)):
            GF2_conv_long(self.x, int(v) % 2)
        elif v is not None:
            ccreadstr(self.x, str(v))

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: str(ntl.GF2(1)) # indirect doctest
            '1'
        """
        return ccrepr(self.x)

    def __reduce__(self):
        """
        Serialize ``self``.

        EXAMPLES::

            sage: a = ntl.GF2(1)
            sage: loads(dumps(a))
            1
        """
        return unpickle_class_value, (ntl_GF2, int(self))

    def __richcmp__(ntl_GF2 self, other, int op):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: a = ntl.GF2(1)
            sage: b = ntl.GF2(0)
            sage: a == a
            True
            sage: a == b
            False
            sage: a == 1
            True
            sage: a < b
            Traceback (most recent call last):
            ...
            TypeError: elements of GF(2) are not ordered
        """
        if op != Py_EQ and op != Py_NE:
            raise TypeError("elements of GF(2) are not ordered")

        cdef ntl_GF2 b
        try:
            b = <ntl_GF2?>other
        except TypeError:
            b = ntl_GF2(other)

        return (op == Py_EQ) == (self.x == b.x)

    def __mul__(self, other):
        """
        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o*o
            1
            sage: o*z
            0
            sage: z*o
            0
            sage: z*z
            0
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        GF2_mul(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __truediv__(self, other):
        """
        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o/o
            1
            sage: o/z
            Traceback (most recent call last):
            ...
            ZeroDivisionError
        """
        cdef ntl_GF2 r
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        if GF2_IsZero((<ntl_GF2>other).x):
            raise ZeroDivisionError
        r = ntl_GF2.__new__(ntl_GF2)
        GF2_div(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __sub__(self, other):
        """
        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o-o
            0
            sage: o-z
            1
            sage: z-o
            1
            sage: z-z
            0
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        GF2_sub(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __add__(self, other):
        """
        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: o+o
            0
            sage: o+z
            1
            sage: z+o
            1
            sage: z+z
            0
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        if not isinstance(self, ntl_GF2):
            self = ntl_GF2(self)
        if not isinstance(other, ntl_GF2):
            other = ntl_GF2(other)
        GF2_add(r.x, (<ntl_GF2>self).x, (<ntl_GF2>other).x)
        return r

    def __neg__(ntl_GF2 self):
        """
        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: -z
            0
            sage: -o
            1
        """
        cdef ntl_GF2 r = ntl_GF2.__new__(ntl_GF2)
        GF2_negate(r.x, self.x)
        return r

    def __pow__(ntl_GF2 self, long e, ignored):
        """
        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: z^2
            0
            sage: o^2
            1
        """
        cdef ntl_GF2 r = ntl_GF2()
        GF2_power(r.x, self.x, e)
        return r

    def __int__(self):
        """
        Return ``self`` as an int.

        EXAMPLES::

            sage: o = ntl.GF2(1)
            sage: z = ntl.GF2(0)
            sage: int(z)
            0
            sage: int(o)
            1
        """
        cdef long l = GF2_conv_to_long(self.x)
        return int(l)


def unpickle_class_value(cls, x):
    """
    Here for unpickling.

    EXAMPLES::

        sage: sage.libs.ntl.ntl_GF2.unpickle_class_value(ntl.GF2,1)
        1
        sage: type(sage.libs.ntl.ntl_GF2.unpickle_class_value(ntl.GF2,1))
        <class 'sage.libs.ntl.ntl_GF2.ntl_GF2'>
    """
    return cls(x)


def unpickle_class_args(cls, x):
    """
    Here for unpickling.

    EXAMPLES::

        sage: sage.libs.ntl.ntl_GF2.unpickle_class_args(ntl.GF2,[1])
        1
        sage: type(sage.libs.ntl.ntl_GF2.unpickle_class_args(ntl.GF2,[1]))
        <class 'sage.libs.ntl.ntl_GF2.ntl_GF2'>
    """
    return cls(*x)
