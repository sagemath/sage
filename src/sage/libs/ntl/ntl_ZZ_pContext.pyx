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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'misc.pxi'
include 'decl.pxi'
import weakref

from sage.ext.cplusplus cimport ccrepr
from sage.rings.integer cimport Integer


cdef class ntl_ZZ_pContext_class():
    def __init__(self, ntl_ZZ v):
        """
        EXAMPLES::

            # You can construct contexts manually.
            sage: c = ntl.ZZ_pContext(11)
            sage: n1 = ntl.ZZ_p(12,c)
            sage: n1
            1

            # or You can construct contexts implicitly.
            sage: n2 = ntl.ZZ_p(12, 7)
            sage: n2
            5
            sage: ntl.ZZ_p(2,3)+ntl.ZZ_p(1,3)
            0
            sage: n2+n1  # Mismatched moduli:  It will go BOOM!
            Traceback (most recent call last):
            ...
            ValueError: You cannot perform arithmetic with elements of different moduli.
        """
        pass

    def __cinit__(self, ntl_ZZ v):
        self.x = ZZ_pContext_c(v.x)
        self.p = v
        self.p_bits = self.p._integer_().nbits()

    def __reduce__(self):
        """
        EXAMPLES::

            sage: c = ntl.ZZ_pContext(13)
            sage: loads(dumps(c)) is c
            True
        """
        return ntl_ZZ_pContext, (self.p,)

    def __repr__(self):
        """
        Return a print representation of ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(7)
            sage: c
            NTL modulus 7
        """
        return "NTL modulus %s"%(self.p)

    def __hash__(self):
        return hash(self.p)

    def modulus(self):
        """
        Return the current modulus associated to this
        context.

        EXAMPLES::

            sage: c = ntl.ZZ_pContext(7)
            sage: c.modulus()
            7

            sage: c = ntl.ZZ_pContext(10^30)
            sage: type(c.modulus())
            <class 'sage.rings.integer.Integer'>
            sage: c.modulus() == 10^30
            True
        """
        return Integer(self.p)

    def restore(self):
        """
        EXAMPLES::

            sage: c1 = ntl.ZZ_p(5,92) ; c2 = ntl.ZZ_p(7,92)
            sage: c1+c2
            12
            sage: d1 = ntl.ZZ_p(38,91) ; d2 = ntl.ZZ_p(3,91)
            sage: d1*d2 ## indirect doctest
            23
        """
        self.restore_c()

    cdef void restore_c(self) noexcept:
        self.x.restore()

    cpdef void _assert_is_current_modulus(self) except *:
        """
        Assert that is currently-set NTL modulus.

        Mostly for debugging purposes. If false, an assertion is raised. This method segfaults if
        the NTL modulus has never been set before.

        EXAMPLES::

            sage: c1 = ntl.ZZ_pContext(7)
            sage: c2 = ntl.ZZ_pContext(5)
            sage: c1.restore()
            sage: c1._assert_is_current_modulus()
            sage: c2._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: 5 != 7
            sage: c2.restore()
            sage: c1._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: 7 != 5
            sage: c2._assert_is_current_modulus()
            sage: ntl.ZZ_pContext(3).restore()
            sage: c1._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: 7 != 3
            sage: c2._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: 5 != 3
        """
        if self.p.x == ntl_ZZ_p_current_modulus():
            return
        raise AssertionError('modulus mismatch: {} != {}'.format(
            self.p,
            ccrepr(ntl_ZZ_p_current_modulus())))


cdef class ntl_ZZ_pContext_factory():

    def __init__(self):
        self.context_dict = {}

    cdef ntl_ZZ_pContext_class make_c(self, ntl_ZZ v):
        """
        Create a new ZZ_pContext.

        INPUT:

        - ``v`` -- an ntl_ZZ
        """
        cdef ntl_ZZ_pContext_class context
        if v in self.context_dict:
            context = <ntl_ZZ_pContext_class> self.context_dict[v]()
            if context is not None:
                return context
        context = ntl_ZZ_pContext_class(v)
        self.context_dict[v] = weakref.ref(context)
        return context


ZZ_pContext_factory = ntl_ZZ_pContext_factory()


def ntl_ZZ_pContext( v ):
    """
    Create a new ZZ_pContext.

    EXAMPLES::

        sage: c = ntl.ZZ_pContext(178)
        sage: n1 = ntl.ZZ_p(212,c)
        sage: n1
        34
    """
    v = ntl_ZZ(v)
    if (v < ntl_ZZ(2)):
        raise ValueError("%s is not a valid modulus." % v)
    return (<ntl_ZZ_pContext_factory>ZZ_pContext_factory).make_c(v)
