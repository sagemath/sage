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

GF2EContextDict = {}


cdef class ntl_GF2EContext_class():
    def __init__(self, ntl_GF2X v):
        """
        EXAMPLES::

            # You can construct contexts manually.
            sage: ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1]))
            sage: n1 = ntl.GF2E([1,1],ctx)
            sage: n1
            [1 1]

            # or You can construct contexts implicitly.
            sage: n2 = ntl.GF2E([0,1], ntl.GF2X([1,1,0,1]))
            sage: n2
            [0 1]
            sage: ntl.GF2E(2, GF(2^8,'a'))+ntl.GF2E([0,1],ctx)
            Traceback (most recent call last):
            ...
            ValueError: You cannot perform arithmetic with elements in different fields.

            sage: n2+n1  # Mismatched moduli:  It will go BOOM!
            [1]
        """
        pass

    def __cinit__(self, ntl_GF2X v):
        self.x = GF2EContext_c(v.x)
        self.m = v

    def __reduce__(self):
        """
        EXAMPLES::

            sage: c = ntl.GF2EContext(GF(2^5,'b'))
            sage: loads(dumps(c)) is c
            True
        """
        return ntl_GF2EContext, (self.m,)

    def __repr__(self):
        """
        Return a print representation of ``self``.

        EXAMPLES::

            sage: c = ntl.GF2EContext(GF(2^16,'a'))
            sage: c
            NTL modulus [1 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 1]
        """
        return "NTL modulus %s"%(self.m)

    def modulus(self):
        """
        Return the current modulus associated to this
        context.

        EXAMPLES::

            sage: c = ntl.GF2EContext(GF(2^7,'foo'))
            sage: c.modulus()
            [1 1 0 0 0 0 0 1]
        """
        return self.m

    def restore(self):
        """
        EXAMPLES::

            sage: c1 = ntl.GF2E([0,1],GF(2^4,'a')) ; c2 = ntl.GF2E([1,0,1],GF(2^4,'a'))
            sage: c1+c2
            [1 1 1]
            sage: d1 = ntl.GF2E([0,1],GF(2^5,'a')) ; d2 = ntl.GF2E([0,0,1],GF(2^5,'a'))
            sage: d1*d2 ## indirect doctest
            [0 0 0 1]
        """
        self.restore_c()

    cdef void restore_c(self) noexcept:
        self.x.restore()


def ntl_GF2EContext( v ):
    """
    Create a new GF2EContext.

    EXAMPLES::

        sage: c = ntl.GF2EContext(GF(2^2,'a'))
        sage: n1 = ntl.GF2E([0,1],c)
        sage: n1
        [0 1]
    """
    v = ntl_GF2X(v)
    if (GF2X_deg((<ntl_GF2X>v).x) < 1):
        raise ValueError("%s is not a valid modulus." % v)
    key = hash(v)
    if key in GF2EContextDict:
        context = GF2EContextDict[key]()
        if context is not None:
            return context
    context = ntl_GF2EContext_class(v)
    GF2EContextDict[key] = weakref.ref(context)
    return context
