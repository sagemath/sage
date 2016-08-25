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

from sage.rings.integer_ring import IntegerRing

ZZ_sage = IntegerRing()


cdef class ntl_ZZ_pContext_class(object):
    def __init__(self, ntl_ZZ v):
        """
        EXAMPLES:
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
            ValueError: You can not perform arithmetic with elements of different moduli.
        """
        pass

    def __cinit__(self, ntl_ZZ v):
        self.x = ZZ_pContext_c(v.x)
        self.p = v
        self.p_bits = self.p._integer_().nbits()

    def __reduce__(self):
        """
        EXAMPLES:
            sage: c = ntl.ZZ_pContext(13)
            sage: loads(dumps(c)) is c
            True
        """
        return ntl_ZZ_pContext, (self.p,)

    def __repr__(self):
        """
        Returns a print representation of self.

        EXAMPLES:
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

        EXAMPLES:
            sage: c = ntl.ZZ_pContext(7)
            sage: c.modulus()
            7

            sage: c = ntl.ZZ_pContext(10^30)
            sage: type(c.modulus())
            <type 'sage.rings.integer.Integer'>
            sage: c.modulus() == 10^30
            True
        """
        return ZZ_sage(self.p)


    def restore(self):
        """
        EXAMPLES:
            sage: c1 = ntl.ZZ_p(5,92) ; c2 = ntl.ZZ_p(7,92)
            sage: c1+c2
            12
            sage: d1 = ntl.ZZ_p(38,91) ; d2 = ntl.ZZ_p(3,91)
            sage: d1*d2 ## indirect doctest
            23
        """
        self.restore_c()

    cdef void restore_c(self):
        self.x.restore()

cdef class ntl_ZZ_pContext_factory(object):
    def __init__(self):
        self.context_dict = {}

    cdef ntl_ZZ_pContext_class make_c(self, ntl_ZZ v):
        """
        Creates a new ZZ_pContext.

        INPUT:
        v -- an ntl_ZZ
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
    EXAMPLES:
        sage: c = ntl.ZZ_pContext(178)
        sage: n1 = ntl.ZZ_p(212,c)
        sage: n1
        34
    """
    v = ntl_ZZ(v)
    if (v < ntl_ZZ(2)):
        raise ValueError("%s is not a valid modulus." % v)
    return (<ntl_ZZ_pContext_factory>ZZ_pContext_factory).make_c(v)
