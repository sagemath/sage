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

include 'misc.pxi'
include 'decl.pxi'

from sage.ext.cplusplus cimport ccrepr
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX


ZZ_pEContextDict = {}


cdef class ntl_ZZ_pEContext_class():
    def __init__(self, ntl_ZZ_pX f):
        """
        EXAMPLES:

        You can construct contexts manually::

            sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([4,1,6],25))
            sage: n1=c.ZZ_pE([10,17,12])
            sage: n1
            [2 15]

        Or you can construct contexts implicitly::

            sage: n2=ntl.ZZ_pE(12, ntl.ZZ_pX([1,1,1],7))
            sage: n2
            [5]
            sage: n2+n1  # Mismatched moduli:  It will go BOOM!
            Traceback (most recent call last):
            ...
            ValueError: You cannot perform arithmetic with elements of different moduli.
        """
        pass

    def __cinit__(self, ntl_ZZ_pX f):
        self.pc = f.c
        self.pc.restore_c()
        self.x = ZZ_pEContext_c(f.x)
        ZZ_pEContextDict[(repr(f), repr(f.c.p))] = self
        self.f = f
        self.ptrs.zzpc = &(self.pc.x)
        self.ptrs.zzpec = &(self.x)

    def __reduce__(self):
        """
        EXAMPLES::

           sage: c=ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1],7))
           sage: loads(dumps(c)) is c
           True
        """
        return ntl_ZZ_pEContext, (self.f,)

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7)); c
            NTL modulus [1 1 1] (mod 7)
        """
        return "NTL modulus %s (mod %s)"%(self.f, self.pc.p)

    def get_pc(self):
        """
        Return the ZZ_pContext contained within ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7)); c
            NTL modulus [1 1 1] (mod 7)
            sage: c.get_pc()
            NTL modulus 7
        """
        return self.pc

    def polynomial(self):
        """
        Return the ZZ_pX polynomial defining ``self``.

        EXAMPLES::

            sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
            sage: c.polynomial()
            [1 1 1]
        """
        return self.f

    def restore(self):
        """
        Manually set the global NTL modulus to be ``self``.

        This should be done automatically by all of the NTL wrapper classes.

        CRUCIAL: If you are writing your own classes that use ZZ_p_c, ZZ_pX_c, ZZ_pE_c, ZZ_pEX_c
        then you MUST restore the context before calling off to NTL for anything.  If the context has been
        switched by other code then behavior of operations is undefined.  See the NTL documentation for
        more details (or the wrappers in sage.libs.ntl)
        """
        self.restore_c()

    cdef void restore_c(self) noexcept:
        """
        Set the global NTL modulus to be ``self``.

        CRUCIAL: If you are writing your own classes that use ZZ_p_c, ZZ_pX_c, ZZ_pE_c, ZZ_pEX_c
        then you MUST restore the context before calling off to NTL for anything.  If the context has been
        switched by other code then behavior of operations is undefined.  See the NTL documentation for
        more details (or the wrappers in sage.libs.ntl)
        """
        self.pc.restore_c()
        self.x.restore()

    #def ZZ_pX(self,v = None):
    #    from ntl_ZZ_pX import ntl_ZZ_pX
    #    return ntl_ZZ_pX(v,modulus=self)

    def ZZ_pE(self, v = None):
        """
        Return a ZZ_pE object with modulus ``self`` out of the data v.

        EXAMPLES::

            sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
            sage: c.ZZ_pE([4,3])
            [4 3]
        """
        from sage.libs.ntl.ntl_ZZ_pE import ntl_ZZ_pE
        return ntl_ZZ_pE(v,modulus=self)

    def ZZ_pEX(self, v = None):
        """
        Return a ZZ_pE object with modulus ``self`` out of the data v.

        EXAMPLES::

            sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7))
            sage: c.ZZ_pEX([4,3])
            [[4] [3]]
        """
        from sage.libs.ntl.ntl_ZZ_pEX import ntl_ZZ_pEX
        return ntl_ZZ_pEX(v, modulus=self)

    cpdef void _assert_is_current_modulus(self) except *:
        """
        Assert that this is currently-set NTL modulus.

        Mostly for debugging purposes. If false, an assertion is raised. This method
        segfaults if the NTL modulus has never been set before.

        EXAMPLES::

            sage: c1 = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 5))
            sage: c2 = ntl.ZZ_pEContext(ntl.ZZ_pX([1,2,1], 5))
            sage: c1.restore()
            sage: c1._assert_is_current_modulus()
            sage: c2._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: [1 2 1] != [1 1 1]
            sage: c2.restore()
            sage: c1._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: [1 1 1] != [1 2 1]
            sage: c2._assert_is_current_modulus()
            sage: ntl.ZZ_pContext(3).restore()
            sage: c1._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: 5 != 3
            sage: c2._assert_is_current_modulus()
            Traceback (most recent call last):
            ...
            AssertionError: modulus mismatch: 5 != 3
        """
        self.pc._assert_is_current_modulus()
        if self.f.x == ZZ_pE_current_modulus().val():
            return
        raise AssertionError('modulus mismatch: {} != {}'.format(
            ccrepr(self.f.x),
            ccrepr(ZZ_pE_current_modulus().val())))


def ntl_ZZ_pEContext( ntl_ZZ_pX f):
    """
    Create an ntl_ZZ_pEContext.

    Such an object must be created before any ZZ_pE or ZZ_pEX objects can be used.

    The context handling should be taken care of by the wrapper classes.

    EXAMPLES::

        sage: c = ntl.ZZ_pEContext(ntl.ZZ_pX([1,1,1], 7)); c
        NTL modulus [1 1 1] (mod 7)
    """
    try:
        return ZZ_pEContextDict[repr(f), repr(f.c.p)]
    except KeyError:
        # Creating the following object caches it.
        return ntl_ZZ_pEContext_class(f)
