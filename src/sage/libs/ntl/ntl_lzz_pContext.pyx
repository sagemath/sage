# distutils: libraries = NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.mpz cimport mpz_get_si
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer cimport Integer

zz_pContextDict = {}

cdef class ntl_zz_pContext_class():
    def __init__(self, long v):
        """
        EXAMPLES::

            # You can construct contexts manually.
            sage: c = ntl.zz_pContext(11)
            sage: n1 = ntl.zz_p(12,c)
            sage: n1
            1

            # or You can construct contexts implicitly.
            sage: n2=ntl.zz_p(12, 7)
            sage: n2
            5
            sage: ntl.zz_p(2,3)+ntl.zz_p(1,3)
            0
            sage: n2+n1  # Mismatched moduli:  It will go BOOM!
            Traceback (most recent call last):
            ...
            ValueError: arithmetic operands must have the same modulus.
        """
        pass

    def __cinit__(self, long v):
        if v > NTL_SP_BOUND:
            raise ValueError("Modulus (=%s) is too big" % v)
        elif v < 2:
            # Issue 13940: only moduli greater than one are supported.
            raise ValueError("Modulus (=%s) is too small" % v)

        self.x = zz_pContext_c(v)
        zz_pContextDict[repr(v)] = self
        self.p = v

    def __reduce__(self):
        """
        EXAMPLES::

            sage: c=ntl.zz_pContext(13)
            sage: loads(dumps(c)) is c
            True
        """
        return ntl_zz_pContext, (self.p,)

    def modulus(self):
        """
        Print the modulus for ``self``.

        EXAMPLES::

            sage: c1 = ntl.zz_pContext(36)
            sage: c1.modulus()
            36
        """
        return self.p

    def restore(self):
        """
        Restore a zz_pContext.

        EXAMPLES::

            sage: c = ntl.zz_pContext(5)
            sage: m = ntl.zz_p(4,7)
            sage: c.restore()
        """
        self.restore_c()

    cdef void restore_c(self) noexcept:
        """
        Actual code for the above.

        EXAMPLES::

            sage: n = ntl.zz_p(3,5)
            sage: m = ntl.zz_p(4,7)
            sage: n*n ## indirect doctest
            4
        """
        self.x.restore()


def ntl_zz_pContext( v ):
    """
    Creation function for a zz_p context.

    EXAMPLES::

        sage: f = ntl.zz_pContext(26)
        sage: f = ntl.zz_pContext(10^100)
        Traceback (most recent call last):
        ...
        ValueError: Modulus (=10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000) is too big
        sage: f = ntl.zz_pContext(1)
        Traceback (most recent call last):
        ...
        ValueError: Modulus (=1) is too small
        sage: f = ntl.zz_pContext(0)
        Traceback (most recent call last):
        ...
        ValueError: Modulus (=0) is too small
    """
    if v > NTL_SP_BOUND:
        raise ValueError("Modulus (=%s) is too big" % v)
    if isinstance(v, Integer):
        v = mpz_get_si((<Integer>v).value)
    try:
        return zz_pContextDict[repr(v)]
    except KeyError:
        return ntl_zz_pContext_class(v)
