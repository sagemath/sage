# distutils: depends = NTL/ZZ.h
# distutils: libraries = NTL_LIBRARIES gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
Conversion between NTL's ``ZZ`` and various other types
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.mpz cimport mpz_init, mpz_clear
from sage.libs.gmp.pylong cimport mpz_set_pylong

cdef extern from "ntlwrap_impl.h":
    void ZZ_to_mpz(mpz_t output, ZZ_c* x)
    void mpz_to_ZZ(ZZ_c *output, mpz_srcptr x)

cdef void PyLong_to_ZZ(ZZ_c* z, value) noexcept:
    """
    Convert ``value`` (which must be a Python ``long``) to NTL.
    """
    cdef mpz_t t
    mpz_init(t)
    mpz_set_pylong(t, value)
    mpz_to_ZZ(z, t)
    mpz_clear(t)
