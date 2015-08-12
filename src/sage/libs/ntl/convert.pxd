# distutils: depends = NTL/ZZ.h

# Copied from ntl_ZZ_decl to avoid circular cimports
cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    ctypedef struct ZZ_c "struct ZZ":
        pass

from sage.libs.gmp.types cimport mpz_t, mpz_srcptr

cdef void ZZ_to_mpz(mpz_t output, ZZ_c* x)
cdef void mpz_to_ZZ(ZZ_c *output, mpz_srcptr x)
cdef void PyLong_to_ZZ(ZZ_c* z, value)
