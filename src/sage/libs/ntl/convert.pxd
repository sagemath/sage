from sage.libs.ntl.types cimport ZZ_c
from sage.libs.gmp.types cimport mpz_t, mpz_srcptr

cdef void ZZ_to_mpz(mpz_t output, ZZ_c* x) noexcept
cdef void mpz_to_ZZ(ZZ_c *output, mpz_srcptr x) noexcept
cdef void PyLong_to_ZZ(ZZ_c* z, value) noexcept
