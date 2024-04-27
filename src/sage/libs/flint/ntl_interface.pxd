# distutils: language = c++
# distutils: libraries = flint
# distutils: depends = flint/NTL-interface.h

from sage.libs.flint.types cimport fmpz_t, fmpz_poly_t

from sage.libs.ntl.ZZ cimport ZZ_c
from sage.libs.ntl.ZZX cimport ZZX_c

# flint/NTL-interface.h
cdef extern from "flint_ntl_wrap.h":
    void fmpz_poly_get_ZZX(ZZX_c output, fmpz_poly_t poly)
    void fmpz_poly_set_ZZX(fmpz_poly_t output, ZZX_c poly)

    void fmpz_get_ZZ(ZZ_c output, fmpz_t f)
    void fmpz_set_ZZ(fmpz_t output, ZZ_c z)
