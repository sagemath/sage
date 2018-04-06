# distutils: libraries = flint

from sage.libs.flint.types cimport fmpz, slong, ulong, fmpz_t

cdef extern from "flint/fmpz_vec.h":
    fmpz * _fmpz_vec_init(slong)
    void _fmpz_vec_clear(fmpz *, slong)
    ulong _fmpz_vec_max_limbs(const fmpz *, slong)
    void _fmpz_vec_scalar_mod_fmpz(fmpz *, fmpz *, long, fmpz_t)
    void _fmpz_vec_scalar_smod_fmpz(fmpz *, fmpz *, long, fmpz_t)
