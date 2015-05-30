from sage.libs.flint.types cimport fmpz, slong, ulong

cdef extern from "flint/fmpz_vec.h":
    ulong _fmpz_vec_max_limbs(const fmpz *, slong)
