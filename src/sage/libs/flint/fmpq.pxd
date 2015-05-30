from sage.libs.gmp.types cimport mpq_t
from sage.libs.flint.types cimport fmpq_t

cdef extern from "flint/fmpq.h":
    # Memory management
    void fmpq_init(fmpq_t)
    void fmpq_clear(fmpq_t)

    # Conversion
    void fmpq_set_mpq(fmpq_t, const mpq_t)
    void fmpq_get_mpq(mpq_t, const fmpq_t)
