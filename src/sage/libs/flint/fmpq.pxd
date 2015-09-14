from sage.libs.gmp.types cimport mpq_t
from sage.libs.flint.types cimport fmpq_t

cdef extern from "flint/fmpq.h":
    # Memory management
    void fmpq_init(fmpq_t)
    void fmpq_clear(fmpq_t)

    # Conversion
    void fmpq_set_mpq(fmpq_t, const mpq_t)
    void fmpq_get_mpq(mpq_t, const fmpq_t)

    # Comparison
    int fmpq_is_zero(const fmpq_t res)
    int fmpq_is_one(const fmpq_t res)
    int fmpq_equal(const fmpq_t x, const fmpq_t y)
    int fmpq_sgn(const fmpq_t)
    int fmpq_cmp(const fmpq_t x, const fmpq_t y)

