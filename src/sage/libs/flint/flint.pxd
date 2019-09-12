# distutils: libraries = flint
# distutils: depends = flint/flint.h flint/fmpz.h

# flint/flint.h
cdef extern from "flint_wrap.h":
    cdef unsigned long FLINT_BIT_COUNT(unsigned long)
    void flint_free(void * ptr)

# flint/fmpz.h
cdef extern from "flint_wrap.h":
    void _fmpz_cleanup()
    void _fmpz_cleanup_mpz_content()
