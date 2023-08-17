# distutils: libraries = flint
# distutils: depends = flint/qsieve.h

from sage.libs.flint.types cimport *

# flint/qsieve.h
cdef extern from "flint_wrap.h":
    void qsieve_factor(fmpz_factor_t, const fmpz_t)
