# sage_setup: distribution = sagemath-linbox
from sage.libs.gmp.types cimport mpz_t

cdef extern from "iml.h":
    cdef enum SOLU_POS:
        LeftSolu
        RightSolu

    cdef long nullspaceMP(long n, long m, const mpz_t *A, mpz_t * *mp_N_pass)
    cdef void nonsingSolvLlhsMM(SOLU_POS solupos, long n, long m, mpz_t *mp_A, mpz_t *mp_B, mpz_t *mp_N, mpz_t mp_D)
