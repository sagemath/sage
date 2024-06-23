# sage_setup: distribution = sagemath-objects
# distutils: language = c++
# distutils: libraries = gmpxx gmp

from sage.libs.gmp.types cimport mpz_t, mpq_t

cdef extern from 'gmpxx.h':
    cdef cppclass mpz_class:
        mpz_class()
        mpz_class(int i)
        mpz_class(mpz_t z)
        mpz_class(mpz_class)
        mpz_t get_mpz_t()
        mpz_class operator%(mpz_class, mpz_class)

    cdef cppclass mpq_class:
        mpq_class()
        mpz_t get_num_mpz_t()
        mpz_t get_den_mpz_t()
