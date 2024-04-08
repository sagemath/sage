# sage_setup: distribution = sagemath-objects
from sage.libs.gmp.types cimport gmp_randstate_t

# The c_random() method on randstate objects gives a value
# 0 <= n <= SAGE_RAND_MAX
cdef extern from *:
    int SAGE_RAND_MAX "(0x7fffffff)"  # 2^31 - 1


cdef class randstate:
    cdef gmp_randstate_t gmp_state
    cdef object _seed
    cdef object _python_random

    cdef object _gap_saved_seed
    cdef object _pari_saved_seed

    cdef object _gp_saved_seeds

    cpdef set_seed_libc(self, bint force)
    cpdef set_seed_ntl(self, bint force)

    cpdef int c_random(self) noexcept
    cpdef double c_rand_double(self) noexcept

    cpdef ZZ_seed(self)
    cpdef long_seed(self)

cpdef randstate current_randstate()
cpdef int random() noexcept
