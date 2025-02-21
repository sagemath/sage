# sage_setup: distribution = sagemath-objects
# distutils: libraries = gmp

cdef extern from "gmp.h":
    void mp_set_memory_functions(
        void *(*) (size_t),
        void *(*) (void *, size_t, size_t),
        void (*) (void *, size_t))
