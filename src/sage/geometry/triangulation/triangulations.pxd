# sage_setup: distribution = sagemath-polyhedra
cdef extern from "triangulations.h":
    ctypedef void* triangulations_ptr
    cdef triangulations_ptr init_triangulations \
        (int n, int d, int star, bint fine, object seed, object flips)
    cdef object next_triangulation(triangulations_ptr)
    cdef void delete_triangulations(triangulations_ptr)
