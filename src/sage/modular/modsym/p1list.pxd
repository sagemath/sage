# sage_setup: distribution = sagemath-schemes


cdef class export:
    cdef int c_p1_normalize_int(self, int N, int u, int v,
                                int* uu, int* vv, int* ss,
                                int compute_s) except -1

    cdef int c_p1_normalize_llong(self, int N, int u, int v,
                                     int* uu, int* vv, int* ss,
                                     int compute_s) except -1


cdef class P1List:
    cdef int __N
    cdef object __list, __end_hash

    cdef int *g
    cdef int *s
    cdef int *t   # xgcd with N table.

    # Here we use a pointer to a function, so the if logic
    # for normalizing an element does not need to be used
    # every time the user calls the normalize function.
    cdef int (*_normalize)(int N, int u, int v,
                            int* uu, int* vv, int* ss,
                            int compute_s) except -1
    cpdef index(self, int u, int v)
    cdef index_and_scalar(self, int u, int v, int* i, int* s)
