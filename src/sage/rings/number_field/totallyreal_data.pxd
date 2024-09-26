cdef double eval_seq_as_poly(int *f, int n, double x) noexcept
cdef double newton(int *f, int *df, int n, double x0, double eps) noexcept
cdef void newton_in_intervals(int *f, int *df, int n, double *beta, double eps, double *rts) noexcept
cpdef lagrange_degree_3(int n, int an1, int an2, int an3)

cimport sage.rings.integer

cdef int eval_seq_as_poly_int(int *f, int n, int x) noexcept

cdef int easy_is_irreducible(int *a, int n) noexcept

cdef class tr_data:

    cdef int n, k
    cdef double B
    cdef double b_lower, b_upper, gamma

    cdef int *a
    cdef int *amax
    cdef double *beta
    cdef int *gnk

    cdef int *df

    cdef void incr(self, int *f_out, int verbose, int haltk, int phc) noexcept

