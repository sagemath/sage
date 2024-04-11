from sage.libs.gsl.interp cimport *

cdef class Spline:
    cdef double *x
    cdef double *y
    cdef gsl_interp_accel *acc
    cdef gsl_spline *spline
    cdef int started
    cdef object v

    cdef start_interp(self)
    cdef stop_interp(self)
