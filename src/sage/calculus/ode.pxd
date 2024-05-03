cdef class ode_system:
   cdef int  c_j(self,double , double *, double *,double *) noexcept

   cdef int c_f(self,double t, double* , double* ) noexcept
