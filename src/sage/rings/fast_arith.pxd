cpdef prime_range(start, stop=*, algorithm=*, bint py_ints=*)

cdef class arith_int:
    cdef int abs_int(self, int x) except -1
    cdef int sign_int(self, int n) except -2
    cdef int c_gcd_int(self, int a, int b) except -1
    cdef int c_xgcd_int(self, int a, int b, int* ss, int* tt) except -1
    cdef int c_inverse_mod_int(self, int a, int m) except -1
    cdef int c_rational_recon_int(self, int a, int m, int* n, int* d) except -1

cdef class arith_llong:
    cdef long long abs_longlong(self, long long x) except -1
    cdef long long sign_longlong(self, long long n) except -2
    cdef long long c_gcd_longlong(self, long long a, long long b) except -1
    cdef long long c_xgcd_longlong(self, long long a, long long b,
                                   long long *ss,
                                   long long *tt) except -1
    cdef long long c_inverse_mod_longlong(self, long long a, long long m) except -1
    cdef long long c_rational_recon_longlong(self, long long a, long long m,
                                             long long *n, long long *d) except -1
