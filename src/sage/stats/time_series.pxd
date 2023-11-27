cdef class TimeSeries:
    cdef double* _values
    cdef Py_ssize_t _length
    cpdef rescale(self, double s) noexcept
    cpdef double sum(self) noexcept
