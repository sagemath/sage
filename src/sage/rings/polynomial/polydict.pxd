cdef class PolyDict:
    cdef object __repn
    cdef object __zero

cdef class ETuple:
    cdef size_t _length
    cdef size_t _nonzero
    cdef int *_data

    cpdef ETuple eadd(ETuple self, ETuple self)
    cpdef ETuple esub(ETuple self, ETuple self)
    cpdef ETuple emul(ETuple self, int factor)
    cpdef ETuple emin(ETuple self, ETuple self)
    cpdef ETuple emax(ETuple self, ETuple self)
    cpdef ETuple eadd_p(ETuple self, int other, int pos)
    cdef ETuple _new(ETuple self)
