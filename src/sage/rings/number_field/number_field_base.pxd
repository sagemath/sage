from sage.rings.ring cimport Field

cdef class NumberField(Field):
    cdef int _embedded_real
    cdef list _gen_approx

    cpdef _get_embedding_approx(self, size_t i)
