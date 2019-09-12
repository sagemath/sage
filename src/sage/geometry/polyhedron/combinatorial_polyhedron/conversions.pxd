from libc.stdint cimport uint64_t

cdef int Vrepr_list_to_bit_repr(tuple Vrepr_list, uint64_t *output,
                                 size_t face_length) except -1

cdef int incidences_to_bit_repr(tuple incidences, uint64_t *output,
                                size_t face_length) except -1

cdef size_t bit_repr_to_Vrepr_list(uint64_t *face, size_t *output,
                                    size_t face_length) except -1
