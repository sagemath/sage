from libc.stdint cimport uint64_t

cdef int Vrep_list_to_bit_rep(tuple Vrep_list, uint64_t *output,
                                 size_t face_length) except -1

cdef int incidences_to_bit_rep(tuple incidences, uint64_t *output,
                                size_t face_length) except -1

cdef size_t bit_rep_to_Vrep_list(uint64_t *face, size_t *output,
                                    size_t face_length) except -1
