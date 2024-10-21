# sage_setup: distribution = sagemath-polyhedra
from sage.geometry.polyhedron.combinatorial_polyhedron.face_list_data_structure cimport face_t

cdef int Vrep_list_to_bit_rep(tuple Vrep_list, face_t output) except -1

cdef int incidences_to_bit_rep(tuple incidences, face_t output) except -1

cdef size_t bit_rep_to_Vrep_list(face_t face, size_t *output) except -1
