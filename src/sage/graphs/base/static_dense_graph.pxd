from sage.data_structures.binary_matrix cimport binary_matrix_t
from libc.stdint cimport uint32_t, uint64_t
from sage.ext.memory cimport check_calloc

cdef dict dense_graph_init(binary_matrix_t m, g, translation = ?)
