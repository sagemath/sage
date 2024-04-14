# sage_setup: distribution = sagemath-graphs
from sage.data_structures.binary_matrix cimport binary_matrix_t
from libc.stdint cimport uint32_t, uint64_t

cdef dict dense_graph_init(binary_matrix_t m, g, translation = ?, force_undirected = ?)
