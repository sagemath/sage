from .c_graph cimport CGraph, CGraphBackend
from .static_sparse_graph cimport short_digraph, ushort

cdef class StaticSparseCGraph(CGraph):
    cdef short_digraph g
    cdef short_digraph g_rev
    cdef bint _directed
    cdef int * number_of_loops

    cpdef int out_degree(self, int u) except -1
    cpdef int in_degree(self, int u) except -1

cdef class StaticSparseBackend(CGraphBackend):
    cdef int _order
    cdef bint _multiedges
    cdef list _vertex_to_labels
    cdef dict _vertex_to_int
