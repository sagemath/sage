# sage_setup: distribution = sagemath-graphs

from libc.stdint cimport uint64_t, uint32_t, INT32_MAX, UINT32_MAX

from sage.data_structures.bitset cimport *
from sage.graphs.base.c_graph cimport CGraph, CGraphBackend
from sage.graphs.base.static_sparse_graph cimport short_digraph, ushort


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
    cdef StaticSparseCGraph _cg
    cdef inline CGraph cg(self):
        return <CGraph> self._cg
    cdef int _use_edge_iterator_on_subgraph(self, CGraphBackend other, object vertices, const int modus) except -1
